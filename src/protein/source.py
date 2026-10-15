"""ENSP -> protein sequence from the Ensembl REST API, with a per-protein disk cache."""

from __future__ import annotations

import logging
from collections.abc import Callable, Iterable
from pathlib import Path

import requests

from src import config
from src.gdc.http import make_session

log = logging.getLogger(__name__)


class ProteinSource:
    def __init__(
        self,
        cache_dir: Path | None = None,
        session: requests.Session | None = None,
        base_url: str | None = None,
    ):
        self.cache_dir = cache_dir or config.PROTEIN_CACHE_DIR
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.session = session or make_session()
        self.base_url = base_url or config.ENSEMBL_BASE_URL
        # IDs whose batch failed after retries (outage, 5xx). Not cached, so a re-run retries them.
        self.failed: set[str] = set()

    def _path(self, ensp: str) -> Path:
        return self.cache_dir / f"{ensp}.txt"

    def cached(self, ensp: str) -> str | None:
        p = self._path(ensp)
        return p.read_text().strip() if p.exists() else None

    def fetch_many(
        self,
        ensps: Iterable[str],
        progress: Callable[[int, int], None] | None = None,
    ) -> dict[str, str]:
        """Return {ensp: sequence}. IDs Ensembl does not know are absent from the result."""
        ids = sorted({e.split(".")[0] for e in ensps if e})
        out: dict[str, str] = {}
        todo = []
        for e in ids:
            seq = self.cached(e)
            if seq:
                out[e] = seq
            else:
                todo.append(e)
        batches = [
            todo[i : i + config.ENSEMBL_BATCH_SIZE]
            for i in range(0, len(todo), config.ENSEMBL_BATCH_SIZE)
        ]
        failures_in_row = 0
        for n, batch in enumerate(batches, start=1):
            if failures_in_row >= config.ENSEMBL_MAX_CONSECUTIVE_FAILURES:
                self.failed.update(batch)  # Ensembl looks down; skip instead of waiting
                fetched = {}
            else:
                try:
                    fetched = self._fetch_batch(batch)
                    failures_in_row = 0
                except requests.RequestException as exc:
                    # Keep partial results: one failed batch should not sink the whole run.
                    log.warning("Ensembl batch of %d failed after retries: %s", len(batch), exc)
                    self.failed.update(batch)
                    failures_in_row += 1
                    fetched = {}
            for ensp, seq in fetched.items():
                self._path(ensp).write_text(seq)
                out[ensp] = seq
            if progress:
                progress(n, len(batches))
        missing = set(ids) - set(out) - self.failed
        if missing:
            log.warning("Ensembl returned no sequence for %d proteins", len(missing))
        return out

    def _fetch_batch(self, ids: list[str]) -> dict[str, str]:
        r = self.session.post(
            f"{self.base_url}/sequence/id",
            params={"type": "protein"},
            json={"ids": ids},
            headers={"Content-Type": "application/json", "Accept": "application/json"},
            timeout=config.HTTP_TIMEOUT_S,
        )
        if r.status_code == 400:
            # One bad ID fails the whole batch; fall back to one-by-one.
            return {k: v for i in ids for k, v in self._fetch_one(i).items()}
        r.raise_for_status()
        return {item["query"].split(".")[0]: item["seq"] for item in r.json() if item.get("seq")}

    def _fetch_one(self, ensp: str) -> dict[str, str]:
        r = self.session.get(
            f"{self.base_url}/sequence/id/{ensp}",
            params={"type": "protein"},
            headers={"Accept": "application/json"},
            timeout=config.HTTP_TIMEOUT_S,
        )
        if r.status_code in (400, 404):
            return {}
        r.raise_for_status()
        return {ensp: r.json()["seq"]}

    def release(self) -> str:
        """Ensembl data release for provenance; 'unknown' if the info endpoint is down."""
        try:
            r = self.session.get(
                f"{self.base_url}/info/data",
                headers={"Accept": "application/json"},
                timeout=config.HTTP_TIMEOUT_S,
            )
            r.raise_for_status()
            return ",".join(str(x) for x in r.json().get("releases", [])) or "unknown"
        except requests.RequestException:
            return "unknown"

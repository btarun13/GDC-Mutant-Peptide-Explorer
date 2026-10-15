"""GDC API calls: projects, open MAF file search, download with md5 check and disk cache."""

from __future__ import annotations

import hashlib
import logging
from collections.abc import Callable, Iterable
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path

import requests

from src import config
from src.gdc.http import make_session

log = logging.getLogger(__name__)

PROJECT_FIELDS = ("project_id", "name", "disease_type", "primary_site", "program.name")
FILE_FIELDS = (
    "file_id",
    "file_name",
    "file_size",
    "md5sum",
    "data_release",
    "analysis.workflow_type",
    "cases.case_id",
    "cases.project.project_id",
)


@dataclass(frozen=True)
class MafFile:
    file_id: str
    file_name: str
    file_size: int
    md5sum: str
    data_release: str
    workflow_type: str
    project_id: str
    case_id: str


def open_maf_filter(project_ids: Iterable[str], workflow_type: str | None = None) -> dict:
    content = [
        {"op": "in", "content": {"field": "cases.project.project_id", "value": list(project_ids)}},
        {"op": "=", "content": {"field": "data_type", "value": config.MAF_DATA_TYPE}},
        {"op": "=", "content": {"field": "data_format", "value": config.MAF_DATA_FORMAT}},
        {"op": "=", "content": {"field": "access", "value": config.MAF_ACCESS}},
    ]
    if workflow_type:
        content.append(
            {"op": "=", "content": {"field": "analysis.workflow_type", "value": workflow_type}}
        )
    return {"op": "and", "content": content}


class GdcClient:
    def __init__(self, session: requests.Session | None = None, base_url: str | None = None):
        self.session = session or make_session()
        self.base_url = base_url or config.GDC_BASE_URL

    def list_projects(self) -> list[dict]:
        hits: list[dict] = []
        start = 0
        while True:
            r = self.session.get(
                f"{self.base_url}/projects",
                params={
                    "size": config.GDC_PAGE_SIZE,
                    "from": start,
                    "fields": ",".join(PROJECT_FIELDS),
                },
                timeout=config.HTTP_TIMEOUT_S,
            )
            r.raise_for_status()
            data = r.json()["data"]
            hits.extend(data["hits"])
            start += len(data["hits"])
            if not data["hits"] or start >= data["pagination"]["total"]:
                return hits

    def workflow_types(self, project_ids: list[str]) -> dict[str, int]:
        """Open MAF counts per analysis.workflow_type, from a facet query."""
        payload = {
            "filters": open_maf_filter(project_ids),
            "facets": "analysis.workflow_type",
            "size": "0",
        }
        r = self.session.post(f"{self.base_url}/files", json=payload, timeout=config.HTTP_TIMEOUT_S)
        r.raise_for_status()
        buckets = r.json()["data"]["aggregations"]["analysis.workflow_type"]["buckets"]
        return {b["key"]: b["doc_count"] for b in buckets}

    def find_open_mafs(
        self, project_ids: list[str], workflow_type: str | None = None
    ) -> list[MafFile]:
        files: list[MafFile] = []
        start = 0
        while True:
            payload = {
                "filters": open_maf_filter(project_ids, workflow_type),
                "fields": ",".join(FILE_FIELDS),
                "format": "JSON",
                "size": str(config.GDC_PAGE_SIZE),
                "from": str(start),
                "sort": "file_name:asc",
            }
            r = self.session.post(
                f"{self.base_url}/files", json=payload, timeout=config.HTTP_TIMEOUT_S
            )
            r.raise_for_status()
            data = r.json()["data"]
            files.extend(_to_maf_file(h) for h in data["hits"])
            start += len(data["hits"])
            if not data["hits"] or start >= data["pagination"]["total"]:
                break
        log.info("Found %d open MAF files for %s", len(files), project_ids)
        return files

    def download_maf(self, maf: MafFile, cache_dir: Path | None = None) -> Path:
        cache_dir = cache_dir or config.MAF_CACHE_DIR
        cache_dir.mkdir(parents=True, exist_ok=True)
        dest = cache_dir / f"{maf.file_id}.maf.gz"
        if dest.exists() and md5_of(dest) == maf.md5sum:
            return dest
        tmp = dest.with_suffix(".part")
        with self.session.get(
            f"{self.base_url}/data/{maf.file_id}",
            stream=True,
            timeout=config.HTTP_DOWNLOAD_TIMEOUT_S,
        ) as r:
            r.raise_for_status()
            with open(tmp, "wb") as fh:
                for chunk in r.iter_content(1 << 20):
                    fh.write(chunk)
        if md5_of(tmp) != maf.md5sum:
            tmp.unlink(missing_ok=True)
            raise ValueError(f"md5 mismatch for GDC file {maf.file_id}")
        tmp.replace(dest)
        return dest

    def download_many(
        self,
        mafs: list[MafFile],
        cache_dir: Path | None = None,
        progress: Callable[[int, int], None] | None = None,
    ) -> tuple[dict[str, Path], dict[str, str]]:
        """Download in parallel. Returns (file_id -> path, file_id -> error message)."""
        paths: dict[str, Path] = {}
        errors: dict[str, str] = {}
        with ThreadPoolExecutor(max_workers=config.GDC_DOWNLOAD_WORKERS) as pool:
            futures = {pool.submit(self.download_maf, m, cache_dir): m for m in mafs}
            for done, fut in enumerate(as_completed(futures), start=1):
                maf = futures[fut]
                try:
                    paths[maf.file_id] = fut.result()
                except (requests.RequestException, ValueError, OSError) as exc:
                    log.warning("Download failed for %s: %s", maf.file_id, exc)
                    errors[maf.file_id] = str(exc)
                if progress:
                    progress(done, len(mafs))
        return paths, errors


def _to_maf_file(hit: dict) -> MafFile:
    cases = hit.get("cases") or [{}]
    return MafFile(
        file_id=hit["file_id"],
        file_name=hit.get("file_name", ""),
        file_size=int(hit.get("file_size") or 0),
        md5sum=hit.get("md5sum", ""),
        data_release=hit.get("data_release", ""),
        workflow_type=(hit.get("analysis") or {}).get("workflow_type", ""),
        project_id=(cases[0].get("project") or {}).get("project_id", ""),
        case_id=cases[0].get("case_id", ""),
    )


def md5_of(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()

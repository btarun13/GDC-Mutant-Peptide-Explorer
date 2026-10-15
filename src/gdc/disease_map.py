"""Free-text disease -> GDC project matches, done locally against the full project list."""

from __future__ import annotations

import difflib
from dataclasses import dataclass, field


@dataclass
class ProjectMatch:
    project_id: str
    name: str
    program: str
    matched_on: list[str] = field(default_factory=list)  # e.g. "disease_type: Nevi and Melanomas"
    fuzzy: bool = False
    score: float = 1.0

    @property
    def label(self) -> str:
        return f"{self.project_id} - {self.name}"


def _searchable_fields(project: dict) -> list[tuple[str, str]]:
    out = [("project_id", project.get("project_id", "")), ("name", project.get("name", ""))]
    for key in ("disease_type", "primary_site"):
        values = project.get(key) or []
        if isinstance(values, str):
            values = [values]
        out.extend((key, v) for v in values if v and v != "Not Reported")
    return out


def match_projects(query: str, projects: list[dict], max_fuzzy: int = 5) -> list[ProjectMatch]:
    """Case-insensitive substring match first; if nothing matches, the closest few by fuzzy
    ratio, marked fuzzy=True so the UI can say no exact match was found."""
    q = query.strip().lower()
    if not q:
        return []
    exact: list[ProjectMatch] = []
    for p in projects:
        hits = [f"{k}: {v}" for k, v in _searchable_fields(p) if q in v.lower()]
        if hits:
            exact.append(_match(p, hits))
    if exact:
        # TCGA first (the default MAF workflow is well characterised there), then by id.
        return sorted(exact, key=lambda m: (m.program != "TCGA", m.project_id))

    scored: list[ProjectMatch] = []
    for p in projects:
        best_score, best_field = 0.0, ""
        for k, v in _searchable_fields(p):
            for token in [v.lower(), *v.lower().split()]:
                s = difflib.SequenceMatcher(None, q, token).ratio()
                if s > best_score:
                    best_score, best_field = s, f"{k}: {v}"
        m = _match(p, [best_field])
        m.fuzzy, m.score = True, best_score
        scored.append(m)
    scored.sort(key=lambda m: -m.score)
    return scored[:max_fuzzy]


def _match(project: dict, hits: list[str]) -> ProjectMatch:
    return ProjectMatch(
        project_id=project["project_id"],
        name=project.get("name", ""),
        program=(project.get("program") or {}).get("name", ""),
        matched_on=hits,
    )

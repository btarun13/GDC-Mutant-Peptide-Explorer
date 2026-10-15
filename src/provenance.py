"""Provenance record attached to every result view and export."""

from __future__ import annotations

import json
import platform
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone

from src import config
from src.gdc.client import MafFile


@dataclass
class Provenance:
    disease_query: str
    project_ids: list[str]
    workflow_type: str
    gdc_files: list[dict]
    genome_build: str = config.GENOME_BUILD
    gdc_gene_model: str = config.GDC_GENE_MODEL
    ensembl_release: str = "unknown"
    predictors: dict[str, str] = field(default_factory=dict)  # name -> version
    parameters: dict = field(default_factory=dict)
    funnel: dict = field(default_factory=dict)
    notes: list[str] = field(
        default_factory=lambda: [config.RESEARCH_DISCLAIMER, config.MASKED_MAF_NOTE]
    )
    timestamp_utc: str = field(
        default_factory=lambda: datetime.now(timezone.utc).isoformat(timespec="seconds")
    )
    python: str = platform.python_version()

    @staticmethod
    def files_from(mafs: list[MafFile]) -> list[dict]:
        return [
            {
                "file_id": m.file_id,
                "file_name": m.file_name,
                "data_release": m.data_release,
                "project_id": m.project_id,
                "case_id": m.case_id,
                "md5sum": m.md5sum,
            }
            for m in mafs
        ]

    def to_json(self) -> str:
        return json.dumps(asdict(self), indent=2, default=str)

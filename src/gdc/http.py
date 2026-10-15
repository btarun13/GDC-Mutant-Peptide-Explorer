"""Shared HTTP session with timeouts and retry/backoff on 429 and 5xx."""

from __future__ import annotations

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

from src import config


def make_session() -> requests.Session:
    retry = Retry(
        total=config.HTTP_RETRIES,
        backoff_factor=config.HTTP_BACKOFF_FACTOR,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=frozenset({"GET", "POST"}),
        respect_retry_after_header=True,
    )
    session = requests.Session()
    adapter = HTTPAdapter(max_retries=retry, pool_maxsize=config.GDC_DOWNLOAD_WORKERS * 2)
    session.mount("https://", adapter)
    session.mount("http://", adapter)
    return session

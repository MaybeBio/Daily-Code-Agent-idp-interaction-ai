import logging
import time
from pathlib import Path

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

logger = logging.getLogger(__name__)


class MobiDBClient:
    """Client for fetching data from the MobiDB API."""

    BASE_URL = "https://mobidb.org/api"

    def __init__(self, cache_dir: Path | None = None, page_size: int = 100):
        self.cache_dir = cache_dir
        self.page_size = page_size
        if self.cache_dir:
            self.cache_dir.mkdir(parents=True, exist_ok=True)
            
        self.session = requests.Session()
        self.session.headers.update({
            "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36 IDPInterp/0.1"
        })
        
        retries = Retry(total=5, backoff_factor=1, status_forcelist=[429, 500, 502, 503, 504])
        self.session.mount("https://", HTTPAdapter(max_retries=retries))

    def fetch_by_filter(
        self,
        feature_filter: str | None = None,
        reviewed_only: bool = False,
        max_length: int | None = None,
        max_entries: int | None = None,
    ) -> list[dict]:
        """
        Fetch proteins matching the given criteria using cursor pagination.
        """
        logger.info(f"MobiDB fetch: filter={feature_filter}, max_len={max_length}, reviewed={reviewed_only}")

        entries = []
        last_id = None
        consecutive_failures = 0
        max_failures = 5

        while True:
            params = {
                "limit": min(self.page_size, 1000),
                "sort_asc": "_id",
                "format": "json"
            }
            if feature_filter:
                params[feature_filter] = "exists"
            if max_length:
                params["length_max"] = max_length
            if last_id:
                params["_id_gt"] = last_id

            try:
                resp = self.session.get(f"{self.BASE_URL}/download_page", params=params, timeout=30)
                resp.raise_for_status()
                data = resp.json()
                consecutive_failures = 0
            except Exception as e:
                consecutive_failures += 1
                logger.error(
                    f"Failed to fetch from MobiDB ({consecutive_failures}/{max_failures}): {e}"
                )
                if consecutive_failures >= max_failures:
                    raise RuntimeError(
                        f"MobiDB fetch failed {max_failures}x in a row "
                        f"(filter={feature_filter}). Giving up."
                    )
                time.sleep(5)
                continue

            batch = data.get("data", [])
            if not batch:
                break

            # Filter locally for reviewed if requested
            if reviewed_only:
                batch = [e for e in batch if e.get("reviewed") is True]

            entries.extend(batch)
            if max_entries and len(entries) >= max_entries:
                entries = entries[:max_entries]
                break

            last_id = data.get("metadata", {}).get("last_id")
            if not last_id:
                break

            if len(entries) % 1000 == 0:
                logger.info(f"Fetched {len(entries)} entries...")

        logger.info(f"Total fetched: {len(entries)}")
        return entries

import os
from pathlib import Path

# Static configuration
STATIC_DATA_DIR = Path(__file__).parent.parent.parent / "static" / "data"
STATIC_BASE_URL = "/static/data"

def get_local_storage_url(filename: str) -> str:
    """
    Returns the local URL path for a given structural file.
    Assumes the file will be placed in backend/static/data/.
    """
    return f"{STATIC_BASE_URL}/{filename}"

def get_local_filepath(filename: str) -> Path:
    """
    Returns the absolute path on disk for a given structural file.
    """
    return STATIC_DATA_DIR / filename

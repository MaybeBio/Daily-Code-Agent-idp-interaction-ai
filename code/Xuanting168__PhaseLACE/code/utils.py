from __future__ import annotations

import hashlib
import json
import os
import threading
import time
from datetime import datetime, timezone
from pathlib import Path


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def atomic_json(path: Path, value: dict) -> None:
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(value, handle, ensure_ascii=False, indent=2, allow_nan=False)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


class Status:
    def __init__(self, path: Path, phase: str) -> None:
        self.path = path
        self.started = time.monotonic()
        self.data = {"status": "RUNNING", "phase": phase, "pid": os.getpid(), "started_at": utc_now()}
        self.lock = threading.Lock()
        self.stop = threading.Event()
        self.update()
        self.thread = threading.Thread(target=self._heartbeat, daemon=True)
        self.thread.start()

    def update(self, **fields: object) -> None:
        with self.lock:
            self.data.update(fields)
            self.data["heartbeat_at"] = utc_now()
            self.data["elapsed_seconds"] = round(time.monotonic() - self.started, 2)
            atomic_json(self.path, self.data)

    def _heartbeat(self) -> None:
        while not self.stop.wait(20):
            self.update()

    def close(self, status: str, **fields: object) -> None:
        self.stop.set()
        self.thread.join(timeout=3)
        self.update(status=status, **fields)


LOCK_EX=1
LOCK_NB=2
LOCK_UN=4


def flock(handle, operation: int) -> None:
    if os.name=='nt':
        import msvcrt
        handle.seek(0,2)
        if handle.tell()==0:
            handle.write(' ');handle.flush()
        handle.seek(0)
        msvcrt.locking(handle.fileno(),msvcrt.LK_UNLCK if operation==LOCK_UN else msvcrt.LK_NBLCK,1)
    else:
        import fcntl
        fcntl.flock(handle,fcntl.LOCK_UN if operation==LOCK_UN else fcntl.LOCK_EX|fcntl.LOCK_NB)


#!/usr/bin/env python3

from __future__ import annotations

import os
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

os.environ.setdefault("CUDA_DEVICE", "3")

from TDP43.tdp43_condensate_simulation.condensate_workflows import main_run_condensate


if __name__ == "__main__":
    main_run_condensate(default_temp_k=260.0)

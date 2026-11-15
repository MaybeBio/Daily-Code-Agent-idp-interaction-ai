from __future__ import annotations

import logging
from pathlib import Path
from datetime import datetime

class Logger:
    def __init__(self, log_dir: str | Path):
        log_dir = Path(log_dir)
        log_dir.mkdir(parents=True, exist_ok=True)

        log_file = log_dir / "log.txt"

        self.logger = logging.getLogger(str(log_file))
        self.logger.setLevel(logging.INFO)

        # 避免重复添加 handler（很关键）
        if not self.logger.handlers:
            handler = logging.FileHandler(log_file, encoding="utf-8")
            formatter = logging.Formatter(
                "[%(asctime)s] %(message)s",
                datefmt="%Y-%m-%d %H:%M:%S"
            )
            handler.setFormatter(formatter)
            self.logger.addHandler(handler)

    def write(self, msg: str) -> None:
        self.logger.info(msg)

    def print(self, msg: str) -> None:
        self.write(msg)  # 写文件
        print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {msg}")  # 控制台输出
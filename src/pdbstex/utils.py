from __future__ import annotations

import hashlib
import os
import re
import sys
import time
from dataclasses import dataclass
from typing import Any, Iterable


@dataclass
class WarningItem:
    level: str
    message: str
    context: dict[str, Any]


class WarningCollector:
    def __init__(self) -> None:
        self._items: list[WarningItem] = []

    def add(self, message: str, level: str = "warning", **context: Any) -> None:
        ctx: dict[str, Any] = {}
        for k, v in context.items():
            ctx[k] = v
        self._items.append(WarningItem(level=level, message=message, context=ctx))

    def items(self) -> list[dict[str, Any]]:
        out: list[dict[str, Any]] = []
        for it in self._items:
            out.append({"level": it.level, "message": it.message, "context": it.context})
        return out

    def merge(self, other: "WarningCollector") -> None:
        for it in other._items:
            self._items.append(it)

    def has_warnings(self) -> bool:
        if len(self._items) > 0:
            return True
        return False


def sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(1024 * 1024)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def is_tty() -> bool:
    try:
        if sys.stdin is not None and sys.stdin.isatty():
            return True
    except Exception:
        return False
    return False


def normalize_pdb_id(value: str) -> str:
    v = value.strip()
    v = v.upper()
    if len(v) == 4 and re.match(r"^[0-9A-Z]{4}$", v):
        return v
    return ""


def ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def safe_relpath(path: str, start: str) -> str:
    try:
        return os.path.relpath(path, start=start)
    except Exception:
        return path


def sleep_backoff(attempt: int, base: float = 0.75, cap: float = 20.0) -> None:
    if attempt < 0:
        attempt = 0
    wait = base * (2.0 ** attempt)
    if wait > cap:
        wait = cap
    time.sleep(wait)

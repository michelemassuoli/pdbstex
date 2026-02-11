from __future__ import annotations

import gzip
import os
from typing import Optional

from .utils import ensure_dir


def write_bytes(path: str, data: bytes) -> None:
    parent = os.path.dirname(path)
    if parent != "":
        ensure_dir(parent)
    with open(path, "wb") as f:
        f.write(data)


def gunzip_bytes(data: bytes) -> Optional[bytes]:
    try:
        return gzip.decompress(data)
    except Exception:
        return None


def gunzip_file(src_path: str, dst_path: str) -> bool:
    try:
        with gzip.open(src_path, "rb") as fin:
            raw = fin.read()
        parent = os.path.dirname(dst_path)
        if parent != "":
            ensure_dir(parent)
        with open(dst_path, "wb") as fout:
            fout.write(raw)
        return True
    except Exception:
        return False

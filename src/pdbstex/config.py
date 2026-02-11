from __future__ import annotations

import json
import os
from dataclasses import dataclass, field
from typing import Any


@dataclass
class Settings:
    out_dir: str = "PDBStEx_output"
    max_entries: int = 0
    include_no_ligand: bool = False
    methods: list[str] = field(default_factory=lambda: ["X-RAY DIFFRACTION"])
    max_resolution: float = 2.5
    include_missing_resolution: bool = False
    bridge_cutoff: float = 3.5
    include_adjacent_polymer_bridge: bool = False
    remove_ions: bool = False
    altloc_mode: str = "split"
    keep_hydrogens: bool = False
    max_rps: float = 2.0
    interactive_rate_limit: bool = True
    cofactor_ids: list[str] = field(default_factory=list)
    ligand_ids: list[str] = field(default_factory=list)
    save_raw_downloads: bool = True


def _as_bool(x: Any, default: bool) -> bool:
    if isinstance(x, bool):
        return x
    if isinstance(x, int):
        if x == 0:
            return False
        return True
    if isinstance(x, str):
        v = x.strip().lower()
        if v in ["true", "1", "yes", "y", "on"]:
            return True
        if v in ["false", "0", "no", "n", "off"]:
            return False
    return default


def _as_int(x: Any, default: int) -> int:
    try:
        return int(x)
    except Exception:
        return default


def _as_float(x: Any, default: float) -> float:
    try:
        return float(x)
    except Exception:
        return default


def _as_str_list(x: Any) -> list[str]:
    out: list[str] = []
    if x is None:
        return out
    if isinstance(x, str):
        s = x.strip()
        if s != "":
            out.append(s)
        return out
    if isinstance(x, list):
        for item in x:
            if isinstance(item, str):
                s = item.strip()
                if s != "":
                    out.append(s)
            else:
                s2 = str(item).strip()
                if s2 != "":
                    out.append(s2)
        return out
    s3 = str(x).strip()
    if s3 != "":
        out.append(s3)
    return out


def load_settings(path: str) -> Settings:
    ext = os.path.splitext(path)[1].lower()
    data: dict[str, Any] = {}
    if ext == ".json":
        with open(path, "r", encoding="utf-8") as f:
            data = json.load(f)
    elif ext == ".toml":
        import tomllib

        with open(path, "rb") as f:
            data = tomllib.load(f)
    else:
        raise ValueError("Config file must be .json or .toml")

    s = Settings()
    if isinstance(data, dict):
        if "out_dir" in data:
            s.out_dir = str(data["out_dir"])
        if "max_entries" in data:
            s.max_entries = _as_int(data["max_entries"], s.max_entries)
        if "include_no_ligand" in data:
            s.include_no_ligand = _as_bool(data["include_no_ligand"], s.include_no_ligand)
        if "methods" in data:
            s.methods = _as_str_list(data["methods"])
        if "max_resolution" in data:
            s.max_resolution = _as_float(data["max_resolution"], s.max_resolution)
        if "include_missing_resolution" in data:
            s.include_missing_resolution = _as_bool(data["include_missing_resolution"], s.include_missing_resolution)
        if "bridge_cutoff" in data:
            s.bridge_cutoff = _as_float(data["bridge_cutoff"], s.bridge_cutoff)
        if "include_adjacent_polymer_bridge" in data:
            s.include_adjacent_polymer_bridge = _as_bool(data["include_adjacent_polymer_bridge"], s.include_adjacent_polymer_bridge)
        if "remove_ions" in data:
            s.remove_ions = _as_bool(data["remove_ions"], s.remove_ions)
        if "altloc_mode" in data:
            s.altloc_mode = str(data["altloc_mode"]).strip()
        if "keep_hydrogens" in data:
            s.keep_hydrogens = _as_bool(data["keep_hydrogens"], s.keep_hydrogens)
        if "max_rps" in data:
            s.max_rps = _as_float(data["max_rps"], s.max_rps)
        if "interactive_rate_limit" in data:
            s.interactive_rate_limit = _as_bool(data["interactive_rate_limit"], s.interactive_rate_limit)
        if "cofactor_ids" in data:
            s.cofactor_ids = _as_str_list(data["cofactor_ids"])
        if "ligand_ids" in data:
            s.ligand_ids = _as_str_list(data["ligand_ids"])
        if "save_raw_downloads" in data:
            s.save_raw_downloads = _as_bool(data["save_raw_downloads"], s.save_raw_downloads)
    return s


def apply_overrides(settings: Settings, overrides: dict[str, Any]) -> Settings:
    s = Settings(**settings.__dict__)
    if "out_dir" in overrides and overrides["out_dir"] is not None:
        s.out_dir = str(overrides["out_dir"])
    if "max_entries" in overrides and overrides["max_entries"] is not None:
        s.max_entries = int(overrides["max_entries"])
    if "include_no_ligand" in overrides and overrides["include_no_ligand"] is not None:
        s.include_no_ligand = bool(overrides["include_no_ligand"])
    if "methods" in overrides and overrides["methods"] is not None:
        s.methods = []
        for m in overrides["methods"]:
            s.methods.append(str(m))
    if "max_resolution" in overrides and overrides["max_resolution"] is not None:
        s.max_resolution = float(overrides["max_resolution"])
    if "include_missing_resolution" in overrides and overrides["include_missing_resolution"] is not None:
        s.include_missing_resolution = bool(overrides["include_missing_resolution"])
    if "bridge_cutoff" in overrides and overrides["bridge_cutoff"] is not None:
        s.bridge_cutoff = float(overrides["bridge_cutoff"])
    if "include_adjacent_polymer_bridge" in overrides and overrides["include_adjacent_polymer_bridge"] is not None:
        s.include_adjacent_polymer_bridge = bool(overrides["include_adjacent_polymer_bridge"])
    if "remove_ions" in overrides and overrides["remove_ions"] is not None:
        s.remove_ions = bool(overrides["remove_ions"])
    if "altloc_mode" in overrides and overrides["altloc_mode"] is not None:
        s.altloc_mode = str(overrides["altloc_mode"])
    if "keep_hydrogens" in overrides and overrides["keep_hydrogens"] is not None:
        s.keep_hydrogens = bool(overrides["keep_hydrogens"])
    if "max_rps" in overrides and overrides["max_rps"] is not None:
        s.max_rps = float(overrides["max_rps"])
    if "interactive_rate_limit" in overrides and overrides["interactive_rate_limit"] is not None:
        s.interactive_rate_limit = bool(overrides["interactive_rate_limit"])
    if "cofactor_ids" in overrides and overrides["cofactor_ids"] is not None:
        s.cofactor_ids = []
        for x in overrides["cofactor_ids"]:
            s.cofactor_ids.append(str(x).upper())
    if "ligand_ids" in overrides and overrides["ligand_ids"] is not None:
        s.ligand_ids = []
        for x in overrides["ligand_ids"]:
            s.ligand_ids.append(str(x).upper())
    if "save_raw_downloads" in overrides and overrides["save_raw_downloads"] is not None:
        s.save_raw_downloads = bool(overrides["save_raw_downloads"])
    return s

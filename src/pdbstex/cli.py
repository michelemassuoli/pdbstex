from __future__ import annotations

import argparse
import sys

from .config import load_settings, apply_overrides, Settings
from .pipeline import run_pipeline
from .utils import normalize_pdb_id


def _parse_args(argv: list[str]) -> argparse.Namespace:
    p = argparse.ArgumentParser(prog="pdbstex", description="PDBStEx: RCSB PDB search, download, and split.")
    p.add_argument("terms", nargs="*", help="Search terms (full-text) or PDB IDs.")
    p.add_argument("--config", dest="config", default=None, help="Config file (.json or .toml). CLI flags override config.")
    p.add_argument("--out-dir", dest="out_dir", default=None, help="Output directory.")
    p.add_argument("--pdb-id", dest="pdb_ids", action="append", default=None, help="Explicit PDB id to process (repeatable).")
    p.add_argument("--include-no-ligand", dest="include_no_ligand", action="store_true", default=None, help="Include entries without ligands (nonpolymer entities).")
    p.add_argument("--method", dest="methods", action="append", default=None, help="Allowed experimental method (repeatable). Default: X-RAY DIFFRACTION.")
    p.add_argument("--max-resolution", dest="max_resolution", type=float, default=None, help="Maximum resolution (Angstrom). Default: 2.5.")
    p.add_argument("--include-missing-resolution", dest="include_missing_resolution", action="store_true", default=None, help="Keep entries even if resolution is missing.")
    p.add_argument("--max-entries", dest="max_entries", type=int, default=None, help="Max entries to process (0 means no limit).")
    p.add_argument("--bridge-cutoff", dest="bridge_cutoff", type=float, default=None, help="Distance cutoff for bridge water detection (heavy atoms). Default: 3.5.")
    p.add_argument("--include-adjacent-polymer-bridge", dest="include_adjacent_polymer_bridge", action="store_true", default=None, help="Keep bridge waters between adjacent polymer residues.")
    p.add_argument("--remove-ions", dest="remove_ions", action="store_true", default=None, help="Remove ions from generated structures.")
    p.add_argument("--altloc-mode", dest="altloc_mode", default=None, choices=["split", "keep", "best", "first"], help="Altloc handling: split, keep, best, first.")
    p.add_argument("--keep-hydrogens", dest="keep_hydrogens", action="store_true", default=None, help="Include hydrogens (default excludes).")
    p.add_argument("--max-rps", dest="max_rps", type=float, default=None, help="Max requests per second to RCSB.")
    p.add_argument("--no-interactive-rate-limit", dest="interactive_rate_limit", action="store_false", default=None, help="Disable interactive prompts on rate limit.")
    p.add_argument("--ligand-id", dest="ligand_ids", action="append", default=None, help="Ligand CCD id (repeatable). If set, cofactors are nonpolymer excluding these ligands by default.")
    p.add_argument("--cofactor-id", dest="cofactor_ids", action="append", default=None, help="Cofactor CCD id (repeatable).")
    p.add_argument("--no-raw", dest="save_raw_downloads", action="store_false", default=None, help="Do not save raw downloads and raw PDB conversion.")
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    if argv is None:
        argv = sys.argv[1:]
    ns = _parse_args(argv)

    settings = Settings()
    if ns.config is not None:
        settings = load_settings(ns.config)

    overrides = {
        "out_dir": ns.out_dir,
        "max_entries": ns.max_entries,
        "include_no_ligand": ns.include_no_ligand,
        "methods": ns.methods,
        "max_resolution": ns.max_resolution,
        "include_missing_resolution": ns.include_missing_resolution,
        "bridge_cutoff": ns.bridge_cutoff,
        "include_adjacent_polymer_bridge": ns.include_adjacent_polymer_bridge,
        "remove_ions": ns.remove_ions,
        "altloc_mode": ns.altloc_mode,
        "keep_hydrogens": ns.keep_hydrogens,
        "max_rps": ns.max_rps,
        "interactive_rate_limit": ns.interactive_rate_limit,
        "cofactor_ids": ns.cofactor_ids,
        "ligand_ids": ns.ligand_ids,
        "save_raw_downloads": ns.save_raw_downloads,
    }

    settings = apply_overrides(settings, overrides)

    pdb_ids: list[str] = []
    if ns.pdb_ids is not None:
        for x in ns.pdb_ids:
            pid = normalize_pdb_id(str(x))
            if pid != "":
                pdb_ids.append(pid)

    terms: list[str] = []
    for t in ns.terms:
        pid = normalize_pdb_id(str(t))
        if pid != "":
            pdb_ids.append(pid)
        else:
            terms.append(str(t))

    if len(terms) == 0 and len(pdb_ids) == 0:
        print("PDBStEx error: provide at least one search term or --pdb-id.", file=sys.stderr)
        raise SystemExit(2)

    manifest = run_pipeline(terms=terms, pdb_ids=pdb_ids, settings=settings)
    if manifest is not None:
        pass

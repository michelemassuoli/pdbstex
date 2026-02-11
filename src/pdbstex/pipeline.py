from __future__ import annotations

import json
import os
from dataclasses import dataclass
from typing import Any, Optional

import gemmi

from .config import Settings
from .io import write_bytes, gunzip_bytes
from .processing import generate_outputs_for_structure, build_metadata
from .rcsb import Requester, search_entries, resolve_obsolete_pdb_id, fetch_entry_metadata
from .utils import WarningCollector, ensure_dir, normalize_pdb_id


@dataclass
class EntryResult:
    pdb_code: str
    status: str
    output_dir: str
    metadata_path: str
    warnings: list[dict[str, Any]]


def _download_and_prepare_structure(requester: Requester, pdb_id: str, out_dir: str, kind: str, warnings: WarningCollector, save_raw: bool) -> Optional[gemmi.Structure]:
    pid = normalize_pdb_id(pdb_id)
    if pid == "":
        return None

    if kind == "asu":
        fname = pid + ".cif.gz"
        raw_base = "asymmetric_unit"
    else:
        fname = pid + "-assembly1.cif.gz"
        raw_base = "biological_assembly_1"

    url = requester.endpoints.files_download + "/" + fname
    data = requester.download_bytes(url, timeout=180.0)
    if data is None:
        warnings.add("Failed to download mmCIF.", pdb_id=pid, kind=kind, url=url)
        return None

    decompressed = gunzip_bytes(data)
    if decompressed is None:
        warnings.add("Failed to decompress mmCIF (.gz).", pdb_id=pid, kind=kind, url=url)
        return None

    raw_dir = os.path.join(out_dir, pid, "raw")
    ensure_dir(raw_dir)

    cif_path = os.path.join(raw_dir, raw_base + ".cif")
    write_bytes(cif_path, decompressed)

    st = None
    try:
        st = gemmi.read_structure(cif_path)
    except Exception as e:
        warnings.add("Failed to parse mmCIF with gemmi.", pdb_id=pid, kind=kind, path=cif_path, error=str(e))
        return None

    if save_raw:
        pdb_path = os.path.join(raw_dir, raw_base + ".pdb")
        try:
            st.write_pdb(pdb_path)
            warnings.add("Converted mmCIF to PDB. Note that PDB format may lose information compared to mmCIF.", pdb_id=pid, kind=kind, path=pdb_path, level="info")
        except Exception as e:
            warnings.add("Failed to write raw PDB conversion.", pdb_id=pid, kind=kind, error=str(e))

    return st


def _passes_entry_filters(entry_meta: Optional[dict[str, Any]], settings: Settings, warnings: WarningCollector, pdb_id: str) -> bool:
    if entry_meta is None:
        warnings.add("Entry metadata not available; cannot apply metadata filters reliably.", pdb_id=pdb_id)
        return True

    methods_ok = True
    if settings.methods is not None and len(settings.methods) > 0:
        methods_ok = False
        mlist: list[str] = []
        try:
            ei = entry_meta.get("rcsb_entry_info", {})
            mval = ei.get("experimental_method", None)
            if isinstance(mval, list):
                for x in mval:
                    mlist.append(str(x))
            elif mval is not None:
                mlist.append(str(mval))
        except Exception:
            mlist = []
        for want in settings.methods:
            for got in mlist:
                if str(got).strip().upper() == str(want).strip().upper():
                    methods_ok = True
                    break
            if methods_ok:
                break
        if methods_ok is False:
            warnings.add("Entry filtered out by experimental method.", pdb_id=pdb_id, methods=mlist, allowed=settings.methods)
            return False

    if settings.max_resolution is not None and settings.max_resolution > 0:
        res = None
        try:
            ei2 = entry_meta.get("rcsb_entry_info", {})
            rval = ei2.get("resolution_combined", None)
            if isinstance(rval, list):
                if len(rval) > 0:
                    try:
                        res = float(rval[0])
                    except Exception:
                        res = None
            elif rval is not None:
                try:
                    res = float(rval)
                except Exception:
                    res = None
        except Exception:
            res = None

        if res is None:
            if settings.include_missing_resolution:
                warnings.add("Resolution missing; keeping entry due to include_missing_resolution.", pdb_id=pdb_id, level="info")
            else:
                warnings.add("Entry filtered out due to missing resolution.", pdb_id=pdb_id)
                return False
        else:
            if float(res) > float(settings.max_resolution):
                warnings.add("Entry filtered out by resolution threshold.", pdb_id=pdb_id, resolution=res, max_resolution=settings.max_resolution)
                return False

    if settings.include_no_ligand is False:
        npl = None
        try:
            ei3 = entry_meta.get("rcsb_entry_info", {})
            nv = ei3.get("nonpolymer_entity_count", None)
            if nv is not None:
                npl = int(nv)
        except Exception:
            npl = None
        if npl is not None:
            if npl <= 0:
                warnings.add("Entry filtered out because it has no nonpolymer entities.", pdb_id=pdb_id)
                return False

    return True


def process_entry(requester: Requester, pdb_id: str, settings: Settings, out_root: str) -> EntryResult:
    warnings = WarningCollector()
    pid0 = normalize_pdb_id(pdb_id)
    pid = resolve_obsolete_pdb_id(requester, pid0, warnings)
    if pid == "":
        return EntryResult(pdb_code=pid0, status="invalid_id", output_dir="", metadata_path="", warnings=warnings.items())

    entry_dir = os.path.join(out_root, pid)
    if os.path.exists(entry_dir):
        warnings.add("Output folder exists; skipping entry to avoid overwrite.", pdb_id=pid, path=entry_dir)
        return EntryResult(pdb_code=pid, status="skipped_exists", output_dir=entry_dir, metadata_path="", warnings=warnings.items())

    ensure_dir(entry_dir)

    entry_meta = fetch_entry_metadata(requester, pid)
    if _passes_entry_filters(entry_meta, settings, warnings, pid) is False:
        return EntryResult(pdb_code=pid, status="filtered", output_dir=entry_dir, metadata_path="", warnings=warnings.items())

    au_st = _download_and_prepare_structure(requester, pid, out_root, "asu", warnings, save_raw=settings.save_raw_downloads)
    ba1_st = _download_and_prepare_structure(requester, pid, out_root, "ba1", warnings, save_raw=settings.save_raw_downloads)

    au_idx = None
    ba1_idx = None

    if au_st is not None:
        au_idx = generate_outputs_for_structure(
            st=au_st,
            out_root=out_root,
            pdb_code=pid,
            dataset="asymmetric_unit",
            ligand_ids=settings.ligand_ids,
            cofactor_ids=settings.cofactor_ids,
            remove_ions=settings.remove_ions,
            altloc_mode=settings.altloc_mode,
            bridge_cutoff=settings.bridge_cutoff,
            include_adjacent_polymer_bridge=settings.include_adjacent_polymer_bridge,
            keep_hydrogens=settings.keep_hydrogens,
            warnings=warnings,
        )

    if ba1_st is not None:
        ba1_idx = generate_outputs_for_structure(
            st=ba1_st,
            out_root=out_root,
            pdb_code=pid,
            dataset="biological_assembly_1",
            ligand_ids=settings.ligand_ids,
            cofactor_ids=settings.cofactor_ids,
            remove_ions=settings.remove_ions,
            altloc_mode=settings.altloc_mode,
            bridge_cutoff=settings.bridge_cutoff,
            include_adjacent_polymer_bridge=settings.include_adjacent_polymer_bridge,
            keep_hydrogens=settings.keep_hydrogens,
            warnings=warnings,
        )

    settings_dict: dict[str, Any] = {}
    for k, v in settings.__dict__.items():
        settings_dict[k] = v

    meta = build_metadata(
        pdb_code=pid,
        entry_meta=entry_meta,
        au_files=au_idx,
        ba1_files=ba1_idx,
        out_root=out_root,
        settings=settings_dict,
        warnings=warnings,
    )

    meta_path = os.path.join(out_root, pid, "metadata.json")
    try:
        with open(meta_path, "w", encoding="utf-8") as f:
            json.dump(meta, f, indent=2, ensure_ascii=False)
    except Exception as e:
        warnings.add("Failed to write metadata.json.", pdb_id=pid, path=meta_path, error=str(e))

    status = "ok"
    if au_st is None and ba1_st is None:
        status = "failed_download"

    return EntryResult(pdb_code=pid, status=status, output_dir=os.path.join(out_root, pid), metadata_path=meta_path, warnings=warnings.items())


def run_pipeline(terms: list[str], pdb_ids: list[str], settings: Settings) -> dict[str, Any]:
    out_root = settings.out_dir
    ensure_dir(out_root)

    warnings = WarningCollector()
    requester = Requester(max_rps=settings.max_rps, interactive=settings.interactive_rate_limit, warnings=warnings)

    all_ids: list[str] = []
    for x in pdb_ids:
        pid = normalize_pdb_id(x)
        if pid != "":
            all_ids.append(pid)
        else:
            warnings.add("Invalid PDB id in --pdb-id.", value=x)

    if len(terms) > 0:
        found = search_entries(
            requester=requester,
            terms=terms,
            methods=settings.methods,
            max_resolution=settings.max_resolution,
            include_no_ligand=settings.include_no_ligand,
            warnings=warnings,
            max_entries=settings.max_entries,
        )
        for pid in found:
            if pid not in all_ids:
                all_ids.append(pid)

    if settings.max_entries > 0:
        limited: list[str] = []
        i = 0
        for pid in all_ids:
            limited.append(pid)
            i += 1
            if i >= settings.max_entries:
                break
        all_ids = limited

    results: list[dict[str, Any]] = []
    idx = 0
    while idx < len(all_ids):
        pid = all_ids[idx]
        print("PDBStEx:", pid, "(" + str(idx + 1) + "/" + str(len(all_ids)) + ")")
        res = process_entry(requester, pid, settings, out_root)
        results.append(
            {
                "PDB_code": res.pdb_code,
                "status": res.status,
                "output_dir": res.output_dir,
                "metadata_json": res.metadata_path,
                "warnings": res.warnings,
            }
        )
        idx += 1

    manifest: dict[str, Any] = {
        "schema_version": "pdbstex_manifest_v1",
        "terms": terms,
        "pdb_ids": pdb_ids,
        "settings": {k: v for k, v in settings.__dict__.items()},
        "run_warnings": warnings.items(),
        "results": results,
    }

    manifest_path = os.path.join(out_root, "manifest.json")
    try:
        if os.path.exists(manifest_path):
            warnings.add("manifest.json already exists; not overwriting.", path=manifest_path)
        else:
            with open(manifest_path, "w", encoding="utf-8") as f:
                json.dump(manifest, f, indent=2, ensure_ascii=False)
    except Exception as e:
        warnings.add("Failed to write manifest.json.", path=manifest_path, error=str(e))

    return manifest

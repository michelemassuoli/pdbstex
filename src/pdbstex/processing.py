from __future__ import annotations

import json
import os
import re
from dataclasses import dataclass
from typing import Any, Optional, Tuple

import gemmi

from .utils import WarningCollector, ensure_dir, safe_relpath, sha256_file, normalize_pdb_id

ION_IDS_DEFAULT = {
    "LI",
    "NA",
    "K",
    "RB",
    "CS",
    "MG",
    "CA",
    "SR",
    "BA",
    "AL",
    "GA",
    "IN",
    "TL",
    "FE",
    "ZN",
    "CU",
    "MN",
    "CO",
    "NI",
    "CD",
    "HG",
    "CL",
    "BR",
    "I",
    "F",
    "AG",
    "AU",
    "PT",
    "IR",
    "OS",
    "PB",
    "SN",
    "SE",
    "CR",
    "V",
    "MO",
    "W",
}

WATER_RESNAMES = {"HOH", "WAT", "H2O", "DOD"}

COMMON_COFACTORS = {
    "ATP",
    "ADP",
    "AMP",
    "GTP",
    "GDP",
    "GMP",
    "CTP",
    "CDP",
    "CMP",
    "UTP",
    "UDP",
    "UMP",
    "NAD",
    "NADP",
    "FAD",
    "FMN",
    "SAM",
    "SAH",
    "COA",
    "ACP",
    "HEM",
    "HEME",
    "PLP",
    "NAP",
    "NAI",
    "NDP",
    "TPP",
    "THM",
    "PMP",
    "P6G",
}

@dataclass
class StructureFiles:
    cif_path: str
    pdb_path: str

@dataclass
class OutputIndex:
    receptor_variants: dict[str, list[StructureFiles]]
    chain_variants: dict[str, dict[str, list[StructureFiles]]]
    ligands: list[StructureFiles]
    cofactors: list[StructureFiles]
    ions: list[StructureFiles]
    metadata_json: str

def _norm_altloc(v: object) -> str:
    s = str(v)
    if s == "\x00":
        return ""
    s = s.strip()
    if s == "":
        return ""
    if s == ".":
        return ""
    if s == "?":
        return ""
    if s == " ":
        return ""
    return s


def _altloc_letters_in_structure(st: gemmi.Structure, keep_hydrogens: bool) -> list[str]:
    letters_set: set[str] = set()
    model = st[0]
    for chain in model:
        for res in chain:
            for atom in res:
                if keep_hydrogens is False:
                    if atom.element.name == "H" or atom.element.name == "D":
                        continue
                al = _norm_altloc(atom.altloc)
                if al != "":
                    letters_set.add(al)
    out: list[str] = []
    for x in sorted(list(letters_set)):
        out.append(x)
    return out

def _choose_best_altloc(st: gemmi.Structure, keep_hydrogens: bool) -> str:
    score: dict[str, float] = {}
    model = st[0]
    for chain in model:
        for res in chain:
            for atom in res:
                if keep_hydrogens is False:
                    if atom.element.name == "H" or atom.element.name == "D":
                        continue
                al = _norm_altloc(atom.altloc)
                if al != "":
                    if al not in score:
                        score[al] = 0.0
                    try:
                        score[al] = score[al] + float(atom.occ)
                    except Exception:
                        score[al] = score[al] + 1.0
    best = ""
    best_val = -1.0
    for k, v in score.items():
        if v > best_val:
            best_val = v
            best = k
    return best

def _make_like_structure(st: gemmi.Structure) -> gemmi.Structure:
    new = gemmi.Structure()
    try:
        new.name = st.name
    except Exception:
        pass
    try:
        new.cell = st.cell
    except Exception:
        pass
    try:
        new.spacegroup_hm = st.spacegroup_hm
    except Exception:
        pass
    try:
        new.resolution = st.resolution
    except Exception:
        pass
    return new

def _clone_residue_filtered(res: gemmi.Residue, altloc_keep: Optional[str], keep_hydrogens: bool) -> Optional[gemmi.Residue]:
    new_res = gemmi.Residue()
    new_res.name = res.name
    new_res.seqid = res.seqid
    try:
        new_res.subchain = res.subchain
    except Exception:
        pass
    try:
        new_res.label_seq = res.label_seq
    except Exception:
        pass
    try:
        new_res.het_flag = res.het_flag
    except Exception:
        pass
    kept_any = False
    for atom in res:
        if keep_hydrogens is False:
            if atom.element.name == "H" or atom.element.name == "D":
                continue
        if altloc_keep is not None and altloc_keep != "":
            al = _norm_altloc(atom.altloc)
            if al != "":
                if al != altloc_keep:
                    continue

        new_res.add_atom(atom.clone())
        kept_any = True
    if kept_any is False:
        return None
    return new_res

def _add_chain_filtered(new_model: gemmi.Model, chain: gemmi.Chain, residues_to_keep: set[Tuple[str, int, str, str]], altloc_keep: Optional[str], keep_hydrogens: bool) -> None:
    new_chain = gemmi.Chain(chain.name)
    added = False
    for res in chain:
        rid = (chain.name, int(res.seqid.num), str(res.seqid.icode), str(res.name))
        if rid in residues_to_keep:
            new_res = _clone_residue_filtered(res, altloc_keep, keep_hydrogens)
            if new_res is not None:
                new_chain.add_residue(new_res)
                added = True
    if added:
        new_model.add_chain(new_chain)

def _res_id_tuple(chain_name: str, res: gemmi.Residue) -> Tuple[str, int, str, str]:
    return (chain_name, int(res.seqid.num), str(res.seqid.icode), str(res.name))

def _classify_residues(st: gemmi.Structure, ligand_ids: list[str], cofactor_ids: list[str], remove_ions: bool, warnings: WarningCollector) -> tuple[set[Tuple[str, int, str, str]], set[Tuple[str, int, str, str]], set[Tuple[str, int, str, str]], set[Tuple[str, int, str, str]]]:
    polymer: set[Tuple[str, int, str, str]] = set()
    ligands: set[Tuple[str, int, str, str]] = set()
    cofactors: set[Tuple[str, int, str, str]] = set()
    waters: set[Tuple[str, int, str, str]] = set()
    lig_set: set[str] = set()
    for x in ligand_ids:
        lig_set.add(str(x).upper())
    cof_set: set[str] = set()
    for x in cofactor_ids:
        cof_set.add(str(x).upper())
    st.setup_entities()
    model = st[0]
    for chain in model:
        for res in chain:
            rid = _res_id_tuple(chain.name, res)
            resname = str(res.name).upper()
            et = None
            try:
                et = res.entity_type
            except Exception:
                et = None
            if resname in WATER_RESNAMES:
                waters.add(rid)
                continue
            if et == gemmi.EntityType.Polymer:
                polymer.add(rid)
                continue
            if remove_ions is True:
                if resname in ION_IDS_DEFAULT:
                    continue
            if len(lig_set) > 0:
                if resname in lig_set:
                    ligands.add(rid)
                else:
                    cofactors.add(rid)
                continue
            if resname in cof_set:
                cofactors.add(rid)
                continue
            if resname in COMMON_COFACTORS:
                cofactors.add(rid)
                continue
            ligands.add(rid)
    if len(ligand_ids) == 0 and len(cofactor_ids) == 0:
        warnings.add("Ligand IDs not provided; cofactors will be inferred with a small built-in list and may be incomplete.")
    if len(ligand_ids) == 0 and len(cofactor_ids) > 0:
        warnings.add("Ligand IDs not provided; ligands are defined as nonpolymer entities excluding waters/ions and listed cofactors.")
    return polymer, ligands, cofactors, waters

def _compute_bridge_waters(st: gemmi.Structure, polymer: set[Tuple[str, int, str, str]], ligands: set[Tuple[str, int, str, str]], cofactors: set[Tuple[str, int, str, str]], waters: set[Tuple[str, int, str, str]], cutoff: float, include_adjacent_polymer_bridge: bool, keep_hydrogens: bool, warnings: WarningCollector) -> set[Tuple[str, int, str, str]]:
    bridge: set[Tuple[str, int, str, str]] = set()
    model = st[0]
    ns = gemmi.NeighborSearch(model, st.cell, cutoff)
    ns.populate(include_h=keep_hydrogens)
    water_oxygen_names = {"O", "OW", "O1", "O2"}
    for chain in model:
        for res in chain:
            rid = _res_id_tuple(chain.name, res)
            if rid not in waters:
                continue
            oxy_atoms: list[gemmi.Atom] = []
            for atom in res:
                if keep_hydrogens is False:
                    if atom.element.name == "H" or atom.element.name == "D":
                        continue
                if atom.element.name == "O":
                    oxy_atoms.append(atom)
                else:
                    an = str(atom.name).strip().upper()
                    if an in water_oxygen_names:
                        oxy_atoms.append(atom)
            if len(oxy_atoms) == 0:
                continue
            partners: set[Tuple[str, int, str, str]] = set()
            partner_polymer: list[Tuple[str, int, str, str]] = []
            partner_nonpolymer: list[Tuple[str, int, str, str]] = []
            for oxy in oxy_atoms:
                marks = ns.find_atoms(oxy.pos, radius=cutoff)
                for m in marks:
                    cra = m.to_cra(model)
                    if cra is None:
                        continue
                    ch_name = cra.chain.name
                    r2 = cra.residue
                    if r2 is None:
                        continue
                    rid2 = _res_id_tuple(ch_name, r2)
                    if rid2 == rid:
                        continue
                    if keep_hydrogens is False:
                        if cra.atom.element.name == "H" or cra.atom.element.name == "D":
                            continue
                    if rid2 in waters:
                        continue
                    if rid2 in polymer:
                        partners.add(rid2)
                        partner_polymer.append(rid2)
                        continue
                    if rid2 in ligands or rid2 in cofactors:
                        partners.add(rid2)
                        partner_nonpolymer.append(rid2)
                        continue
            if len(partners) < 2:
                continue
            if include_adjacent_polymer_bridge is False:
                if len(partner_nonpolymer) == 0 and len(partner_polymer) >= 2:
                    exclude = False
                    i = 0
                    while i < len(partner_polymer):
                        j = i + 1
                        while j < len(partner_polymer):
                            p1 = partner_polymer[i]
                            p2 = partner_polymer[j]
                            if p1[0] == p2[0]:
                                try:
                                    n1 = int(p1[1])
                                    n2 = int(p2[1])
                                except Exception:
                                    n1 = 0
                                    n2 = 0
                                if abs(n1 - n2) <= 1:
                                    exclude = True
                                    break
                            j += 1
                        if exclude:
                            break
                        i += 1
                    if exclude:
                        continue
            bridge.add(rid)
    return bridge

def _write_structure_files(st: gemmi.Structure, cif_path: str, pdb_path: str, warnings: WarningCollector) -> None:
    try:
        doc = st.make_mmcif_document()
        ensure_dir(os.path.dirname(cif_path))
        doc.write_file(cif_path)
    except Exception as e:
        warnings.add("Failed to write mmCIF.", path=cif_path, error=str(e))
    try:
        ensure_dir(os.path.dirname(pdb_path))
        st.write_pdb(pdb_path)
    except Exception as e:
        warnings.add("Failed to write PDB.", path=pdb_path, error=str(e))

def _extract_structure_variant(st: gemmi.Structure, residues_to_keep: set[Tuple[str, int, str, str]], altloc_keep: Optional[str], keep_hydrogens: bool) -> gemmi.Structure:
    new_st = _make_like_structure(st)
    model = gemmi.Model("1")
    src_model = st[0]
    for chain in src_model:
        _add_chain_filtered(model, chain, residues_to_keep, altloc_keep, keep_hydrogens)
    new_st.add_model(model)
    new_st.setup_entities()
    return new_st

def _extract_residue_as_structure(st: gemmi.Structure, rid: Tuple[str, int, str, str], altloc_keep: Optional[str], keep_hydrogens: bool) -> Optional[gemmi.Structure]:
    chain_name = rid[0]
    model = st[0]
    res_found = None
    chain_found = None
    for chain in model:
        if chain.name == chain_name:
            chain_found = chain
            break
    if chain_found is None:
        return None
    for res in chain_found:
        if int(res.seqid.num) == int(rid[1]) and str(res.seqid.icode) == str(rid[2]) and str(res.name) == str(rid[3]):
            res_found = res
            break
    if res_found is None:
        return None
    new_st = _make_like_structure(st)
    new_model = gemmi.Model("1")
    new_chain = gemmi.Chain(chain_name)
    new_res = _clone_residue_filtered(res_found, altloc_keep, keep_hydrogens)
    if new_res is None:
        return None
    new_chain.add_residue(new_res)
    new_model.add_chain(new_chain)
    new_st.add_model(new_model)
    new_st.setup_entities()
    return new_st

def _altloc_plan(st: gemmi.Structure, altloc_mode: str, keep_hydrogens: bool, warnings: WarningCollector) -> list[Optional[str]]:
    mode = str(altloc_mode).strip().lower()
    letters = _altloc_letters_in_structure(st, keep_hydrogens=keep_hydrogens)
    if mode == "keep":
        return [None]
    if mode == "split":
        if len(letters) == 0:
            return [None]
        out: list[Optional[str]] = []
        for x in letters:
            out.append(x)
        return out
    if mode == "best":
        best = _choose_best_altloc(st, keep_hydrogens=keep_hydrogens)
        if best == "":
            return [None]
        return [best]
    if mode == "first":
        if len(letters) == 0:
            return [None]
        if "A" in letters:
            return ["A"]
        return [letters[0]]
    warnings.add("Unknown altloc_mode; using keep.", altloc_mode=altloc_mode)
    return [None]

def _sanitize_filename(s: str) -> str:
    out = ""
    for ch in s:
        if re.match(r"[A-Za-z0-9._-]", ch):
            out += ch
        else:
            out += "_"
    if out == "":
        out = "x"
    return out

def _residue_file_tag(rid: Tuple[str, int, str, str]) -> str:
    chain = _sanitize_filename(rid[0])
    num = str(rid[1])
    icode = _sanitize_filename(rid[2])
    name = _sanitize_filename(rid[3])
    if icode == "" or icode == " ":
        return chain + "_" + num + "_" + name
    return chain + "_" + num + "_" + icode + "_" + name

def generate_outputs_for_structure(st: gemmi.Structure, out_root: str, pdb_code: str, dataset: str, ligand_ids: list[str], cofactor_ids: list[str], remove_ions: bool, altloc_mode: str, bridge_cutoff: float, include_adjacent_polymer_bridge: bool, keep_hydrogens: bool, warnings: WarningCollector) -> OutputIndex:
    ensure_dir(out_root)
    ensure_dir(os.path.join(out_root, pdb_code))
    entry_root = os.path.join(out_root, pdb_code, dataset)
    ensure_dir(entry_root)
    polymer, ligands, cofactors, waters = _classify_residues(st, ligand_ids, cofactor_ids, remove_ions, warnings)
    if len(ligands) == 0:
        warnings.add("No ligands detected in structure.", label=pdb_code)
    bridge_waters = _compute_bridge_waters(st=st, polymer=polymer, ligands=ligands, cofactors=cofactors, waters=waters, cutoff=float(bridge_cutoff), include_adjacent_polymer_bridge=bool(include_adjacent_polymer_bridge), keep_hydrogens=bool(keep_hydrogens), warnings=warnings)
    altloc_choices = _altloc_plan(st, altloc_mode, keep_hydrogens=keep_hydrogens, warnings=warnings)
    receptor_variants: dict[str, list[StructureFiles]] = {}
    chain_variants: dict[str, dict[str, list[StructureFiles]]] = {}
    lig_files: list[StructureFiles] = []
    cof_files: list[StructureFiles] = []
    ion_files: list[StructureFiles] = []
    ions_set: set[Tuple[str, int, str, str]] = set()
    if remove_ions is False:
        model2 = st[0]
        for chain2 in model2:
            for res2 in chain2:
                rid2 = _res_id_tuple(chain2.name, res2)
                resname2 = str(res2.name).upper()
                if resname2 in ION_IDS_DEFAULT:
                    ions_set.add(rid2)
    alt_label_list: list[str] = []
    for alt in altloc_choices:
        if alt is None:
            alt_label_list.append("altNONE")
        else:
            alt_label_list.append("alt" + str(alt))
    idx_alt = 0
    while idx_alt < len(altloc_choices):
        alt = altloc_choices[idx_alt]
        alt_tag = alt_label_list[idx_alt]
        receptor_only = _extract_structure_variant(st, polymer, altloc_keep=alt, keep_hydrogens=keep_hydrogens)
        receptor_plus = set()
        for rid3 in polymer:
            receptor_plus.add(rid3)
        for rid3 in cofactors:
            receptor_plus.add(rid3)
        for rid3 in ions_set:
            receptor_plus.add(rid3)
        receptor_plus_struct = _extract_structure_variant(st, receptor_plus, altloc_keep=alt, keep_hydrogens=keep_hydrogens)
        receptor_bridge = set()
        for rid3 in receptor_plus:
            receptor_bridge.add(rid3)
        for rid3 in bridge_waters:
            receptor_bridge.add(rid3)
        receptor_bridge_struct = _extract_structure_variant(st, receptor_bridge, altloc_keep=alt, keep_hydrogens=keep_hydrogens)
        key_only = pdb_code + "_" + dataset + "_receptor_only_" + alt_tag
        key_plus = pdb_code + "_" + dataset + "_receptor_plus_cofactors_" + alt_tag
        key_bridge = pdb_code + "_" + dataset + "_receptor_bridge_waters_" + alt_tag
        out_dir = os.path.join(out_root, pdb_code, dataset, "receptor")
        ensure_dir(out_dir)
        files_only = StructureFiles(cif_path=os.path.join(out_dir, key_only + ".cif"), pdb_path=os.path.join(out_dir, key_only + ".pdb"))
        files_plus = StructureFiles(cif_path=os.path.join(out_dir, key_plus + ".cif"), pdb_path=os.path.join(out_dir, key_plus + ".pdb"))
        files_bridge = StructureFiles(cif_path=os.path.join(out_dir, key_bridge + ".cif"), pdb_path=os.path.join(out_dir, key_bridge + ".pdb"))
        _write_structure_files(receptor_only, files_only.cif_path, files_only.pdb_path, warnings)
        _write_structure_files(receptor_plus_struct, files_plus.cif_path, files_plus.pdb_path, warnings)
        _write_structure_files(receptor_bridge_struct, files_bridge.cif_path, files_bridge.pdb_path, warnings)
        if "receptor_only" not in receptor_variants:
            receptor_variants["receptor_only"] = []
        receptor_variants["receptor_only"].append(files_only)
        if "receptor_plus_cofactors" not in receptor_variants:
            receptor_variants["receptor_plus_cofactors"] = []
        receptor_variants["receptor_plus_cofactors"].append(files_plus)
        if "receptor_bridge_waters" not in receptor_variants:
            receptor_variants["receptor_bridge_waters"] = []
        receptor_variants["receptor_bridge_waters"].append(files_bridge)
        idx_alt += 1
    model = st[0]
    for chain in model:
        chain_id = chain.name
        if chain_id == "":
            chain_id = "?"
        chain_polymer: set[Tuple[str, int, str, str]] = set()
        for res in chain:
            rid = _res_id_tuple(chain.name, res)
            if rid in polymer:
                chain_polymer.add(rid)
        if len(chain_polymer) == 0:
            continue
        chain_variants[chain_id] = {}
        idx_alt2 = 0
        while idx_alt2 < len(altloc_choices):
            alt2 = altloc_choices[idx_alt2]
            alt_tag2 = alt_label_list[idx_alt2]
            ch_only = _extract_structure_variant(st, chain_polymer, altloc_keep=alt2, keep_hydrogens=keep_hydrogens)
            ch_plus_set: set[Tuple[str, int, str, str]] = set()
            for rid in chain_polymer:
                ch_plus_set.add(rid)
            for rid in cofactors:
                ch_plus_set.add(rid)
            for rid in ions_set:
                ch_plus_set.add(rid)
            ch_plus = _extract_structure_variant(st, ch_plus_set, altloc_keep=alt2, keep_hydrogens=keep_hydrogens)
            ch_bridge_set: set[Tuple[str, int, str, str]] = set()
            for rid in ch_plus_set:
                ch_bridge_set.add(rid)
            for rid in bridge_waters:
                ch_bridge_set.add(rid)
            ch_bridge = _extract_structure_variant(st, ch_bridge_set, altloc_keep=alt2, keep_hydrogens=keep_hydrogens)
            chain_dir = os.path.join(out_root, pdb_code, dataset, "chains")
            ensure_dir(chain_dir)
            base_name = pdb_code + "_" + dataset + "_chain_" + _sanitize_filename(chain_id) + "_" + alt_tag2
            f_only = StructureFiles(cif_path=os.path.join(chain_dir, base_name + "_only.cif"), pdb_path=os.path.join(chain_dir, base_name + "_only.pdb"))
            f_plus = StructureFiles(cif_path=os.path.join(chain_dir, base_name + "_plus_cofactors.cif"), pdb_path=os.path.join(chain_dir, base_name + "_plus_cofactors.pdb"))
            f_bridge = StructureFiles(cif_path=os.path.join(chain_dir, base_name + "_bridge_waters.cif"), pdb_path=os.path.join(chain_dir, base_name + "_bridge_waters.pdb"))
            _write_structure_files(ch_only, f_only.cif_path, f_only.pdb_path, warnings)
            _write_structure_files(ch_plus, f_plus.cif_path, f_plus.pdb_path, warnings)
            _write_structure_files(ch_bridge, f_bridge.cif_path, f_bridge.pdb_path, warnings)
            if "chain_only" not in chain_variants[chain_id]:
                chain_variants[chain_id]["chain_only"] = []
            chain_variants[chain_id]["chain_only"].append(f_only)
            if "chain_plus_cofactors" not in chain_variants[chain_id]:
                chain_variants[chain_id]["chain_plus_cofactors"] = []
            chain_variants[chain_id]["chain_plus_cofactors"].append(f_plus)
            if "chain_bridge_waters" not in chain_variants[chain_id]:
                chain_variants[chain_id]["chain_bridge_waters"] = []
            chain_variants[chain_id]["chain_bridge_waters"].append(f_bridge)
            idx_alt2 += 1
    lig_dir = os.path.join(out_root, pdb_code, dataset, "ligands")
    cof_dir = os.path.join(out_root, pdb_code, dataset, "cofactors")
    ion_dir = os.path.join(out_root, pdb_code, dataset, "ions")
    ensure_dir(lig_dir)
    ensure_dir(cof_dir)
    ensure_dir(ion_dir)
    ligand_rids_sorted = sorted(list(ligands), key=lambda x: (x[3], x[0], x[1], x[2]))
    cofactor_rids_sorted = sorted(list(cofactors), key=lambda x: (x[3], x[0], x[1], x[2]))
    ion_rids_sorted = sorted(list(ions_set), key=lambda x: (x[3], x[0], x[1], x[2]))
    for rid in ligand_rids_sorted:
        res_st = _extract_residue_as_structure(st, rid, altloc_keep=None, keep_hydrogens=keep_hydrogens)
        if res_st is None:
            warnings.add("Failed to extract ligand residue.", rid=str(rid))
            continue
        altloc_choices_res = _altloc_plan(res_st, altloc_mode, keep_hydrogens=keep_hydrogens, warnings=warnings)
        idx_a = 0
        while idx_a < len(altloc_choices_res):
            alt = altloc_choices_res[idx_a]
            alt_tag = "altNONE"
            if alt is not None:
                alt_tag = "alt" + str(alt)
            res_st2 = _extract_residue_as_structure(st, rid, altloc_keep=alt, keep_hydrogens=keep_hydrogens)
            if res_st2 is None:
                idx_a += 1
                continue
            tag = _residue_file_tag(rid) + "_" + alt_tag
            f = StructureFiles(cif_path=os.path.join(lig_dir, "ligand_" + tag + ".cif"), pdb_path=os.path.join(lig_dir, "ligand_" + tag + ".pdb"))
            _write_structure_files(res_st2, f.cif_path, f.pdb_path, warnings)
            lig_files.append(f)
            idx_a += 1
    for rid in cofactor_rids_sorted:
        res_st = _extract_residue_as_structure(st, rid, altloc_keep=None, keep_hydrogens=keep_hydrogens)
        if res_st is None:
            warnings.add("Failed to extract cofactor residue.", rid=str(rid))
            continue
        altloc_choices_res = _altloc_plan(res_st, altloc_mode, keep_hydrogens=keep_hydrogens, warnings=warnings)
        idx_a = 0
        while idx_a < len(altloc_choices_res):
            alt = altloc_choices_res[idx_a]
            alt_tag = "altNONE"
            if alt is not None:
                alt_tag = "alt" + str(alt)
            res_st2 = _extract_residue_as_structure(st, rid, altloc_keep=alt, keep_hydrogens=keep_hydrogens)
            if res_st2 is None:
                idx_a += 1
                continue
            tag = _residue_file_tag(rid) + "_" + alt_tag
            f = StructureFiles(cif_path=os.path.join(cof_dir, "cofactor_" + tag + ".cif"), pdb_path=os.path.join(cof_dir, "cofactor_" + tag + ".pdb"))
            _write_structure_files(res_st2, f.cif_path, f.pdb_path, warnings)
            cof_files.append(f)
            idx_a += 1
    for rid in ion_rids_sorted:
        res_st = _extract_residue_as_structure(st, rid, altloc_keep=None, keep_hydrogens=keep_hydrogens)
        if res_st is None:
            warnings.add("Failed to extract ion residue.", rid=str(rid))
            continue
        altloc_choices_res = _altloc_plan(res_st, altloc_mode, keep_hydrogens=keep_hydrogens, warnings=warnings)
        idx_a = 0
        while idx_a < len(altloc_choices_res):
            alt = altloc_choices_res[idx_a]
            alt_tag = "altNONE"
            if alt is not None:
                alt_tag = "alt" + str(alt)
            res_st2 = _extract_residue_as_structure(st, rid, altloc_keep=alt, keep_hydrogens=keep_hydrogens)
            if res_st2 is None:
                idx_a += 1
                continue
            tag = _residue_file_tag(rid) + "_" + alt_tag
            f = StructureFiles(cif_path=os.path.join(ion_dir, "ion_" + tag + ".cif"), pdb_path=os.path.join(ion_dir, "ion_" + tag + ".pdb"))
            _write_structure_files(res_st2, f.cif_path, f.pdb_path, warnings)
            ion_files.append(f)
            idx_a += 1
    meta_path = os.path.join(out_root, pdb_code, dataset, "index.json")
    idx = OutputIndex(receptor_variants=receptor_variants, chain_variants=chain_variants, ligands=lig_files, cofactors=cof_files, ions=ion_files, metadata_json=meta_path)
    return idx

def build_metadata(pdb_code: str, entry_meta: Optional[dict[str, Any]], au_files: Optional[OutputIndex], ba1_files: Optional[OutputIndex], out_root: str, settings: dict[str, Any], warnings: WarningCollector) -> dict[str, Any]:
    pid = normalize_pdb_id(pdb_code)
    title = ""
    deposit_date = ""
    release_date = ""
    revision_date = ""
    experimental_method = []
    resolution = None
    struct_keywords = ""
    polymer_entity_count = None
    nonpolymer_entity_count = None
    if entry_meta is not None and isinstance(entry_meta, dict):
        try:
            if "struct" in entry_meta and isinstance(entry_meta["struct"], dict):
                if "title" in entry_meta["struct"] and entry_meta["struct"]["title"] is not None:
                    title = str(entry_meta["struct"]["title"])
        except Exception:
            pass
        try:
            if "rcsb_accession_info" in entry_meta and isinstance(entry_meta["rcsb_accession_info"], dict):
                ai = entry_meta["rcsb_accession_info"]
                if "deposit_date" in ai and ai["deposit_date"] is not None:
                    deposit_date = str(ai["deposit_date"])
                if "initial_release_date" in ai and ai["initial_release_date"] is not None:
                    release_date = str(ai["initial_release_date"])
                if "revision_date" in ai and ai["revision_date"] is not None:
                    revision_date = str(ai["revision_date"])
        except Exception:
            pass
        try:
            if "rcsb_entry_info" in entry_meta and isinstance(entry_meta["rcsb_entry_info"], dict):
                ei = entry_meta["rcsb_entry_info"]
                if "experimental_method" in ei and ei["experimental_method"] is not None:
                    experimental_method = []
                    if isinstance(ei["experimental_method"], list):
                        for m in ei["experimental_method"]:
                            experimental_method.append(str(m))
                    else:
                        experimental_method.append(str(ei["experimental_method"]))
                if "resolution_combined" in ei and ei["resolution_combined"] is not None:
                    if isinstance(ei["resolution_combined"], list):
                        if len(ei["resolution_combined"]) > 0:
                            try:
                                resolution = float(ei["resolution_combined"][0])
                            except Exception:
                                resolution = None
                    else:
                        try:
                            resolution = float(ei["resolution_combined"])
                        except Exception:
                            resolution = None
                if "polymer_entity_count" in ei and ei["polymer_entity_count"] is not None:
                    try:
                        polymer_entity_count = int(ei["polymer_entity_count"])
                    except Exception:
                        polymer_entity_count = None
                if "nonpolymer_entity_count" in ei and ei["nonpolymer_entity_count"] is not None:
                    try:
                        nonpolymer_entity_count = int(ei["nonpolymer_entity_count"])
                    except Exception:
                        nonpolymer_entity_count = None
        except Exception:
            pass
        try:
            if "struct_keywords" in entry_meta and isinstance(entry_meta["struct_keywords"], dict):
                sk = entry_meta["struct_keywords"]
                if "pdbx_keywords" in sk and sk["pdbx_keywords"] is not None:
                    struct_keywords = str(sk["pdbx_keywords"])
        except Exception:
            pass
    files: dict[str, Any] = {"asymmetric_unit": {}, "biological_assembly_1": {}}
    def collect_index(idx: OutputIndex) -> dict[str, Any]:
        out: dict[str, Any] = {}
        out["receptor_variants"] = {}
        for k, lst in idx.receptor_variants.items():
            out["receptor_variants"][k] = []
            for f in lst:
                out["receptor_variants"][k].append({"cif": safe_relpath(f.cif_path, out_root), "pdb": safe_relpath(f.pdb_path, out_root), "sha256_cif": sha256_file(f.cif_path) if os.path.exists(f.cif_path) else "", "sha256_pdb": sha256_file(f.pdb_path) if os.path.exists(f.pdb_path) else ""})
        out["chains"] = {}
        for chain_id, variants in idx.chain_variants.items():
            out["chains"][chain_id] = {}
            for k2, lst2 in variants.items():
                out["chains"][chain_id][k2] = []
                for f2 in lst2:
                    out["chains"][chain_id][k2].append({"cif": safe_relpath(f2.cif_path, out_root), "pdb": safe_relpath(f2.pdb_path, out_root), "sha256_cif": sha256_file(f2.cif_path) if os.path.exists(f2.cif_path) else "", "sha256_pdb": sha256_file(f2.pdb_path) if os.path.exists(f2.pdb_path) else ""})
        out["ligands"] = []
        for f3 in idx.ligands:
            out["ligands"].append({"cif": safe_relpath(f3.cif_path, out_root), "pdb": safe_relpath(f3.pdb_path, out_root), "sha256_cif": sha256_file(f3.cif_path) if os.path.exists(f3.cif_path) else "", "sha256_pdb": sha256_file(f3.pdb_path) if os.path.exists(f3.pdb_path) else ""})
        out["cofactors"] = []
        for f4 in idx.cofactors:
            out["cofactors"].append({"cif": safe_relpath(f4.cif_path, out_root), "pdb": safe_relpath(f4.pdb_path, out_root), "sha256_cif": sha256_file(f4.cif_path) if os.path.exists(f4.cif_path) else "", "sha256_pdb": sha256_file(f4.pdb_path) if os.path.exists(f4.pdb_path) else ""})
        out["ions"] = []
        for f5 in idx.ions:
            out["ions"].append({"cif": safe_relpath(f5.cif_path, out_root), "pdb": safe_relpath(f5.pdb_path, out_root), "sha256_cif": sha256_file(f5.cif_path) if os.path.exists(f5.cif_path) else "", "sha256_pdb": sha256_file(f5.pdb_path) if os.path.exists(f5.pdb_path) else ""})
        return out
    if au_files is not None:
        files["asymmetric_unit"] = collect_index(au_files)
    if ba1_files is not None:
        files["biological_assembly_1"] = collect_index(ba1_files)
    meta: dict[str, Any] = {"schema_version": "pdbstex_metadata_v1", "Name": title, "PDB_code": pid, "deposit_date": deposit_date, "release_date": release_date, "revision_date": revision_date, "experimental_method": experimental_method, "resolution_angstrom": resolution, "keywords": struct_keywords, "polymer_entity_count": polymer_entity_count, "nonpolymer_entity_count": nonpolymer_entity_count, "settings": settings, "files": files, "warnings": warnings.items()}
    return meta

#!/usr/bin/env python3
"""Chemistry-aware conformer deduplication with RDKit graph matching."""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdMolAlign


def read_sdf(paths: list[Path]) -> list[tuple[str, int, Chem.Mol]]:
    mols: list[tuple[str, int, Chem.Mol]] = []
    for path in paths:
        supplier = Chem.SDMolSupplier(str(path), removeHs=False, sanitize=True)
        for idx, mol in enumerate(supplier):
            if mol is None:
                raise ValueError(f"RDKit could not parse molecule {idx} in {path}")
            mols.append((str(path), idx, mol))
    return mols


def full_match_to_template(template: Chem.Mol, mol: Chem.Mol) -> tuple[int, ...]:
    if template.GetNumAtoms() != mol.GetNumAtoms():
        raise ValueError("atom count differs from first molecule")
    match = mol.GetSubstructMatch(template)
    if len(match) != template.GetNumAtoms():
        raise ValueError("molecule graph does not fully match the first molecule")
    return match


def renumber_to_template(template: Chem.Mol, mol: Chem.Mol) -> Chem.Mol:
    match = full_match_to_template(template, mol)
    return Chem.RenumberAtoms(mol, list(match))


def mol_name(mol: Chem.Mol, fallback: str) -> str:
    return mol.GetProp("_Name") if mol.HasProp("_Name") else fallback


def main() -> int:
    parser = argparse.ArgumentParser(description="Deduplicate SDF conformers using RDKit graph-aware RMSD.")
    parser.add_argument("sdf", nargs="+", type=Path, help="Input SDF file(s).")
    parser.add_argument("--tolerance", type=float, default=1e-3, help="RMSD duplicate tolerance in Angstrom.")
    parser.add_argument("--write-unique", type=Path, help="Optional output SDF for unique representatives.")
    args = parser.parse_args()

    if args.tolerance < 0:
        raise ValueError("--tolerance must be non-negative")

    records = read_sdf(args.sdf)
    if not records:
        raise ValueError("no molecules were read")

    template = records[0][2]
    unique: list[tuple[int, str, int, Chem.Mol]] = []
    duplicates: list[tuple[int, str, int, int, float]] = []

    start = time.perf_counter()
    for input_index, (path, mol_index, mol) in enumerate(records):
        try:
            aligned_mol = renumber_to_template(template, mol)
        except ValueError as exc:
            raise ValueError(f"{path} molecule {mol_index}: {exc}") from exc

        best_unique_idx = -1
        best_rmsd = float("inf")
        for unique_idx, (_, _, _, representative) in enumerate(unique):
            rmsd = rdMolAlign.GetBestRMS(aligned_mol, representative)
            if rmsd < best_rmsd:
                best_rmsd = rmsd
                best_unique_idx = unique_idx

        if best_unique_idx >= 0 and best_rmsd <= args.tolerance:
            representative_input_index = unique[best_unique_idx][0]
            duplicates.append((input_index, path, mol_index, representative_input_index, best_rmsd))
        else:
            unique.append((input_index, path, mol_index, aligned_mol))

    elapsed = time.perf_counter() - start

    if args.write_unique:
        writer = Chem.SDWriter(str(args.write_unique))
        for unique_index, (input_index, path, mol_index, mol) in enumerate(unique):
            out = Chem.Mol(mol)
            out.SetProp("_Name", mol_name(out, f"unique_{unique_index + 1:04d}"))
            out.SetIntProp("source_input_index", input_index)
            out.SetProp("source_path", path)
            out.SetIntProp("source_molecule_index", mol_index)
            writer.write(out)
        writer.close()

    print(f"input_count {len(records)}")
    print(f"unique_count {len(unique)}")
    print(f"duplicate_count {len(duplicates)}")
    print(f"tolerance_angstrom {args.tolerance:.10f}")
    print(f"elapsed_seconds {elapsed:.10f}")
    print(f"rdkit_version {Chem.rdBase.rdkitVersion}")
    for unique_index, (input_index, path, mol_index, mol) in enumerate(unique):
        print(
            "unique"
            f" {unique_index} input_index {input_index} path {path} molecule_index {mol_index}"
            f" name {mol_name(mol, '')}"
        )
    for input_index, path, mol_index, representative_input_index, rmsd in duplicates:
        print(
            "duplicate"
            f" input_index {input_index} path {path} molecule_index {mol_index}"
            f" representative_input_index {representative_input_index} rmsd_angstrom {rmsd:.10f}"
        )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"error: {exc}", file=sys.stderr)
        raise SystemExit(2)

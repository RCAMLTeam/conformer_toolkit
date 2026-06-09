#!/usr/bin/env python3
"""Benchmark conformer_identical against RDKit rdMolAlign.GetBestRMS."""

from __future__ import annotations

import argparse
import math
import re
import subprocess
import time
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign


def write_xyz(mol: Chem.Mol, path: Path) -> None:
    conf = mol.GetConformer()
    lines = [str(mol.GetNumAtoms()), path.stem]
    for atom in mol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        lines.append(f"{atom.GetSymbol()} {pos.x:.10f} {pos.y:.10f} {pos.z:.10f}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def transformed_copy(mol: Chem.Mol) -> Chem.Mol:
    # Apply only a rigid transform. If this benchmark reports nonzero RMSD, the
    # transform or atom ordering changed internal geometry by mistake.
    out = Chem.Mol(mol)
    conf = out.GetConformer()
    angle = math.radians(37.0)
    cos_a = math.cos(angle)
    sin_a = math.sin(angle)
    for i in range(out.GetNumAtoms()):
        pos = conf.GetAtomPosition(i)
        x = cos_a * pos.x - sin_a * pos.y + 1.25
        y = sin_a * pos.x + cos_a * pos.y - 0.50
        z = pos.z + 0.75
        conf.SetAtomPosition(i, (x, y, z))
    return out


def parse_avg_microseconds(output: str) -> float:
    match = re.search(r"^avg_microseconds ([0-9.eE+-]+)$", output, re.MULTILINE)
    if not match:
        raise RuntimeError(f"C++ benchmark output did not contain avg_microseconds:\n{output}")
    return float(match.group(1))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repeat", type=int, default=100_000)
    parser.add_argument("--binary", type=Path, default=Path("src/conformer_identical"))
    parser.add_argument("--outdir", type=Path, default=Path("src/benchmark_data"))
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    mol1 = Chem.AddHs(Chem.MolFromSmiles("CCO"))
    status = AllChem.EmbedMolecule(mol1, randomSeed=7)
    if status != 0:
        raise RuntimeError("RDKit failed to embed ethanol")
    AllChem.UFFOptimizeMolecule(mol1)
    mol2 = transformed_copy(mol1)

    xyz1 = args.outdir / "ethanol_a.xyz"
    xyz2 = args.outdir / "ethanol_rot_trans.xyz"
    write_xyz(mol1, xyz1)
    write_xyz(mol2, xyz2)

    cpp = subprocess.run(
        [str(args.binary), str(xyz1), str(xyz2), "--repeat", str(args.repeat)],
        check=True,
        text=True,
        capture_output=True,
    )
    cpp_avg_us = parse_avg_microseconds(cpp.stdout)

    # RDKit does graph/symmetry work in GetBestRMS, so this benchmark measures
    # fixed-order RMSD speed versus a more general chemistry-aware operation.
    start = time.perf_counter()
    rmsd = 0.0
    for _ in range(args.repeat):
        rmsd = rdMolAlign.GetBestRMS(mol1, mol2)
    elapsed = time.perf_counter() - start
    rdkit_avg_us = elapsed * 1e6 / args.repeat

    print(f"repeat {args.repeat}")
    print(f"rdkit_version {Chem.rdBase.rdkitVersion}")
    print(f"molecule ethanol_explicit_h atoms {mol1.GetNumAtoms()}")
    print(f"cpp_avg_microseconds {cpp_avg_us:.6f}")
    print(f"rdkit_getbestrms_avg_microseconds {rdkit_avg_us:.6f}")
    print(f"ratio_rdkit_over_cpp {rdkit_avg_us / cpp_avg_us:.3f}")
    print(f"rdkit_rmsd_angstrom {rmsd:.10f}")
    print("cpp_output_begin")
    print(cpp.stdout.strip())
    print("cpp_output_end")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

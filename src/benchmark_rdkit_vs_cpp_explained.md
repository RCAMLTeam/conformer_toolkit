# `benchmark_rdkit_vs_cpp.py` Debug Notes

This script compares the fixed-order C++ RMSD loop against RDKit `rdMolAlign.GetBestRMS`.

## What It Benchmarks

1. Build an explicit-H ethanol molecule with RDKit.
2. Generate a rotated/translated copy.
3. Write both structures to XYZ for the C++ tool.
4. Run:

```bash
src/conformer_identical ethanol_a.xyz ethanol_rot_trans.xyz --repeat N
```

5. Run `rdMolAlign.GetBestRMS(mol1, mol2)` `N` times in Python.
6. Print average microseconds per comparison.

## Interpretation

The C++ path is expected to be faster because it assumes fixed atom order and only solves rigid alignment RMSD.

RDKit is expected to be slower because `GetBestRMS` handles graph/symmetry-aware matching.

The timings are therefore not an apples-to-apples chemistry benchmark. They answer:

```text
How much faster is the fixed-order C++ RMSD loop when its assumptions are valid?
```

## Common Debug Points

- If the C++ binary is missing, rebuild it from the package root.
- If RDKit import fails, use the local environment at `./conformer_toolkit/bin/python`.
- If RMSD is not near zero, check whether the test molecule transform changed internal geometry rather than only rotating/translating.

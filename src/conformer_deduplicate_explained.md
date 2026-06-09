# `conformer_deduplicate.cpp` Debug Notes

This file implements batch deduplication for XYZ conformers.

## Main Flow

1. Parse CLI options in `main`.
2. Read the first XYZ as the reference format.
3. For each input conformer:
   - strict mode: require same atom symbols in the same order
   - `--allow-reorder`: require same formula, then search graph-preserving atom mappings
4. Compare the conformer against each unique representative.
5. If best RMSD is within tolerance, record it as duplicate.
6. Otherwise keep it as a new unique representative.
7. Optionally write unique XYZ files with `--write-unique`.

## Strict Mode

Strict mode uses:

```cpp
assert_same_atom_indexing(reference, mol, path)
aligned_rmsd(candidate.mol, mol)
```

This is the fastest path. It assumes atom index `i` in every XYZ refers to the same chemical atom.

## `--allow-reorder` Mode

This mode is for XYZ files where atom order differs. XYZ has no bond graph, so the code infers one:

```cpp
infer_graph(mol, bond_scale)
```

The default `bond_scale` is `1.1`. Two atoms are bonded when:

```text
distance <= bond_scale * (radius_a + radius_b)
```

Then each atom gets a local signature:

```text
element | degree | sorted neighbor elements
```

Example:

```text
C|degree=4|neighbors=C,H,H,H,
C|degree=2|neighbors=C,C,
O|degree=1|neighbors=C,
```

Candidate atom mappings are only allowed when signatures match. During recursive search, the code also checks that assigned atoms preserve every bond/non-bond relationship.

The key function is:

```cpp
best_graph_preserving_reordered_rmsd(...)
```

## RMSD Algorithm

`aligned_rmsd`:

1. Centers both coordinate sets.
2. Builds the covariance terms.
3. Uses Horn's quaternion formulation.
4. Estimates the largest eigenvalue with power iteration.
5. Converts that into optimal rigid-alignment RMSD.

If RMSD results look suspicious, test `aligned_rmsd` first with translated/rotated copies before debugging mapping logic.

## Common Debug Points

- If `--allow-reorder` says no valid mapping exists, inspect inferred bonds. The default `--bond-scale 1.1` may be too strict for stretched bonds.
- If too many mappings are checked, the molecule is probably symmetric or has many repeated local signatures. Use RDKit SDF mode when possible.
- If an expected duplicate is missed, compare:
  - formula match
  - inferred graph signatures
  - `--bond-scale`
  - final RMSD tolerance
- If a chemically invalid mapping is accepted, the local graph signatures may be too coarse. The next improvement would be deeper graph invariants or RDKit-generated mappings.

## Important Limitations

This C++ tool infers connectivity from XYZ distances. It does not know bond order, formal charge, aromaticity, stereochemistry, or resonance. For chemically robust matching, prefer `deduplicate_rdkit.py` with SDF/MOL input.

# `deduplicate_rdkit.py` Debug Notes

This file implements chemistry-aware deduplication for SDF inputs using RDKit.

## Main Flow

1. Read one or more SDF files with:

```python
Chem.SDMolSupplier(..., removeHs=False, sanitize=True)
```

2. Use the first molecule as the reference graph.
3. For each molecule:
   - find a full graph match to the reference
   - renumber atoms to reference order
   - compare against existing unique representatives with `rdMolAlign.GetBestRMS`
4. Keep the first representative of each unique conformer group.
5. Optionally write unique representatives to SDF.

## Atom Mapping

The mapping function is:

```python
full_match_to_template(template, mol)
```

It calls:

```python
mol.GetSubstructMatch(template)
```

The returned tuple maps:

```text
template atom index -> mol atom index
```

Then:

```python
Chem.RenumberAtoms(mol, list(match))
```

creates a molecule whose atom order follows the template.

## RMSD

The duplicate test uses:

```python
rdMolAlign.GetBestRMS(aligned_mol, representative)
```

RDKit handles graph/symmetry matching internally for RMSD, so this path is more chemically robust than XYZ-only matching.

## Common Debug Points

- If parsing fails, inspect the SDF for invalid valence, missing coordinate blocks, or unsupported records.
- If graph matching fails, the molecule may not be the same chemical graph as the first SDF molecule.
- If hydrogens differ, matching can fail because `removeHs=False` preserves explicit hydrogens. Make sure all conformers use the same hydrogen representation.
- If deduplication is slower than C++, that is expected. RDKit is doing graph-aware work.

## Known Limitation

`GetSubstructMatch` returns one valid match, not necessarily the lowest-RMSD mapping among all symmetry-allowed mappings. `GetBestRMS` still performs symmetry-aware RMSD for duplicate classification, but the SDF atom renumbering step may not always choose the same mapping that minimizes RMSD. If exact reference-order output matters for highly symmetric molecules, extend this script to enumerate `GetSubstructMatches(..., uniquify=False)`, score each mapping with RMSD, and renumber using the best mapping.

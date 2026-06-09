# Conformer Toolkit

Tools for checking whether conformers are duplicates and reducing a conformer set to unique representatives.

The toolkit has two deduplication paths:

- Fast C++ XYZ path for ordered conformers, with an optional geometry-only atom reordering mode.
- RDKit SDF path for chemistry-aware atom matching when molecular connectivity is available.

## Algorithms

### Ordered XYZ RMSD

`conformer_identical` and the default `conformer_deduplicate` mode assume every XYZ file uses the same atom index order.

For each comparison, the code:

1. Reads atom symbols and Cartesian coordinates from XYZ.
2. Verifies atom count and atom symbols match in the same order.
3. Removes translation by centering both conformers.
4. Computes the optimal rigid rotation using Horn's quaternion RMSD formulation.
5. Classifies conformers as duplicates when aligned RMSD is less than or equal to the tolerance.

This is the fastest and safest mode when conformers came from the same molecule file or workflow and atom ordering is preserved.

### XYZ Geometry-Only Reordering

`conformer_deduplicate --allow-reorder` supports XYZ files whose atoms are not in the same order.

For each comparison, the code:

1. Verifies the two conformers have the same element counts.
2. Searches element-preserving atom mappings.
3. Computes aligned RMSD for each mapping.
4. Uses the lowest RMSD mapping for duplicate detection.

This mode is useful for XYZ-only data, but it has no bond graph. Symmetric molecules or many repeated atoms can produce ambiguous mappings. If mapping search becomes too large, increase `--max-mappings` or use the RDKit path with connectivity.

### RDKit Graph-Aware Deduplication

`deduplicate_rdkit.py` is the preferred path when you have SDF/MOL-style connectivity.

For each molecule, the script:

1. Reads molecules with explicit hydrogens preserved.
2. Uses RDKit substructure graph matching against the first molecule.
3. Renumbers atoms to the first molecule's atom order.
4. Compares conformers with `rdMolAlign.GetBestRMS`.
5. Keeps the first representative of each unique conformer group.

This is slower than the fixed-order C++ path, but it is chemically safer when atom order is inconsistent.

## Build

From the package root:

```bash
g++ -O3 -std=c++17 -Wall -Wextra -pedantic src/conformer_identical.cpp -o src/conformer_identical
g++ -O3 -std=c++17 -Wall -Wextra -pedantic src/conformer_deduplicate.cpp -o src/conformer_deduplicate
```

The RDKit script uses the existing local Python environment:

```bash
./conformer_toolkit/bin/python -m py_compile src/deduplicate_rdkit.py
```

## Usage

### Compare Two Ordered XYZ Conformers

```bash
src/conformer_identical conformer_a.xyz conformer_b.xyz
src/conformer_identical conformer_a.xyz conformer_b.xyz 0.0005
```

Default tolerance is `1e-3` Angstrom.

Exit codes:

- `0`: identical
- `1`: different
- `2`: usage or input error

### Deduplicate Ordered XYZ Conformers

```bash
src/conformer_deduplicate --tolerance 0.001 conformer_*.xyz
```

Write unique representatives:

```bash
src/conformer_deduplicate --tolerance 0.001 --write-unique unique_xyz conformer_*.xyz
```

Default behavior requires the same atom symbols in the same order for every input. If atom order differs, the program exits with an error.

### Deduplicate XYZ Conformers With Reordered Atoms

```bash
src/conformer_deduplicate --allow-reorder --tolerance 0.001 conformer_*.xyz
```

Limit or expand the mapping search:

```bash
src/conformer_deduplicate --allow-reorder --max-mappings 10000000 conformer_*.xyz
```

Use this only when you have XYZ files and no connectivity. For symmetric or highly repetitive structures, prefer RDKit with SDF input.

### Deduplicate SDF With RDKit

```bash
./conformer_toolkit/bin/python src/deduplicate_rdkit.py --tolerance 0.001 input.sdf
```

Write unique representatives to an SDF:

```bash
./conformer_toolkit/bin/python src/deduplicate_rdkit.py --tolerance 0.001 --write-unique unique.sdf input.sdf
```

Multiple SDF files can be passed:

```bash
./conformer_toolkit/bin/python src/deduplicate_rdkit.py --tolerance 0.001 batch_1.sdf batch_2.sdf
```

## Output

`conformer_deduplicate` prints a line-based summary:

```text
input_count 3
unique_count 2
duplicate_count 1
tolerance_angstrom 0.0010000000
unique 0 input_index 0 path conformer_001.xyz
duplicate input_index 1 path conformer_002.xyz representative_input_index 0 representative_path conformer_001.xyz rmsd_angstrom 0.0000000000
```

`deduplicate_rdkit.py` prints similar counts and includes RDKit version plus SDF molecule indices.

## Benchmark

Benchmark the C++ fixed-order RMSD loop against RDKit `GetBestRMS`:

```bash
./conformer_toolkit/bin/python src/benchmark_rdkit_vs_cpp.py --repeat 100000
```

In the local ethanol test, the C++ fixed-order loop was about 9x faster than `rdMolAlign.GetBestRMS`. RDKit does more chemistry-aware work, so use the faster C++ path only when its atom-order assumptions are valid.

## Limitations

- XYZ files do not contain reliable connectivity. Geometry-only reordering can be ambiguous.
- The C++ tools compare conformers by RMSD only; they do not check molecular charge, bond order, stereochemistry, or energy.
- `--allow-reorder` performs exhaustive element-wise mapping search and can become expensive for many identical atoms.
- RDKit graph-aware mode requires inputs that RDKit can parse with connectivity, such as SDF/MOL.

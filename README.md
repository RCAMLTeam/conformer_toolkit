# Conformer Toolkit

Tools for checking whether conformers are duplicates and reducing a conformer set to unique representatives.

The toolkit has two deduplication paths:

- Fast C++ XYZ path for ordered conformers, with an optional graph-aware atom reordering mode from inferred bonds.
- RDKit SDF path for chemistry-aware atom matching when molecular connectivity is available.

## Algorithms

### Ordered XYZ RMSD

`conformer_identical` and the default `conformer_deduplicate` mode assume every XYZ file uses the same atom index order. RDKit C++ is used by default to infer and verify molecular chemistry before RMSD comparison.

For each comparison, the code:

1. Reads atom symbols and Cartesian coordinates from XYZ.
2. Uses RDKit `DetermineBonds`/xyz2mol to infer connectivity, bond orders, formal charges, and 3D stereochemistry.
3. Verifies atom count and atom symbols match in the same order.
4. Requires a full RDKit graph match with chirality enabled.
5. Removes translation and computes the optimal proper rotation using Horn's quaternion RMSD formulation.
6. Classifies conformers as duplicates only when chemistry and stereochemistry match and RMSD is less than or equal to the tolerance.

This is the fastest and safest mode when conformers came from the same molecule file or workflow and atom ordering is preserved.

### RDKit XYZ Graph-Aware Reordering

`conformer_deduplicate --allow-reorder` supports XYZ files whose atoms are not in the same order.

For each comparison, the code:

1. Verifies the two conformers have the same element counts.
2. Uses RDKit `DetermineBonds` to infer full molecular graphs from XYZ coordinates.
3. Uses RDKit substructure matching to enumerate topology-preserving atom mappings.
4. Chooses the lowest-RMSD topology mapping to normalize atom order.
5. Performs a second RDKit match with chirality enabled during duplicate classification.
6. Keeps opposite tetrahedral or double-bond stereoisomers as distinct records.

Topology normalization intentionally ignores stereochemistry so stereoisomers can share a consistent atom order. Duplicate classification does not ignore stereochemistry.

### RDKit Graph-Aware Deduplication

`deduplicate_rdkit.py` remains available for SDF/MOL inputs that already contain explicit connectivity.

For each molecule, the script:

1. Reads molecules with explicit hydrogens preserved.
2. Uses RDKit substructure graph matching against the first molecule.
3. Renumbers atoms to the first molecule's atom order.
4. Compares conformers with `rdMolAlign.GetBestRMS`.
5. Keeps the first representative of each unique conformer group.

This is slower than the fixed-order C++ path, but it is chemically safer when atom order is inconsistent.

## Build

The C++ tools require RDKit C++ development headers/libraries, pybind11 headers for the Python extension, and the RDKit `DetermineBonds` module.

```bash
sudo apt-get install cmake g++ librdkit-dev pybind11-dev python3-dev
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

Some packaged RDKit builds, including Debian's, omit the `DetermineBonds` library and header. Use a matching RDKit source checkout in that case:

```bash
git clone --depth 1 --branch Release_2025_03_1 https://github.com/rdkit/rdkit.git /path/to/rdkit-source
cmake -S . -B build \
  -DRDKIT_ROOT=/usr \
  -DRDKIT_SOURCE_ROOT=/path/to/rdkit-source \
  -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
```

For a non-system RDKit installation, set `RDKIT_ROOT` to its prefix. The helper builds the Python extension through the same CMake configuration:

```bash
python3 src/build_pybind.py \
  --rdkit-root /path/to/rdkit-prefix \
  --rdkit-source-root /path/to/matching-rdkit-source
```

## Usage

The C++ batch tool is built around a `Conformer_Batch` class. It can load conformers from explicit XYZ file paths, all `.xyz` files under a directory, or a concatenated multi-XYZ file. The class exposes `index_cleanup` for atom-order normalization and `remove_duplicates`, which can optionally run `index_cleanup` before RMSD deduplication.

### Python Interface

After building the pybind11 extension, import it with `src` on `PYTHONPATH`:

```bash
PYTHONPATH=src ./conformer_toolkit/bin/python
```

```python
from conformer_toolkit_cpp import Conformer_Batch

batch = Conformer_Batch.from_xyz_files([
    "conf_001.xyz",
    "conf_002.xyz",
])

result = batch.remove_duplicates(
    tolerance=0.001,
    run_index_cleanup=True,
    charge=0,
)

print(len(result.unique), len(result.duplicates))
for duplicate in result.duplicates:
    print(duplicate.path, duplicate.representative_path, duplicate.rmsd)
```

Other loaders:

```python
batch = Conformer_Batch.from_xyz_directory("conformer_dir")
batch = Conformer_Batch.from_multi_xyz("conformers.xyz")
```

Run index cleanup on its own:

```python
batch.index_cleanup(max_mappings=1_000_000, bond_scale=1.3, charge=0)
batch.write_records("cleaned_xyz", "cleaned")
```

### Compare Two Ordered XYZ Conformers

```bash
src/conformer_identical conformer_a.xyz conformer_b.xyz
src/conformer_identical conformer_a.xyz conformer_b.xyz 0.0005
src/conformer_identical conformer_a.xyz conformer_b.xyz --charge -1
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

Load all `.xyz` files under a directory:

```bash
src/conformer_deduplicate --xyz-dir conformer_dir --tolerance 0.001
```

Load a concatenated multi-XYZ file:

```bash
src/conformer_deduplicate --multi-xyz conformers.xyz --tolerance 0.001
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

`--allow-reorder` calls `Conformer_Batch::remove_duplicates` with its cleanup option enabled, so atom indexing is normalized before duplicate removal.

Run index cleanup on its own and write reordered XYZ files:

```bash
src/conformer_deduplicate --multi-xyz conformers.xyz --index-cleanup-only --write-cleaned cleaned_xyz
```

Write cleaned conformers while also removing duplicates:

```bash
src/conformer_deduplicate --allow-reorder --write-cleaned cleaned_xyz --write-unique unique_xyz conformer_*.xyz
```

Control inferred bonding and mapping search:

```bash
src/conformer_deduplicate --allow-reorder --bond-scale 1.3 conformer_*.xyz
src/conformer_deduplicate --allow-reorder --max-mappings 10000000 conformer_*.xyz
src/conformer_deduplicate --allow-reorder --charge 1 conformer_*.xyz
```

`--bond-scale` is passed to RDKit `DetermineBonds` as the covalent-radius multiplier. The default is `1.3`.

XYZ does not encode total molecular charge. `--charge` therefore defaults to `0` and must be provided for ions. One charge value applies to the complete input batch.

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

- XYZ files do not encode bonds, bond orders, formal charges, isotopes, or stereo annotations. RDKit reconstructs these from coordinates and the supplied total charge, so unusual, organometallic, highly distorted, or ambiguous structures can still be misassigned.
- Total charge cannot be inferred reliably from XYZ and must be passed with `--charge`.
- Stereo comparison is only as reliable as RDKit's inferred graph and 3D stereo assignment. Nearly planar or otherwise geometrically ambiguous centers may remain unspecified.
- `--allow-reorder` can become expensive for highly symmetric molecules. `--max-mappings` bounds RDKit's match enumeration.
- When multiple topology-compatible mappings exist, index cleanup chooses the lowest-RMSD mapping; duplicate classification subsequently requires stereochemistry to match.
- The C++ tools do not compare conformer energies.

## Source Debug Notes

Debugging notes are embedded directly in the source files as `//` comments for C++ and `#` comments for Python.

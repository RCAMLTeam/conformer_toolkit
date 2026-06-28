# Conformer Toolkit

Tools for checking whether conformers are duplicates and reducing a conformer set to unique representatives.

The toolkit has two deduplication paths:

- Fast C++ XYZ path for ordered conformers, with an optional graph-aware atom reordering mode from inferred bonds.
- RDKit SDF path for chemistry-aware atom matching when molecular connectivity is available.

All commands below are run from the repository root. CMake writes the two
executables and the optional Python extension to `src/`.

## Repository Layout

```text
.
├── CMakeLists.txt
├── README.md
└── src/
    ├── conformer_deduplicate.cpp
    ├── conformer_identical.cpp
    ├── conformer_toolkit_bindings.cpp
    ├── conformer_toolkit/
    │   └── deduplicate_rdkit.py
    └── testdata/
```

## Algorithms

### Ordered XYZ RMSD

`conformer_identical` and the default `conformer_deduplicate` mode assume every XYZ file uses the same atom index order. RDKit C++ is used by default to infer and verify molecular chemistry before RMSD comparison.

For each comparison, the code:

1. Reads atom symbols and Cartesian coordinates from XYZ.
2. Uses RDKit `DetermineBonds`/xyz2mol to infer connectivity, bond orders, formal charges, and 3D stereochemistry.
3. Verifies atom count and atom symbols match in the same order.
4. Requires a full RDKit graph match with chirality enabled.
5. Uses RDKit `MolAlign::alignMol` with an explicit same-index atom map to compute the optimal aligned RMSD.
6. Classifies conformers as duplicates only when chemistry and stereochemistry match and RMSD is less than or equal to the tolerance.

This is the simplest mode when conformers came from the same molecule file or workflow and atom ordering is preserved.

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

This remains available when SDF/MOL inputs already carry explicit connectivity.

## Build

Install the Python package from the repository checkout with:

```bash
python -m pip install .
```

This installs the RDKit SDF utilities and the
`conformer-deduplicate-sdf` command. The C++ binaries and optional pybind11
extension continue to use the CMake build described below.

The recommended build path uses Pixi, which installs Python, NumPy, RDKit, RDKit C++ headers/libraries, pybind11, CMake, and the compiler toolchain from conda-forge:

```bash
pixi install
pixi run smoke
```

Useful Pixi tasks:

```bash
pixi run configure
pixi run build
pixi run test
pixi run python-smoke
```

Pixi builds into `build-pixi/` and writes the command-line tools plus Python extension to `src/`.

The manual system build path requires a C++20-capable compiler, RDKit C++ development headers/libraries, pybind11 headers for the Python extension, and the RDKit `DetermineBonds` module.

```bash
sudo apt-get install cmake g++ librdkit-dev pybind11-dev python3-dev
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

After a successful build, the command-line tools are available as
`src/conformer_identical` and `src/conformer_deduplicate`. If Python bindings
are enabled, the `conformer_toolkit_cpp` extension is also written to `src/`.

Some packaged RDKit builds, including Debian's, omit the `DetermineBonds` library and header. Use a matching RDKit source checkout in that case:

```bash
git clone --depth 1 --branch Release_2025_03_1 https://github.com/rdkit/rdkit.git /path/to/rdkit-source
cmake -S . -B build \
  -DRDKIT_ROOT=/usr \
  -DRDKIT_SOURCE_ROOT=/path/to/rdkit-source \
  -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
```

On Debian 13 arm64, the RDKit and pybind11 packages can also be kept entirely under this repository:

```bash
mkdir -p third_party/debs third_party/deb-root
cd third_party/debs
apt-get download \
  librdkit-dev librdkit1t64 pybind11-dev \
  libboost-dev libboost1.83-dev \
  libboost-python1.83.0 libboost-serialization1.83.0 \
  libcoordgen3 libinchi1.07 libmaeparser1
for pkg in *.deb; do dpkg-deb -x "$pkg" ../deb-root; done
cd ../..
git clone --depth 1 --branch Release_2025_03_1 \
  https://github.com/rdkit/rdkit.git third_party/rdkit-source
cmake -S . -B build \
  -DRDKIT_ROOT="$PWD/third_party/deb-root/usr" \
  -DRDKIT_SOURCE_ROOT="$PWD/third_party/rdkit-source" \
  -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

For a non-system RDKit installation, set `RDKIT_ROOT` to its prefix. The helper builds the Python extension through the same CMake configuration:

```bash
python3 src/build_pybind.py \
  --rdkit-root /path/to/rdkit-prefix \
  --rdkit-source-root /path/to/matching-rdkit-source
```

## Usage

The Python extension exposes one C++ container:

- `Conformer_Group` stores all conformers for one molecule. It loads conformers from explicit XYZ file paths, all `.xyz` files under a directory, or a concatenated multi-XYZ file. It exposes `index_cleanup` for atom-order normalization, `remove_duplicates` for RMSD deduplication, and `detect_rings` for RDKit-derived ring information.

`Conformer_Batch` remains available as a Python compatibility alias for `Conformer_Group`.

### Python Interface

After building the pybind11 extension, import it with `src` on `PYTHONPATH`:

```bash
PYTHONPATH=src python3
```

```python
from conformer_toolkit_cpp import Conformer_Group

batch = Conformer_Group.from_xyz_files([
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
batch = Conformer_Group.from_xyz_directory("conformer_dir")
batch = Conformer_Group.from_multi_xyz("conformers.xyz")
```

XYZ comment fields can be parsed into per-conformer string properties by passing
a template. Placeholder names beginning with `_` are matched but omitted from
the stored properties, which is useful for fields such as `{_framenumber}`.
Energy operations convert the selected property to a number:

```python
def load_group():
    return Conformer_Group.from_multi_xyz(
        "conformers.xyz",
        comment_template="frame = {frame} energy = {energy}",
    )

group = load_group()
print(group.records()[0].properties["frame"])

# Match the frame field without retaining it as a property.
energy_only_group = Conformer_Group.from_multi_xyz(
    "conformers.xyz",
    comment_template="frame = {_framenumber} energy = {energy}",
)
print(energy_only_group.records()[0].properties)  # {"energy": "..."}

group = load_group()
group.sort_by_energy()                         # ascending; stable for ties

group = load_group()
group.filter_by_maximum_energy(10.0)           # inclusive upper limit

group = load_group()
group.retain_lowest_energy_percent(25.0)       # ceil(N * 25 / 100) records

group = load_group()
group.filter_by_boltzmann_population_ratio(
    minimum_ratio=0.01,
    temperature_kelvin=298.15,
    energy_to_joules_per_mole=1000.0,          # default assumes kJ/mol
)
```

Each filtering method changes the group in place. The Boltzmann population ratio
is `exp(-(E - E_min) / RT)` and is therefore relative to the lowest-energy
conformer. The default conversion assumes energies are in kJ/mol; use
`energy_to_joules_per_mole=1.0` for J/mol. All methods accept
`energy_property="another_field"` when the template uses a different name.

### Ring Detection

Use the same `Conformer_Group` instance when you want to keep all conformers together and store ring metadata for the molecule:

```python
from conformer_toolkit_cpp import Conformer_Group

group = Conformer_Group.from_xyz_files([
    "benzene_conf_001.xyz",
    "benzene_conf_002.xyz",
])
group.detect_rings(charge=0)

print("molecule adjacency:")
for atom, neighbors in enumerate(group.adjacency_table()):
    print(atom, list(neighbors))

print(len(group), len(group.rings()))
for ring in group.rings():
    print("ring atoms:", list(ring.atoms))
    for atom_adjacency in ring.adjacency_list:
        print(atom_adjacency.atom, list(atom_adjacency.adjacent_atoms))
```

`detect_rings` infers chemistry from the first conformer with RDKit. Because all conformers in a `Conformer_Group` are expected to share connectivity, the inferred bonds from conformer 0 are also stored as the group-level `adjacency_table()`.

Each stored `Ring_Record` contains:

- `atoms`: RDKit's ordered atom indices for the ring.
- `adjacency_list`: one `Ring_Atom_Adjacency` per ring atom, where `atom` is the atom index and `adjacent_atoms` contains the previous and next ring neighbors followed by any directly bonded non-ring atoms.

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

`--allow-reorder` calls `Conformer_Group::remove_duplicates` with its cleanup option enabled, so atom indexing is normalized before duplicate removal.

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
python3 src/deduplicate_rdkit.py --tolerance 0.001 input.sdf
```

Write unique representatives to an SDF:

```bash
python3 src/deduplicate_rdkit.py --tolerance 0.001 --write-unique unique.sdf input.sdf
```

Multiple SDF files can be passed:

```bash
python3 src/deduplicate_rdkit.py --tolerance 0.001 batch_1.sdf batch_2.sdf
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

Benchmark the C++ command-line path against Python RDKit `GetBestRMS`:

```bash
python3 src/benchmark_rdkit_vs_cpp.py --repeat 100000
```

The C++ path uses RDKit C++ for chemistry inference, graph matching, and aligned RMSD. The Python comparison uses RDKit `rdMolAlign.GetBestRMS`, which also considers symmetry.

## Limitations

- XYZ files do not encode bonds, bond orders, formal charges, isotopes, or stereo annotations. RDKit reconstructs these from coordinates and the supplied total charge, so unusual, organometallic, highly distorted, or ambiguous structures can still be misassigned.
- Total charge cannot be inferred reliably from XYZ and must be passed with `--charge`.
- Stereo comparison is only as reliable as RDKit's inferred graph and 3D stereo assignment. Nearly planar or otherwise geometrically ambiguous centers may remain unspecified.
- `--allow-reorder` can become expensive for highly symmetric molecules. `--max-mappings` bounds RDKit's match enumeration.
- When multiple topology-compatible mappings exist, index cleanup chooses the lowest-RMSD mapping; duplicate classification subsequently requires stereochemistry to match.
- The C++ tools do not compare conformer energies.

## Source Debug Notes

Debugging notes are embedded directly in the source files as `//` comments for C++ and `#` comments for Python.

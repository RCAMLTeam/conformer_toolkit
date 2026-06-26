# Tutorials

This folder contains notebook tutorials for the Python-facing conformer toolkit functions.

Run notebooks from the repository root after building the optional pybind11 extension:

```bash
pixi run build
PYTHONPATH=src pixi run python -m notebook tutorials
```

If you open a notebook directly from `tutorials/`, its setup cell still locates the repository root and adds `src/` to `sys.path`.

## Notebooks

- `01_cpp_python_interface.ipynb`: load XYZ conformers with `Conformer_Group`, remove duplicates, normalize atom indices, write cleaned XYZ files, and detect rings.
- `02_knn_k_optimizer.ipynb`: build KNN graphs from an RMSD matrix, select the smallest connected K, and inspect diagnostics.
- `03_rdkit_sdf_deduplication.ipynb`: use the RDKit SDF helper functions and command-line script for graph-aware SDF deduplication.


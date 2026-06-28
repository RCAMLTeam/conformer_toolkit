"""Public Python API for spivey-conformer-toolkit."""

from ._native import (
    Conformer_Batch,
    Conformer_Group,
    Conformer_Record,
    DuplicateEntry,
    Molecule,
    Remove_Duplicates_Result,
    Ring_Atom_Adjacency,
    Ring_Record,
    UniqueEntry,
    Vec3,
)
from .validation import validate_rmsd_matrix

__all__ = [
    "Conformer_Batch",
    "Conformer_Group",
    "Conformer_Record",
    "DuplicateEntry",
    "Molecule",
    "Remove_Duplicates_Result",
    "Ring_Atom_Adjacency",
    "Ring_Record",
    "UniqueEntry",
    "Vec3",
    "validate_rmsd_matrix",
]
__version__ = "0.1.0"

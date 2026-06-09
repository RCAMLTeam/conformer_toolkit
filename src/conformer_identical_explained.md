# `conformer_identical.cpp` Debug Notes

This file compares two ordered XYZ conformers.

## Main Flow

1. Read both XYZ files.
2. Require equal atom count.
3. Require identical atom symbols in identical order.
4. Center both coordinate sets.
5. Compute optimal rigid-alignment RMSD.
6. Print `identical` if RMSD is within tolerance, otherwise `different`.

## Repeat Mode

`--repeat n` runs the RMSD comparison in-process `n` times and reports average timing.

This avoids measuring process startup overhead when benchmarking the RMSD core.

## Debug Points

- Use this tool first when validating RMSD math.
- A pure translation or rotation should give RMSD close to zero.
- A changed bond length should give a nonzero RMSD.
- This tool does not reorder atoms and does not use chemistry. It is only correct when atom indexing already matches.

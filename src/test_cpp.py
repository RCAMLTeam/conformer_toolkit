#!/usr/bin/env python3
"""Regression tests for the RDKit-backed XYZ C++ tools."""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path


def run(command: list[str], expected_code: int, contains: str) -> None:
    result = subprocess.run(command, text=True, capture_output=True)
    if result.returncode != expected_code:
        raise AssertionError(
            f"expected exit {expected_code}, got {result.returncode}\n"
            f"command: {' '.join(command)}\nstdout:\n{result.stdout}\nstderr:\n{result.stderr}"
        )
    if contains not in result.stdout:
        raise AssertionError(f"missing {contains!r} in output:\n{result.stdout}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--identical", type=Path, required=True)
    parser.add_argument("--deduplicate", type=Path, required=True)
    args = parser.parse_args()

    data = Path(__file__).resolve().parent / "testdata"
    run(
        [str(args.identical), str(data / "h2_a.xyz"), str(data / "h2_translated.xyz")],
        0,
        "identical",
    )
    run(
        [str(args.identical), str(data / "h2_a.xyz"), str(data / "h2_stretched.xyz")],
        1,
        "different",
    )
    run(
        [
            str(args.deduplicate),
            "--allow-reorder",
            str(data / "benzene_ordered.xyz"),
            str(data / "benzene_reordered.xyz"),
        ],
        0,
        "unique_count 1",
    )
    run(
        [
            str(args.deduplicate),
            "--allow-reorder",
            str(data / "chiral_bromochlorofluoromethane_r.xyz"),
            str(data / "chiral_bromochlorofluoromethane_r_reordered.xyz"),
        ],
        0,
        "unique_count 1",
    )
    run(
        [
            str(args.deduplicate),
            "--allow-reorder",
            str(data / "chiral_bromochlorofluoromethane_r.xyz"),
            str(data / "chiral_bromochlorofluoromethane_s.xyz"),
        ],
        0,
        "unique_count 2",
    )
    run(
        [
            str(args.deduplicate),
            "--allow-reorder",
            str(data / "dichloroethene_e.xyz"),
            str(data / "dichloroethene_z.xyz"),
        ],
        0,
        "unique_count 2",
    )
    run(
        [
            str(args.identical),
            str(data / "chiral_bromochlorofluoromethane_r.xyz"),
            str(data / "chiral_bromochlorofluoromethane_s.xyz"),
        ],
        1,
        "reason chemistry_or_stereochemistry",
    )
    print("all C++ regression tests passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

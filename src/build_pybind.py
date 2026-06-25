#!/usr/bin/env python3
"""Configure and build the RDKit-backed C++ tools and pybind11 extension."""

from __future__ import annotations

import argparse
import os
import subprocess
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--build-dir", type=Path, default=Path("build"))
    parser.add_argument("--rdkit-root", type=Path, default=os.environ.get("RDKIT_ROOT"))
    parser.add_argument("--rdkit-source-root", type=Path, default=os.environ.get("RDKIT_SOURCE_ROOT"))
    args = parser.parse_args()

    root = Path(__file__).resolve().parents[1]
    build_dir = args.build_dir.resolve()
    configure = [
        "cmake",
        "-S",
        str(root),
        "-B",
        str(build_dir),
        "-DCMAKE_BUILD_TYPE=Release",
        "-DBUILD_PYTHON_BINDINGS=ON",
    ]
    if args.rdkit_root:
        configure.append(f"-DRDKIT_ROOT={args.rdkit_root.resolve()}")
    if args.rdkit_source_root:
        configure.append(f"-DRDKIT_SOURCE_ROOT={args.rdkit_source_root.resolve()}")

    subprocess.run(configure, check=True)
    subprocess.run(
        ["cmake", "--build", str(build_dir), "--parallel", "--target", "conformer_toolkit_cpp"],
        check=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

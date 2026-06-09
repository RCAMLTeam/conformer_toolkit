#!/usr/bin/env python3
"""Build the pybind11 extension for the C++ conformer batch API."""

from __future__ import annotations

import shlex
import subprocess
import sys
import sysconfig
from pathlib import Path

import pybind11


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    src = root / "src" / "conformer_toolkit_bindings.cpp"
    suffix = sysconfig.get_config_var("EXT_SUFFIX")
    if not suffix:
        raise RuntimeError("Python did not report an extension suffix")

    output = root / "src" / f"conformer_toolkit_cpp{suffix}"
    python_include = sysconfig.get_paths()["include"]

    # Keep this command explicit so build failures are easy to copy, inspect,
    # and rerun outside the helper.
    cmd = [
        "g++",
        "-O3",
        "-std=c++17",
        "-Wall",
        "-Wextra",
        "-shared",
        "-fPIC",
        f"-I{python_include}",
        f"-I{pybind11.get_include()}",
        str(src),
        "-o",
        str(output),
    ]

    print(" ".join(shlex.quote(part) for part in cmd))
    subprocess.run(cmd, check=True)
    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

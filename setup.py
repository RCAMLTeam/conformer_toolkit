"""Setuptools bridge for building the CMake-based Python extension."""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext


class CMakeExtension(Extension):
    def __init__(self, name: str) -> None:
        super().__init__(name, sources=[])


class CMakeBuild(build_ext):
    def build_extension(self, ext: Extension) -> None:
        # Prevent a legacy top-level extension from leaking into wheels when a
        # developer rebuilds without first deleting setuptools' build tree.
        build_lib = Path(self.build_lib)
        for pattern in ("conformer_toolkit_cpp*.so", "conformer_toolkit_cpp*.pyd"):
            for legacy_extension in build_lib.glob(pattern):
                legacy_extension.unlink()

        source_dir = Path(__file__).parent.resolve()
        output_dir = Path(self.get_ext_fullpath(ext.name)).parent.resolve()
        build_dir = Path(self.build_temp) / ext.name
        build_dir.mkdir(parents=True, exist_ok=True)

        rdkit_root = os.environ.get("RDKIT_ROOT", sys.prefix)
        configure_args = [
            "cmake",
            "-S",
            str(source_dir),
            "-B",
            str(build_dir),
            f"-DCMAKE_BUILD_TYPE={'Debug' if self.debug else 'Release'}",
            "-DBUILD_PYTHON_BINDINGS=ON",
            f"-DPython3_EXECUTABLE={sys.executable}",
            f"-DRDKIT_ROOT={rdkit_root}",
            f"-DCONFORMER_PYTHON_OUTPUT_DIRECTORY={output_dir}",
        ]
        rdkit_source_root = os.environ.get("RDKIT_SOURCE_ROOT")
        if rdkit_source_root:
            configure_args.append(f"-DRDKIT_SOURCE_ROOT={rdkit_source_root}")

        subprocess.run(configure_args, check=True)
        subprocess.run(
            [
                "cmake",
                "--build",
                str(build_dir),
                "--config",
                "Debug" if self.debug else "Release",
                "--target",
                "conformer_toolkit_cpp",
                "--parallel",
            ],
            check=True,
        )


setup(
    ext_modules=[CMakeExtension("conformer_toolkit._native")],
    cmdclass={"build_ext": CMakeBuild},
)

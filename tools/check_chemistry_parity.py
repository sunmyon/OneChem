#!/usr/bin/env python3
"""Compare shared runtime sources with an explicit AREPO checkout (read only)."""

import argparse
from pathlib import Path
import sys


CORE_DIRS = ("cooling", "include", "reaction_network", "solver", "tables", "thermodynamics")


def core_files(root):
    files = {p.relative_to(root) for p in root.iterdir()
             if p.is_file() and p.suffix in (".c", ".cc", ".cpp", ".h", ".hpp")}
    for directory in CORE_DIRS:
        files.update(p.relative_to(root) for p in (root / directory).rglob("*")
                     if p.is_file() and p.suffix in (".c", ".cc", ".cpp", ".h", ".hpp"))
    return files


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("arepo", type=Path, help="AREPO checkout root")
    args = parser.parse_args()
    local = Path(__file__).resolve().parents[1] / "src"
    remote = args.arepo.resolve() / "src"
    pairs = ((local / "chemistry", remote / "chemistry_ONECHEM", core_files),
             (local / "radiation", remote / "radiation",
              lambda root: {p.relative_to(root) for p in root.rglob("*")
                            if p.is_file() and p.suffix in (".c", ".cc", ".cpp", ".h", ".hpp")}))
    differences = 0
    checked = 0
    for left, right, collect in pairs:
        if not left.is_dir() or not right.is_dir():
            parser.error(f"missing source directory: {left} or {right}")
        for path in sorted(collect(left) | collect(right)):
            a, b = left / path, right / path
            checked += 1
            if not a.is_file() or not b.is_file():
                print(f"MISSING: {a if not a.is_file() else b}")
                differences += 1
            elif a.read_bytes() != b.read_bytes():
                print(f"DIFF: {left.name}/{path}")
                differences += 1
    print(f"Shared chemistry/radiation: {checked} files, {differences} differences")
    return int(differences != 0)


if __name__ == "__main__":
    sys.exit(main())

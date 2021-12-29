#!/usr/bin/env python3
from pathlib import Path
import sys

SEPS = (
    "\n",
    "\t",
    ":",
    ";",
    "=",
)


def main() -> None:
    files = sys.argv[1:]
    if len(files) != 2:
        print(f"too many files: {files}")
        exit(1)
    f1 = Path(files[0])
    f2 = Path(files[1])
    lines1 = decomp(f1.read_text())
    lines2 = decomp(f2.read_text())
    list_cmp(lines1, lines2)
    if len(lines1) != len(lines2):
        print(f"length mismatch: {f1} {len(lines1)} vs {f2} {len(lines2)}")
        exit(1)

    print(f"{f1} and {f2} match!")


def decomp(val: str):
    for sep in SEPS:
        if sep in val:
            return val.split(sep)
    raise ValueError(f"Cannot decompose {val}")


def list_cmp(l1: list[str], l2: list[str]):
    if len(l1) != len(l2):
        raise ValueError(f"Length mismatch: {l1} vs {l2}")

    for i in range(len(l1)):
        if l1[i] == l2[i]:
            continue
        else:
            try:
                f1 = float(l1[i])
                f2 = float(l2[i])
            except ValueError:
                f1 = None
                f2 = None
            if f1 is not None and f2 is not None and f1 == f2:
                continue
        try:
            list_cmp(decomp(l1[i]), decomp(l2[i]))
        except ValueError as e:
            print(e)
            raise ValueError(f"parse error on idx {i}\n{l1[i]}\n{l2[i]}")


###

if __name__ == "__main__":
    main()

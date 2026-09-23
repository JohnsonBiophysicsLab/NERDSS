#!/usr/bin/env python3
"""Compare two JSON restart files value by value.

    compare_restart_json.py <A.json> <B.json> [--ignore KEY]... [--abs-tol X]
                            [--max-report N]

Walks both documents together and prints the path of every value that differs
(the first --max-report of them, default 20) followed by a count, and exits
non-zero if anything differed.  --ignore KEY leaves out every value stored
under that key at any depth: the continuation check ignores currSimTime, which
a restart reaches through a different but equivalent expression, and the format
round trip ignores stateChangeIface, which the .dat format never carried.

Numbers are compared exactly unless --abs-tol is given, which lets two numbers
within that distance of each other pass: the format round trip allows the
1e-20 that the .dat format's twenty decimal places keep.  A JSON null stands
for a NaN (see maybe_null() in write_restart.cpp) and equals another null.
"""
import argparse
import json
import sys


def is_number(x):
    # bool is an int in Python; keep True distinct from 1 as JSON does.
    return isinstance(x, (int, float)) and not isinstance(x, bool)


def walk(a, b, path, ignore, abs_tol, diffs):
    if isinstance(a, dict) and isinstance(b, dict):
        for key in sorted(set(a) | set(b)):
            if key in ignore:
                continue
            sub = f"{path}.{key}" if path else key
            if key not in a or key not in b:
                diffs.append(f"{sub}: only in {'B' if key not in a else 'A'}")
            else:
                walk(a[key], b[key], sub, ignore, abs_tol, diffs)
    elif isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            diffs.append(f"{path}: length {len(a)} vs {len(b)}")
            return
        for i, (x, y) in enumerate(zip(a, b)):
            walk(x, y, f"{path}[{i}]", ignore, abs_tol, diffs)
    elif is_number(a) and is_number(b):
        if a != b and abs(a - b) > abs_tol:
            diffs.append(f"{path}: {a!r} vs {b!r}")
    else:
        if isinstance(a, bool) != isinstance(b, bool) or a != b:
            diffs.append(f"{path}: {a!r} vs {b!r}")


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[1])
    parser.add_argument("a")
    parser.add_argument("b")
    parser.add_argument("--ignore", action="append", default=[], metavar="KEY")
    parser.add_argument("--abs-tol", type=float, default=0.0, metavar="X")
    parser.add_argument("--max-report", type=int, default=20)
    args = parser.parse_args()

    with open(args.a) as f:
        a = json.load(f)
    with open(args.b) as f:
        b = json.load(f)

    diffs = []
    walk(a, b, "", set(args.ignore), args.abs_tol, diffs)
    for line in diffs[: args.max_report]:
        print(line)
    if diffs:
        print(f"{len(diffs)} value(s) differ")
        return 1
    notes = []
    if args.ignore:
        notes.append(f"ignoring {', '.join(args.ignore)}")
    if args.abs_tol:
        notes.append(f"numbers within {args.abs_tol:g}")
    print("identical" + (f" ({'; '.join(notes)})" if notes else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main())

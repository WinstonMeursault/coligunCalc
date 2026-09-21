#!/usr/bin/env python3
"""Verify the installed public header split for CPU-only and CUDA packages."""

from pathlib import Path
import argparse
import sys


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("prefix", type=Path)
    parser.add_argument("--cuda", action="store_true")
    args = parser.parse_args()

    include = args.prefix / "include/coilgun"
    required = [include / "coilgun.hpp"]
    optimization = include / "optimization"
    required.extend(path for path in optimization.glob("*.hpp")
                    if path.name != "cuda_batch_evaluator.hpp")
    missing = [str(path.relative_to(args.prefix)) for path in required if not path.is_file()]
    cuda_header = optimization / "cuda_batch_evaluator.hpp"
    if args.cuda and not cuda_header.is_file():
        missing.append(str(cuda_header.relative_to(args.prefix)))
    if not args.cuda and cuda_header.exists():
        print(f"unexpected CUDA-only header in CPU install: {cuda_header}")
        return 1
    if missing:
        print("missing installed headers:")
        print("\n".join(f"- {item}" for item in missing))
        return 1
    mode = "CUDA" if args.cuda else "CPU-only"
    print(f"install header split: PASS ({mode})")
    return 0


if __name__ == "__main__":
    sys.exit(main())

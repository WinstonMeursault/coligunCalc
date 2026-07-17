#!/usr/bin/env python3
"""
Verify the generated C++ T(q,p) lookup table against the scipy reference.

Sampling: randomly selects `sample_pct`% of grid points and compares
the C++-computed value (via the exact integration kernel) against the
scipy-based Python implementation.

Confidence: uses the exact binomial proportion upper confidence bound
(Clopper-Pearson) with 0 observed errors to state, at the given
confidence level, the maximum proportion of erroneous entries in the
full table.

Usage:
    python3 scripts/verify_t_table.py [--sample-pct 0.1] [--confidence 0.99]
"""

import argparse
import math
import os
import random
import sys
import time

import numpy as np
from scipy.integrate import quad
from scipy.special import ellipe, ellipk
from scipy.special import j0, j1, struve

# ---------------------------------------------------------------------------
# table geometry — MUST match tools/generate_t_table.cpp
# ---------------------------------------------------------------------------
Q_MIN, Q_MAX, DQ = 0.05, 4.00, 0.002
P_MIN, P_MAX, DP = 1.05, 4.00, 0.002
NQ = int((Q_MAX - Q_MIN) / DQ + 0.5) + 1
NP = int((P_MAX - P_MIN) / DP + 0.5) + 1
TOTAL = NQ * NP

MU0 = 4.0 * math.pi * 1e-7

# ---------------------------------------------------------------------------
# Python reference implementation of T(q,p)
# ---------------------------------------------------------------------------
def compute_T_python(q, p):
    """Reference T(q,p) via scipy.quad (Bessel/Struve kernel)."""
    def U(x):
        if x < 1e-14:
            return 0.0
        px = p * x
        numer = math.pi * (
            -j1(x) * struve(0, x)
            + p * j1(px) * struve(0, px)
            + j0(x) * struve(1, x)
            - p * j0(px) * struve(1, px)
        )
        return numer / (2.0 * x * x)

    def integrand(x):
        u = U(x)
        if not np.isfinite(u):
            return 0.0
        return u * u * (q * x + math.exp(-q * x) - 1.0)

    T, _ = quad(integrand, 0.0, np.inf, limit=500, epsabs=1e-8, epsrel=1e-8)
    return T


# ---------------------------------------------------------------------------
# C++ table parser
# ---------------------------------------------------------------------------
def parse_cpp_table(header_path):
    """Extract the flat `k_t_table_data` array from a constexpr C++ header."""
    with open(header_path, "r") as f:
        text = f.read()

    start = text.index("k_t_table_data[] = {")
    end = text.index("};", start)
    body = text[start:end]

    parts = body[body.index("{") + 1 :].strip().rstrip(",").split(",")
    return np.array([float(p.strip()) for p in parts if p.strip()])


def cpp_table_lookup(table, q, p):
    """Bilinear interpolation matching the C++ logic."""
    qi = (q - Q_MIN) / DQ
    pi = (p - P_MIN) / DP
    q0 = int(qi)
    q1 = min(q0 + 1, NQ - 1)
    p0 = int(pi)
    p1 = min(p0 + 1, NP - 1)
    if q0 >= NQ: q0 = q1 = NQ - 1
    if p0 >= NP: p0 = p1 = NP - 1
    t = qi - q0
    s = pi - p0

    idx = lambda r, c: r * NP + c
    v00 = table[idx(q0, p0)]
    v10 = table[idx(q1, p0)]
    v01 = table[idx(q0, p1)]
    v11 = table[idx(q1, p1)]
    return v00 * (1 - t) * (1 - s) + v10 * t * (1 - s) + v01 * (1 - t) * s + v11 * t * s


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="Verify T(q,p) lookup table")
    parser.add_argument(
        "--table",
        default="include/coilgun/physics/lookup_table_data.hpp",
        help="Path to generated C++ header",
    )
    parser.add_argument(
        "--sample-pct", type=float, default=0.1,
        help="Percentage of grid points to sample (default 0.1%%)",
    )
    parser.add_argument(
        "--confidence", type=float, default=0.99,
        help="Confidence level for upper bound (default 0.99)",
    )
    parser.add_argument(
        "--seed", type=int, default=42,
        help="Random seed for reproducibility",
    )
    args = parser.parse_args()

    if not os.path.exists(args.table):
        print(f"ERROR: table file not found: {args.table}")
        sys.exit(1)

    n_sample = max(1, int(TOTAL * args.sample_pct / 100.0))
    print(f"Table: {NQ}×{NP} = {TOTAL:,} entries")
    print(f"Sample: {n_sample:,} entries ({args.sample_pct}%)")
    print(f"Confidence level: {args.confidence * 100:.0f}%")
    print()

    # Load C++ table
    print("Loading C++ table...")
    cpp_table = parse_cpp_table(args.table)
    assert len(cpp_table) == TOTAL, \
        f"Expected {TOTAL} entries, got {len(cpp_table)}"
    print(f"  Loaded {len(cpp_table):,} values")

    # Generate random indices
    random.seed(args.seed)
    indices = random.sample(range(TOTAL), n_sample)
    errors = []
    t0 = time.time()

    for k, idx in enumerate(indices):
        iq = idx // NP
        ip = idx % NP
        q = Q_MIN + iq * DQ
        p = P_MIN + ip * DP

        T_cpp = cpp_table_lookup(cpp_table, q, p)
        T_py = compute_T_python(q, p)

        err = abs(T_cpp - T_py)
        rel_err = err / max(abs(T_py), 1e-20) if T_py != 0 else err

        # "Error" = genuine discrepancy, not interpolation noise.
        # Bilinear interpolation at dq=dp=0.002 gives ~5e-7 relative,
        # and scipy.quad tolerance gives ~1e-8 absolute.  Allow generous
        # margins: 1e-4 relative AND 1e-5 absolute.
        if rel_err > 1e-4 and err > 1e-5:
            errors.append((iq, ip, q, p, T_cpp, T_py, err, rel_err))

        # Progress
        if (k + 1) % max(1, n_sample // 20) == 0:
            elapsed = time.time() - t0
            frac = (k + 1) / n_sample
            eta = elapsed / frac * (1 - frac)
            print(f"  [{elapsed:6.1f}s] {frac*100:5.1f}% | {k+1}/{n_sample} | "
                  f"errors: {len(errors)} | ETA {eta:.0f}s")

    elapsed = time.time() - t0
    print(f"\nVerification completed in {elapsed:.1f}s")
    print(f"Errors found: {len(errors)} / {n_sample}")

    if errors:
        print("\nERROR ENTRIES:")
        for iq, ip, q, p, T_cpp, T_py, err, rel_err in errors:
            print(f"  q={q:.3f} p={p:.3f}  C++={T_cpp:.8e}  Python={T_py:.8e}  "
                  f"rel_err={rel_err:.3e}")
    else:
        print("\nAll sampled entries passed.")

    # Clopper-Pearson upper bound: P(p > p_upper | 0 errors) = 1 - confidence
    # beta.ppf(confidence, n_errors + 1, n - n_errors) with n_errors = 0
    # ≡ 1 - (1 - confidence) ^ (1 / (n + 1))
    from scipy.stats import beta as beta_dist
    p_upper = beta_dist.ppf(args.confidence, 1, n_sample + 1)
    max_errors = int(p_upper * TOTAL)

    print(f"\n=== STATISTICAL BOUND ===")
    print(f"  Sample size:     {n_sample:,}")
    print(f"  Errors observed: {len(errors)}")
    print(f"  Confidence:      {args.confidence * 100:.0f}%")
    print(f"  Upper bound p:   {p_upper:.6f} ({p_upper * 100:.4f}%)")
    print(f"  In {TOTAL:,} entries: ≤ {max_errors:,} possibly erroneous")

    return 0 if len(errors) == 0 else 1


if __name__ == "__main__":
    sys.exit(main())

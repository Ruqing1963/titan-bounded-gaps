#!/usr/bin/env python3
"""
generate_primes.py — Compute primes of Q_q(n) = n^q - (n-1)^q for 18 exponents.

This script reproduces the full dataset used in:
  "Bounded Gaps for High-Degree Polynomial Primes: A Cyclotomic Maynard-Tao Roadmap"
  Ruqing Chen (2025)

Output: data/prime_data.pkl  (~200 MB, dict mapping q -> sorted list of primes)
Runtime: approximately 2-4 hours on a modern machine with gmpy2.

Requirements:
    pip install gmpy2 sympy numpy
"""

import pickle, time, os, sys
from math import gcd

try:
    from gmpy2 import is_prime
    print("Using gmpy2 for primality testing (fast).")
except ImportError:
    from sympy import isprime as is_prime
    print("WARNING: gmpy2 not found, falling back to sympy (much slower).")
    print("Install gmpy2 for ~10x speedup: pip install gmpy2")

# Exponents (all prime)
Q_VALUES = [3, 5, 7, 11, 13, 17, 19, 23, 31, 37, 41, 43, 47, 53, 61, 71, 83, 167]
N_MAX = 10**8  # Search range for n


def titan_primes(q, n_max):
    """Find all primes of Q_q(n) = n^q - (n-1)^q for n = 2, ..., n_max."""
    primes = []
    count = 0
    for n in range(2, n_max + 1):
        val = n**q - (n - 1)**q
        if val > n_max:
            # Once Q_q(n) > n_max, check all remaining
            pass
        if is_prime(int(val)):
            primes.append(int(val))
            count += 1
        if n % 10_000_000 == 0:
            print(f"    q={q}: n={n:,}/{n_max:,}, found {count:,} primes so far")
    return primes


def main():
    output_path = os.path.join(os.path.dirname(__file__), '..', 'data', 'prime_data.pkl')
    results = {}

    print(f"Generating Titan polynomial primes for {len(Q_VALUES)} exponents, N <= {N_MAX:,}")
    print("=" * 70)

    t_total = time.time()
    for q in Q_VALUES:
        t0 = time.time()
        print(f"\n  q = {q} (degree {q-1})...")

        # For large q, Q_q(n) grows very fast, so we need fewer n values
        # Q_q(n) ~ q * n^(q-1), so Q_q(n) <= N_MAX implies n <= (N_MAX/q)^(1/(q-1))
        n_limit = min(N_MAX, int((N_MAX / q) ** (1.0 / (q - 1))) + 100)
        primes = titan_primes(q, n_limit)

        dt = time.time() - t0
        results[q] = sorted(primes)
        print(f"    Done: {len(primes):,} primes in {dt:.1f}s")

    dt_total = time.time() - t_total
    print(f"\n{'=' * 70}")
    print(f"Total: {sum(len(v) for v in results.values()):,} primes in {dt_total:.0f}s")

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, 'wb') as f:
        pickle.dump(results, f)
    print(f"Saved to {output_path} ({os.path.getsize(output_path) / 1e6:.1f} MB)")


if __name__ == '__main__':
    main()

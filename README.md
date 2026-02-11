# Bounded Gaps for High-Degree Polynomial Primes

**A Cyclotomic Maynard–Tao Roadmap**

Ruqing Chen · GUT Geoservice Inc., Montreal  
Paper II in the Titan Polynomial series

---

## Overview

This repository contains the paper, data, and reproducible code for:

> **Bounded Gaps for High-Degree Polynomial Primes: A Cyclotomic Maynard–Tao Roadmap**

We propose a concrete roadmap for establishing bounded gaps among primes of the **Titan polynomial family**:

$$Q_q(n) = n^q - (n-1)^q$$

Leveraging the *Arithmetic Shielding* property proved in [Paper I](https://doi.org/10.5281/zenodo.18582880), we show that the Bombieri–Vinogradov error term vanishes identically on a density-1 set of moduli. The bounded-gaps problem reduces to a single exponential-sum estimate on a sparse, algebraically structured set.

### Key Results

| Result | Status | Data |
|--------|--------|------|
| **Null-Sparse Decomposition**: BV error = 0 on M_null | Proved (Theorem 1) | 83× error concentration verified |
| **Massive Admissibility**: k_max = 2,173 for q = 167 | Proved (Theorem 2) | 345/345 null primes below p₁ |
| **p^½ Cancellation**: exponential sum bounds | Heuristic (Hypothesis 1) | 40/40 tests, max ratio 1.93 |
| **Conditional Bounded Gaps** | Theorem 3 (conditional on Sparse BV_q) | k·M_k ≈ 16,689 ≫ 4 |

## Repository Structure

```
titan-bounded-gaps/
├── README.md
├── LICENSE
├── .gitignore
├── paper/
│   ├── bounded_gaps_FINAL.pdf      # Main paper (8 pages)
│   └── bounded_gaps_FINAL.tex      # LaTeX source
├── figures/
│   ├── fig1_BV_error.pdf           # Figure 1: BV error concentration
│   ├── fig2_admissibility.pdf      # Figure 2: Admissibility advantage
│   ├── fig3_expsum.pdf             # Figure 3: Exponential sum cancellation
│   ├── fig4_gaps.pdf               # Figure 4: Gap distribution (Poisson)
│   └── fig5_roadmap.pdf            # Figure 5: Clustering + roadmap
├── data/
│   ├── prime_data.pkl.gz           # Full dataset (~70 MB compressed)
│   └── prime_statistics.csv        # Summary statistics for all 18 exponents
└── scripts/
    ├── generate_primes.py          # Compute Titan primes (reproduces data/)
    └── generate_figures.py         # Generate all 5 figures (reproduces figures/)
```

## Quick Start

### View the paper

The compiled paper is at [`paper/bounded_gaps_FINAL.pdf`](paper/bounded_gaps_FINAL.pdf).

### Reproduce the figures

```bash
# Decompress the dataset
cd data && gunzip -k prime_data.pkl.gz && cd ..

# Generate all figures
pip install numpy matplotlib scipy sympy
python scripts/generate_figures.py
```

### Reproduce the dataset from scratch

```bash
pip install gmpy2 numpy sympy
python scripts/generate_primes.py    # ~2-4 hours
```

## Dataset

The full dataset contains **~44 million primes** across 18 exponents (q = 3, 5, 7, ..., 167), computed for n ≤ 10⁸.

| q | degree | π_Q(10⁸) | k_max | Sophie Germain |
|---|--------|----------|-------|----------------|
| 3 | 2 | 9,389,636 | 5 | yes |
| 47 | 46 | 1,083,547 | 237 | no |
| 83 | 82 | 347,353 | 85 | yes |
| 167 | 166 | 493,941 | **2,173** | no |

Full statistics for all 18 exponents are in [`data/prime_statistics.csv`](data/prime_statistics.csv).

## Related Work

- **Paper I** (Arithmetic Shielding): [Zenodo doi:10.5281/zenodo.18582880](https://doi.org/10.5281/zenodo.18582880)
- Zhang, Y. (2014). *Bounded gaps between primes*. Ann. Math. 179(3).
- Maynard, J. (2015). *Small gaps between primes*. Ann. Math. 181(1).

## Citation

```bibtex
@article{chen2025bounded,
  title   = {Bounded Gaps for High-Degree Polynomial Primes:
             A Cyclotomic Maynard--Tao Roadmap},
  author  = {Chen, Ruqing},
  year    = {2025},
  note    = {GitHub: https://github.com/Ruqing1963/titan-bounded-gaps}
}
```

## License

This project is licensed under the MIT License — see [LICENSE](LICENSE).

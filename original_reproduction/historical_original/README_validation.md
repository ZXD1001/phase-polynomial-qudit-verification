# Matrix validation for the intersection-block construction

This directory contains a standalone consistency check for the phase-polynomial examples in the manuscript:

- validate_phase_polynomial_blocks.py

## What the script checks

For each example, the script constructs independently:

1. the cyclic generator projectors \(P_v=d^{-1}\sum_{r=0}^{d-1}g_v^r\);
2. the algebraic block \(\prod_{v\in A}P_v\);
3. the branch-summed effect of the two-layer adaptive measurement;
4. the verification operator obtained by uniformly sampling the stated blocks.

It then checks the operator two-norm residual between items 2 and 3 and diagonalizes item 4 to recover the spectral gap.

## Requirements

- Python 3.10 or later
- NumPy

## Run

From this directory:

    python validate_phase_polynomial_blocks.py

Representative output:

    d=2: residual=2.763e-16, gap=0.333333333333
    d=3: residual=5.473e-16, gap=0.333333333333
    d=4: residual=1.732e-16, gap=0.333333333333
    d=6: residual=6.523e-16, gap=0.500000000000

The \(d=4\) and \(d=6\) cases test composite local dimensions. Floating-point residuals are numerical consistency checks and are not substitutes for the analytical proofs.

## Archival step before submission

Deposit this script, this README and the exact submitted manuscript version in a permanent public archive. Replace the manuscript’s Code availability placeholder with the resulting DOI or stable URL.


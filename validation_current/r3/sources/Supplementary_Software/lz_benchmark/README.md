> Version note: references to original Eqs. (38)/(39) identify the previous submitted manuscript. See the package-level MANUSCRIPT_MAPPING.md for this revision. Local audit paths are not required to run the software.

# Li–Zhu Z/X variants: d=6 oriented-path benchmark

This directory implements the uniform MUB, CMUB, SMUB and SCMUB protocols of Li–Zhu, Quantum **10**, 2011 (2026), §§4.1–4.3, Eqs. (30)–(31), and Appendix D.1. It is a small-system operator audit for

\[
|\psi_n\rangle=6^{-n/2}\sum_{\boldsymbol x\in\mathbb Z_6^n}
\exp\!\left[\frac{2\pi i}{6}\sum_{v=0}^{n-2}(2x_v^2+3x_v)x_{v+1}\right]|\boldsymbol x\rangle.
\]

The path edges are oriented left to right; the left variable is quadratic. We computed n=2,3,4. These finite-size spectra do not prove asymptotic scaling or an advantage over every Li–Zhu protocol.

## Reproduce

Requirements: Python 3.12 and NumPy 2.3.5 were used; SciPy is not required. From this directory:

```sh
python3 benchmark_lz.py --n 2 3 4
python3 verify_exports.py
```

Use the Python interpreter in your prepared environment. The independent check must print `"status": "PASS"`; no machine-specific interpreter path is required.

## Conventions

- Local dimensions are all 6. `Z` means the computational basis. `X` means the Fourier eigenbasis with columns `F[x,o]=exp(+2*pi*i*x*o/6)/sqrt(6)`; its bra therefore has the minus sign.
- Every possible terminal site is evaluated. Terminal numbers in protocol names are one based; array metadata are zero based. The other sites are listed in increasing order. The full state and the front-outcome arrays use C order, with the final index fastest.
- MUB uses every Z/X front-basis string with weight `1/2^(n−1)`. CMUB uses only the all-Z and all-X strings, each with weight 1/2. SMUB/SCMUB additionally average uniformly over all terminal sites. No basis-probability optimization is performed.
- All front measurements can be parallel; the terminal performs the normalized conditional-state projection. This is a two-layer implementation. The number of adaptive test templates is distinct from the number of physical conditional terminal bases.
- Target-impossible front outcomes are determined exactly in the cyclotomic integer ring and **rejected**. No identity is inserted on those branches. This is an explicit implementation completion of the normalized-conditional-state formula; the text around Eq. (30) does not prescribe a unique normalized state on a zero-probability branch.
- Two mutually unbiased Z/X bases exist at d=6. We neither assume a complete set of MUB nor implement 3MUB, SD/CSD, or all conceivable changes of local bases.

## Results

| n | MUB, last site terminal | Best MUB over tested terminal sites | CMUB, last site terminal | Best CMUB over tested terminal sites | SMUB | SCMUB |
|---|---:|---:|---:|---:|---:|---:|
| 2 | 0.146446609407 | 0.500000000000 | 0.146446609407 | 0.500000000000 | 0.444206225118 | 0.444206225118 |
| 3 | 0.062591798343 | 0.172318849222 | 0 (sixfold unit eigenvalue) | 0.146446609407 | 0.250000000000 | 0.166666666667 |
| 4 | 0.036939006785 | 0.139571515673 | 0 (28-fold unit eigenvalue) | 0 (all terminal choices) | 0.239275596483 | 0.196528262671 |

These are spectral gaps from full dense Hermitian diagonalization. “0” denotes numerical zero, supported by target-orthogonal unit-eigenvalue witnesses and residuals. The unrounded eigenvalues and raw gaps, including harmless negative values around 10^−15, are retained in JSON/NPZ.

The fixed-terminal two-test CMUB protocol fails on several tested instances, whereas SMUB and SCMUB pass the same uniqueness checks. This is **not** a failure of all Li–Zhu protocols. Already at n=2 a suitable terminal makes the uniform MUB/CMUB gap 1/2, so no blanket strict advantage should be claimed from this benchmark.

## Files

- `benchmark_lz.py`: exact-support branch construction, verification operators, full spectra and zero-gap witnesses.
- `verify_exports.py`: independent floating-amplitude contraction versus exact-support branches; operator-average, acceptance and witness checks.
- `inputs_and_environment.json`: complete family/basis/weight conventions, environment and script hash.
- `results_n*.json`: per-test rank/branch counts and all per-protocol spectral/acceptance diagnostics.
- `conditional_branches_n*.npz`: target vector, phase table, basis masks, terminal/front order, normalized conditional vectors, exact cyclotomic integer pairs, exact probability numerators and denominators.
- `operators_and_spectra_n*.npz`: full test/verification matrices at n≤3; full eigenvalue arrays for every n; target-orthogonal accepted witnesses where the gap vanishes. At n=4 the dense matrices were diagonalized but are not stored; they are reproducible from the branch export/script.
- `verification_report.json`: independent validation residuals and PASS status.
- `technical_evidence.md`: source mapping, construction, fair-comparison limits and detailed terminal-dependent values.

The dense implementation is a reference calculation. It does not imply that the Li–Zhu protocol needs an exponentially long state vector online: this path family has a constant-bond transfer representation, allowing the terminal conditional vector to be evaluated with O(n) fixed-d arithmetic. Local control/precision resources must be priced in the same model for both protocols.

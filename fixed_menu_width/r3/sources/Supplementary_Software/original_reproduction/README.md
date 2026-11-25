> Version note: references to original Eqs. (38)/(39) identify the previous submitted manuscript. See the package-level MANUSCRIPT_MAPPING.md for this revision. Local audit paths are not required to run the software.

# Original manuscript matrix reproduction

This package reproduces the finite examples in the revised manuscript **Eq. (47) and Eq. (48)** (formerly Eqs. (38) and (39)). The d=2 and d=3 cases are explicitly labelled regression instances. It is separate from the proposed six-level path-family benchmark and does not certify every manuscript theorem.

The author’s historical script and README are preserved byte-for-byte in `historical_original/`. They were actually executed for this delivery. `provenance.json` records their source paths and SHA-256 hashes, together with hashes of the referenced main text and SI. The manuscript files are not modified or copied into this package.

## Run

Tested environment: Python 3.12.14, NumPy 2.3.5, macOS 15.3, arm64, complex128. Detailed numerical-library configuration is recorded in `results/environment.json`.

With Python 3.12 and NumPy installed:

```sh
python3 reproduce.py
```

For a fresh environment, install the pinned dependency first:

```sh
python3 -m pip install -r requirements.txt
python3 reproduce.py
```

Use the Python interpreter in your prepared environment; no machine-specific interpreter path is required.

The runner first checks the preserved historical script against its recorded hash, runs it unchanged, then performs the separate implementation. It exits nonzero on a failed check and writes machine-readable results. To preserve the supplied result files during a new run:

```sh
python3 reproduce.py --input inputs.json --output-dir rerun_results
```

The preserved script alone remains runnable:

```sh
python3 historical_original/validate_phase_polynomial_blocks.py
```

## Exact inputs and conventions

`inputs.json` specifies every coefficient, exponent, vertex, color block, expected spectral gap and tolerance.

| Case | Polynomial modulo d | Blocks | Role |
|---|---|---|---|
| d=4, n=4 | x0²x2 + 2x1³x2x3² + x2x3 | {0,1}, {2}, {3} | Revised Main Eq. (47), formerly Eq. (38) |
| d=6, n=3 | 2x0²x2 + 3x1³x2 | {0,1}, {2} | Revised Main Eq. (48), formerly Eq. (39) |
| d=2, n=4 | Same integer expression as the d=4 row, evaluated modulo 2 | Same three blocks | Regression |
| d=3, n=4 | Same integer expression as the d=4 row, evaluated modulo 3 | Same three blocks | Regression |

Sampling is uniform over the listed blocks. Vertex 0 is the leftmost tensor factor and most significant base-d digit. The inverse Fourier convention is `F_minus[x,j] = exp(-2*pi*i*x*j/d)/sqrt(d)`. Polynomial values are reduced modulo d before taking phases; no modular inverses are used.

The independent implementation uses absolute tolerance **1e-12**. The unmodified historical script uses **1e-10**. Residuals are matrix operator two-norms unless explicitly labelled as a vector norm or spectral absolute error.

## What the two implementations check

The preserved historical script constructs shifts, conjugated generators, cyclic generator-projector averages, their block products and independently branch-summed accepting effects. It then compares accepting operators and diagonalizes the verification operator.

The new implementation does not import that code. It constructs ideal joint-sector projectors directly as D-conjugated tensor Fourier projectors. Separately, it obtains conditional phases through finite differences of the **complete polynomial function**, forms the full adaptive measurement matrix, and groups its branches into every joint outcome effect. This avoids reusing the historical per-monomial conditional-phase formula or cyclic projector construction.

For each listed block it checks:

- Exact modular factorization for every computational input.
- Orthonormality of each conditional local basis and of the full branch measurement.
- Completeness of the full outcome measurement.
- Every joint-sector effect against the corresponding directly conjugated ideal projector.
- Accepting-effect idempotence, Hermiticity and target acceptance.

It then compares the full verification spectrum, multiplicities, spectral gap and target acceptance with the character-sector prediction. These finite-dimensional checks support the implementation identities; the analytical proofs remain necessary for arbitrary n or d.

## Actual delivery results

All four cases passed. The largest error among all independent checks was approximately **1.999e-15**, below the 1e-12 tolerance.

| Case | Historical accepting residual | Independent accepting residual | Independent largest joint-sector residual | Gap |
|---|---:|---:|---:|---:|
| d=2 regression | 2.763e-16 | 4.756218e-16 | 4.756218e-16 | 1/3 |
| d=3 regression | 5.680e-16 | 9.425599e-16 | 1.047943e-15 | 1/3 |
| d=4 revised Eq. (47) | 1.732e-16 | 1.224647e-16 | 1.731912e-16 | 1/3 |
| d=6 revised Eq. (48) | 6.505e-16 | 8.329022e-16 | 1.097546e-15 | 1/2 |

Historical residuals above retain the script’s printed precision. Full-precision independent results and spectral multiplicities are in `results/reproduction_results.json`. The decimal residuals can change with numerical libraries and evaluation order; see `NUMERICAL_COMPARISON.md`.

## Files and limits

- `historical_original/`: unmodified historical script and its historical README. Any proposed archival action in that README remains historical advice, not an action performed here.
- `inputs.json`, `reproduce.py`, `requirements.txt`: portable input specification and execution entry point.
- `provenance.json`: source identities and hashes. Source manuscript files remain outside the package and are not needed to execute it.
- `results/`: actual historical stdout/stderr, environment, runner log and machine-readable checks.
- `NUMERICAL_COMPARISON.md`: manuscript numbers versus this run, interpretation and limitations.

This is a **local reproducibility package**, not a public repository or permanent archival deposit. It does not establish physical noise robustness, scalable hardware synthesis, novelty, comparison with all universal protocols, or correctness of every main/SI theorem. No software licence is invented for the author’s historical source; public distribution should use an author-approved licence and archive.

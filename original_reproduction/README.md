# Original manuscript matrix reproduction

This module checks the d=4 and d=6 examples in Main Methods Sec. 4.4, “Composite-dimensional checks”, currently **Eqs. (49) and (50)** (`eq:d4-example`, `eq:d6-example`), and SI Notes 6.4.4 and 13.1 (Table S5). The d=2 and d=3 cases are regression instances. References to Eqs. (38)/(39) in preserved scripts, input descriptions and old records identify the same d=4/d=6 examples in the previous submission. See [MANUSCRIPT_MAPPING.md](../MANUSCRIPT_MAPPING.md) for the complete mapping.

The historical script and README are preserved byte for byte in `historical_original/`. `provenance.json` records their source identities and hashes. Historical source paths identify provenance; the runner does not need those files outside this package. Each run executes the preserved script and a separate full-sector implementation.

## Run

For the complete manuscript checks, follow the three-command quick start in the [package README](../README.md). It installs all dependencies in an environment beside the package and writes this module's results under `author/original_reproduction/results/` in your new destination.

To run **only this module**, first create that environment, then, from the repository root, use:

```sh
../cp-venv/bin/python original_reproduction/reproduce.py --output-dir ../original-rerun
```

Windows PowerShell:

```powershell
..\cp-venv\Scripts\python.exe original_reproduction/reproduce.py --output-dir ..\original-rerun
```

Use a fresh output directory to preserve prior results. This module itself requires only NumPy 2.3.5, listed in its local `requirements.txt`; the root environment already provides it. The runner checks the historical script's hash, executes it unchanged, and runs the independent construction. It exits nonzero if a check fails. Read `reproduction_results.json` for `passed: true` and the per-case results; `environment.json` records the actual Python, NumPy, numerical libraries and `complex128` arithmetic. `historical_stdout.txt` preserves the original script's printed output.

The reference residuals below were recorded with Python 3.12.14 and NumPy 2.3.5 on macOS 15.3 (arm64). That is the recorded environment, not an operating-system requirement. Different numerical libraries can change the final digits while all checks remain within tolerance.

## Exact inputs and conventions


`inputs.json` specifies every coefficient, exponent, vertex, color block, expected spectral gap and tolerance.

| Case | Polynomial modulo d | Blocks | Role |
|---|---|---|---|
| d=4, n=4 | x0²x2 + 2x1³x2x3² + x2x3 | {0,1}, {2}, {3} | Main Eq. (49), `eq:d4-example`; historical Eq. (38) |
| d=6, n=3 | 2x0²x2 + 3x1³x2 | {0,1}, {2} | Main Eq. (50), `eq:d6-example`; historical Eq. (39) |
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

## Preserved reference results

The supplied `results/` records the 2026-09-28 reference execution. All four cases passed in that run. The largest error among all independent checks was approximately **1.999e-15**, below the 1e-12 tolerance.

| Case | Historical accepting residual | Independent accepting residual | Independent largest joint-sector residual | Gap |
|---|---:|---:|---:|---:|
| d=2 regression | 2.763e-16 | 4.756218e-16 | 4.756218e-16 | 1/3 |
| d=3 regression | 5.680e-16 | 9.425599e-16 | 1.047943e-15 | 1/3 |
| d=4, current Eq. (49) | 1.732e-16 | 1.224647e-16 | 1.731912e-16 | 1/3 |
| d=6, current Eq. (50) | 6.505e-16 | 8.329022e-16 | 1.097546e-15 | 1/2 |

Historical residuals above retain the script’s printed precision. Full-precision independent results and spectral multiplicities are in `results/reproduction_results.json`. The decimal residuals can change with numerical libraries and evaluation order; the preserved `NUMERICAL_COMPARISON.md` explains the historical comparison.

## Files and limits

- `historical_original/`: unmodified historical script and its historical README. Any proposed archival action in that README remains historical advice, not an action performed here.
- `inputs.json`, `reproduce.py`, `requirements.txt`: portable input specification and execution entry point.
- `provenance.json`: source identities and hashes. Source manuscript files remain outside the package and are not needed to execute it.
- `results/`: actual historical stdout/stderr, environment, runner log and machine-readable checks.
- `NUMERICAL_COMPARISON.md`: historical audit of the 2026-09-28 reproduction; its old equation numbers, page numbers and revision advice are retained as history.

This is a **local reproducibility package**, not a public repository or permanent archival deposit. It does not establish physical noise robustness, scalable hardware synthesis, novelty, comparison with all universal protocols, or correctness of every main/SI theorem. No software licence is invented for the author’s historical source; public distribution should use an author-approved licence and archive.

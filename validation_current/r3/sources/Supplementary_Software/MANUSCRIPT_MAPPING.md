# Manuscript-to-software mapping

This mapping applies to the 2026-09-28 revision of **Local verification of nonlinear phase-polynomial qudit states**. Equation, theorem and table numbers were checked against the revised LaTeX sources and their compiled label records. Paths below are relative to `Supplementary_Software/`; manuscript labels are included to make later renumbering traceable.

## Original examples and numbering

The original examples and the new directed path are different targets.

| Target | Previous submitted Main | Revised Main | SI | Input and generated evidence |
|---|---|---|---|---|
| d=4, n=4: x0²x2 + 2x1³x2x3² + x2x3 mod 4 | Eq. (38) | Eq. (47), `eq:d4-example` | Eq. (S68); Note 12.1; Table S4 | `original_reproduction/inputs.json`; `results/reproduction_results.json` under that directory, case `d4_eq38` |
| d=6, n=3: 2x0²x2 + 3x1³x2 mod 6 | Eq. (39) | Eq. (48), `eq:d6-example` | Eq. (S69); Note 12.1; Table S4 | The same input/result files, case `d6_eq39` |
| d=2,3 reductions of the d=4 integer polynomial | Additional regression cases, not separate numbered examples | Methods, numerical reproduction | Note 12.1; Table S4 | Cases `d2_regression`, `d3_regression` |
| New directed path: Σv(2xv²+3xv)xv+1 mod 6, n≥2 | Not the original Eq. (39) | Eq. (20), `eq:path-family` | Eq. (S87) | `structured_benchmark/`, `lz_benchmark/`, and `resources/` |

`original_reproduction/reproduce.py` runs the unmodified historical script and an independent full-sector implementation. Table S4 is produced from accepting-effect and all-sector residuals, expected/observed spectral multiplicities and gaps. Historical stdout/stderr and environment information are retained in `original_reproduction/results/`. Existing source metadata and historical files may retain the former equation numbers; their bytes have not been rewritten to conceal that provenance.

## New family and restricted comparisons

| Revised manuscript result | SI location | Code and artifacts | Evidence and scope |
|---|---|---|---|
| Eq. (20): target; Eq. (21): local menu; Theorem 2 and Eq. (22): two-color gap 1/2 and resources | Note 10.1; Eqs. (S87)–(S91) | `structured_benchmark/run_benchmark.py`, `benchmark_results.json` → `paths`; `resources/verify_resource_model.py` | The all-n statement is proved analytically. Physical conditional effects are cross-checked against conjugated projectors for n=2,3,4; full gap diagonalization is performed for n=2,3. |
| Endpoint rank and ordinary Z6 multiaffine comparison; local CRT factorization | Note 10.2; Eqs. (S92), (S93) | Analytical proof in SI; `benchmark_results.json` → `n2_direct_six_level`; CRT basis checks in `resources/resource_checks.json` | The all-n endpoint-rank/LU separation is an analytic statement, not inferred from finite enumerations. The two-site state is attributed to Appel et al. |
| Eq. (23): independent HT-type matching; shared-index comparison | Note 10.3; Eqs. (S94), (S95) | `run_benchmark.py` → `path_check`, `physical`; `benchmark_results.json` → `paths` | Direct physical matching versus product-average identities, saturation witnesses, and independent/shared operator differences at n=2,3,4. The equal 1/6 gaps for all n≥2 have separate analytic proofs. |
| Proposition 3: matching-only class optimum | Note 10.4; Eq. (S96) for n=2 | SI proof; `benchmark_results.json` → `paths[0].optimized_four_setting_gap` | Arbitrary independent sets, correlated preselected menu indices and weights; all unmatched terminals are ignored. Optimum 1/4 at n=2 and 1/6 at n≥3. This is not an optimum over all fixed local measurements. |
| Table 2: exact n=2 CRT product-Pauli optimum 3/11 | Note 10.5; Eqs. (S97)–(S104) | `structured_benchmark/exact_n2_certificate.py`; `benchmark_results.json` → `n2_exact_qutrit_certificate`, `n2_direct_six_level` | Exact rational arithmetic in Q(ω3), including support decisions and the primal/dual identities. Thirty settings attain the optimum in the full prescribed 144-setting class. |
| Table 2: numerical n=3 CRT comparison | Note 10.6; Eqs. (S105), (S106) | `structured_benchmark/pauli_search.py`; `benchmark_results.json` → `n3_numeric` | All 64 qutrit product-Pauli bases are considered, with CRT reduction covering 1728 six-level settings. Weights, complement basis, dual density matrix and convergence history are stored. These are floating-point primal/dual candidates, not outward-rounded certified bounds. |
| Table 2: Li–Zhu MUB, CMUB, SMUB, SCMUB rows | Note 11.1; Eqs. (S107)–(S109); Tables S2, S3 | `lz_benchmark/benchmark_lz.py`, `results_n2.json` through `results_n4.json`, `conditional_branches_n*.npz`, `operators_and_spectra_n*.npz` | All terminal choices; fixed Z/F front bases and uniform probabilities. Table S2 uses terminal numbers 1,…,n, whereas array metadata use 0,…,n−1. Numerical zero gaps retain witness vectors and residuals. |
| Independent check of exported Li–Zhu operators | Note 12.2 | `lz_benchmark/verify_exports.py`, `verification_report.json` | Independent floating contraction versus exact support, branch probabilities, target acceptance, spectra and witness checks. |

In the exact certificate, Eq. (S99) specifies the weights, Eqs. (S100) and (S101) the annihilating polynomial and target-kernel projector, Eq. (S102) the dual state, Eq. (S103) its sixteen rejection expectations, and Eq. (S104) the restricted optimum. The software verifies the explicit identities; the SI gives their implication for optimality, including arbitrary correlated setting weights. These tests do not optimize arbitrary local PVMs, POVMs or LOCC protocols.

## Computation, control and statistics

| Revised manuscript result | SI location | Code and artifacts | Interpretation |
|---|---|---|---|
| Eq. (24): left/right conditional-state contraction | Note 11.2; Eqs. (S111)–(S113) | `resources/verify_resource_model.py`; `resource_checks.json` → `cross_checks`, `large_n_transfer` | Exactly 7,368 n=2,3,4 branches agree in integer coefficients and probabilities. The n=500 result is one observed branch, not a spectrum or complete branch enumeration. O(n) counts fixed-dimensional ring operations; integer bit lengths grow. |
| Complete three-basis compilation and single-ray comparison | Note 11.3; Eq. (S114) | `resources/resource_checks.json`; `compiled_three_bases.json` | Explicit gate matrices, reconstruction residuals, at most 15 two-level rotations plus a diagonal layer for a full unitary, or at most five rotations for a ray. These are constructive upper bounds in the stated continuous-control model. |
| Gate error budget η≤2Gζ | Note 11.3; Eq. (S115), with Note 8 | Analytic argument in SI; gate counts in `resource_checks.json` | No simulation or certificate of an arbitrary hardware noise model, readout noise, correlated memory or fault tolerance is claimed. |
| Main Eq. (3) and numerical copy-count illustration | Eq. (S110); independent-copy assumptions in Main and SI | `summarize_results.py`; `comparison.csv`, `comparison.json` | Uses N=ceil[ln δ/ln(1−νε)] at ε=0.01, δ=0.05, with settings independently resampled on each copy. Rounded N ratios are not exactly gap ratios. Correlated-source retained-output guarantees are different tasks. |
| Data and Code Availability | Note 12 | `README.md`, `run_all.py`, `requirements.txt`, `validation_receipt.json` | The package is accompanying software, with an independently executed six-stage reproduction receipt; no public DOI or repository accession is asserted. |

## Portable reproduction and evidence preservation

From this directory, in an environment satisfying `requirements.txt`:

```sh
python3 run_all.py --destination ../software-rerun
```

The destination must be new. The six stages are original-example reproduction, Li–Zhu generation, independent export checks, structured matching/CRT checks, resource checks, and common statistical summarization. The runner copies source/input files and writes new logs/results only to that destination. It verifies that packaged source bytes remain unchanged.

For the exact certificate alone:

```sh
python3 structured_benchmark/exact_n2_certificate.py
```

Historical files in `original_reproduction/historical_original/` are preserved byte for byte, including their original README. Provenance, source receipts and recorded numerical outputs are not rewritten when documentation is updated. Consult `validation_receipt.json` for the actual independent rerun environment and per-stage exit status.

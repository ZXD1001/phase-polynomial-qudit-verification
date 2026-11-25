# Reproduction guide

For installation and the command to run all checks, see the [quick start](README.md#reproduce-the-results).

Choose a **new result directory outside this package**. If `revision-rerun` already exists, choose another name. The runner preserves the supplied sources and recorded results. Do not use Python `-O` or `-OO`, because the exact checks require assertions. You do not need to run the individual modules first.

Missing dependencies are reported before a result directory is created, together with an installation command for the selected Python interpreter. The main runner uses UTF-8 for its logs and all scientific child processes, including on Windows; paths may contain spaces or non-ASCII characters without requiring manual encoding settings.

## Check the result

A successful run ends with `PASS`. Open `RUN_RECEIPT.json` in your **new result directory** and check that `status`, `scientific_checks` and `import_closure` are all `PASS`. The receipt records the actual Python and library versions, any `dependency_deviations`, commands, exit codes, source hashes and result hashes. `source_bytes_unchanged: true` confirms that the supplied package was preserved. `IMPORT_CLOSURE.json` records the separate check that all scientific imports came from the package or its execution copy.

The main result files, relative to that new directory, are:

| Result | Where to look |
|---|---|
| Original examples: residuals, spectra, gaps and tolerances | `author/original_reproduction/results/reproduction_results.json`; environment and original-script output in the same folder |
| Common-model sample counts, including 598 versus 1796 | `author/comparison.csv` and `author/comparison.json` |
| Matching/CRT and uniform Li–Zhu checks | `author/structured_benchmark/benchmark_results.json` and `author/lz_benchmark/` |
| Exact optimized-comparison certificates and separate floating candidate | `r1/lz_n4_exact_upper_result.json`, `r1/lz_fullmenu_exact_half_result.json`, `r1/n4_primal_recomputed.json` |
| Finite width certificates and supporting checks | `r2/results/`, `r3/results/`, `r4/results/` and their execution receipts |
| Exact formula values for the width table | `tables/claim_values.json` and `tables/width_rows.tex` |

See [MANUSCRIPT_MAPPING.md](MANUSCRIPT_MAPPING.md) for the corresponding manuscript sections and the meaning of each check. Floating residuals can vary in their last digits across numerical libraries; compare the checks with their recorded tolerances, not by requiring identical decimal strings.

If a run fails, inspect the first nonzero `exit_code` in the receipt's `steps` and the associated files in `logs/`. A failed or interrupted run does not count as a completed reproduction. No fixed per-step time limit is imposed; a slower computer can continue until completion. You can interrupt with Ctrl+C and rerun with a new destination.

## Time and tested environments

On 25 November 2026 the scientific pipeline took about five minutes on macOS 15.3 (arm64) with Python 3.12.14 or 3.14.6 and the pinned dependencies. These timings exclude dependency installation. Native runs on Windows and Linux have not been verified. The calculations use a CPU. The n=100 and n=1000 table entries evaluate proved formulas rather than construct large many-body matrices.

## Contents


| Path | Purpose |
|---|---|
| `original_reproduction/`, `structured_benchmark/`, `lz_benchmark/`, `resources/` | Preserved original-example, matching/CRT, uniform Li–Zhu, and resource regressions |
| `run_all.py` | Original six-stage entry point, called by the new runner |
| `revision_benchmarks/` | R1 optimized comparisons, exact certificates, operator identity and ray count checks |
| `fixed_menu_width/r2/` | Tight finite fixed-library width certificates, including n=6 width-one companion |
| `fixed_menu_width/r3/` | Exact-width graph/local checks and odd-factor extension checks |
| `fixed_menu_width/r4/` | Near-optimal witness certificates, combinatorics and independent native local checks |
| `manuscript_tables/` | Exact formula evaluation and six manuscript rows; no large-system matrix simulation |
| `validation_history/` | Previous documentation and receipt retained with an explicit historical identity |
| `validation_current/` | Byte-preserved historical execution from 2026-09-30 |

Round-specific modules with the same name are deliberately isolated. They must not be merged by filename: different versions can have different hashes and assumptions. Unused R2 exploratory searches and unfinished n=7 proposals are excluded. The R3 generator's n=7,8,100 files are formula instances, not claims of exhaustive native checks at those sizes. The packaged historical author snapshot and R4 author archive support provenance checks only; they are not the current software distribution. Personal home paths in environment metadata are redacted and affected manifests refreshed, as disclosed in `METADATA_REDACTIONS.json`; scientific algorithms, certificates and inputs are unchanged.

## Scientific scope


The fixed-menu results concern the **complete** four-PVM library `{Z,B0,B1,B2}` for the native six-dimensional path and a two-layer protocol. The waiting set and first-layer bases are fixed before that copy's outcomes; all six-valued first-layer records and arbitrary correlated choices are allowed. The waiting set need not be independent. The second layer is parallel with no within-layer feed-forward; each site is measured once; tests are perfectly complete. Width counts sites held for the second layer, not controllers, peak memory or latency. Ideal gap loss is not measurement noise or state infidelity.

The n=3 optimized Li–Zhu construction yields the same two-colour operator. For n=4, the 32-setting Z/F-front comparison has an exact upper bound `197537/500000` and a separately marked floating feasible candidate near `0.3950735468`. The 23-setting exact half-gap certificate allows unrestricted terminal rays, so it is **outside** the class restricting every complete terminal PVM to the four-PVM library. The 23-setting support is not asserted minimal. The ray inventory describes this implementation only. Uniform comparison results remain regressions, not optimized method limits.

The tight finite fixed-library bounds are `2/5 <= nu_1*(4) <= 29/70`, `1/3 <= nu_1*(5) <= 26915/68048`, and `7/16 <= nu_2*(6) <= 10049/22500`. The n=6 width-one companion has raw maximum `1319/3000` and is necessary for covering **at most** two sites. Individual native checkers enumerate an **exact** waiting-set size. Raw witness bounds are retained separately from the gap bound capped by the independent menu ceiling `1/2`; for example the R4 n=4 width-two raw maximum is `10/19`, not `1/2`.

Exact all-length thresholds and near-optimal bounds are proved in the SI. Finite programs are implementation checks, not a substitute for these proofs. The odd-factor checks apply to native dimension `2d`, odd `d >= 3`, with a library of `d+1` PVMs; the six-dimensional R4 constants are not extrapolated to that family.

## Precision, evidence and publication status


Each run step is identified as `EXACT`, `NUMERIC`, `FORMULA_EVALUATION`, `PROVENANCE`, or explicitly mixed. Exact rational/cyclotomic certificates never use floating feasibility as proof. The n=3 CRT primal/dual values remain floating candidates. Numerical tolerances and certificate bytes are preserved; the added n=4 numerical re-evaluation uses `1e-12`. The table generator evaluates the stated formulas with `Fraction`; n=100/1000 rows are not dense-spectrum computations. The n=500 transfer example processes one record, not an entire protocol spectrum.

Tested runtime and dependency versions, commands, exit codes, source hashes and generated-file hashes are in that run’s receipt. Absolute temporary execution paths in raw receipts merely identify the actual run and are not required inputs.

## Execution records and import checks

The entry point runs the original six stages, six R1 checks (five preserved checks plus numerical re-evaluation of the supplied n=4 weights), R2's 13 stages, R3's 11 stages, R4's 21 stages, and exact table arithmetic. R3 repeats the historical author's six stages; this is regression redundancy, not additional independent evidence. R4's first stage verifies packaged historical provenance, while the top-level `AUTHOR` group executes the baseline in the new run.

Despite its name, `validation_current/` preserves the **2026-09-30 historical execution**, including its original receipt and import observations. It is not a receipt for the current package or for your run. Every fresh execution writes its own `RUN_RECEIPT.json` and `IMPORT_CLOSURE.json` to the destination you select. No separate author-side audit directory is needed. `SOURCE_PROVENANCE.json` retains the original source hashes and annotates recorded adjustments to execution wrappers. `REVIEWER_REPRODUCTION_UPDATE_20261010.json` separately records the later changes to execution wrappers, file output, documentation and tests, including their before/after hashes; those changes preserve the scientific algorithms, certificates, inputs and historical receipts. `SHA256SUMS.txt` covers the distributed files except itself.

The 2026-10-06 observer repair uses an explicit script launcher; it does not depend on Python loading `sitecustomize`. Every scientific Python child is registered before launch, explicitly loads its observer, and produces a completion record. The launcher propagates through the original runners to all grandchildren and deeper descendants without editing their scientific source. An independent Python subprocess audit hook detects direct launches that bypass the launcher. Empty observations, missing completion records, missing descendant branches, unregistered records or unexpected module paths cannot yield import-provenance `PASS`.

Coverage is checked against this entry point's fixed process topology: **11 top-level Python processes and 60 descendants (71 total)**. The eight non-leaf script paths and their expected child counts are declared in `run_revision_checks.py`; all other scripts must be leaves. If the pipeline is deliberately extended, update those declarations and rerun the tests. The supported scientific commands invoke an existing `.py` script through the selected `sys.executable`, without interpreter flags or `shell=True`; unsupported Python command forms fail rather than silently escape observation. The observer only instruments process launch and records module filenames at completion; it does not alter scientific functions, certificate inputs or numerical tolerances. It is a reproducibility check for the supplied pipeline, not a security sandbox for hostile programs.

The optional observer regression tests are separate from the scientific run. From the package directory, use the same environment:

```sh
../cp-venv/bin/python portability/test_import_audit.py
../cp-venv/bin/python portability/test_reviewer_portability.py
```

On Windows, replace `../cp-venv/bin/python` with `..\cp-venv\Scripts\python.exe`.

The import-audit tests use temporary scripts and check normal child/grandchild coverage, absent and partial observations, bypassed launching, abrupt termination, and preloaded or unexpected `sitecustomize` modules. The preserved `portability/sitecustomize.py` is the legacy observer and is no longer used by the official entry point. See `IMPORT_AUDIT_REPAIR_20261006.json` for that repair's scope. The reviewer-portability tests check the later execution-wrapper changes. These optional tests are not required before running the main entry point.


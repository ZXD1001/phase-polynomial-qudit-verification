# Fixed complete-library width checks

Keep the three round namespaces separate. `r2/reproduce.py` runs 13 checks; `r3/reproduce.py` runs 11 stages including a six-stage historical-author regression; `r4/reproduce.py` runs 21 checks. Each requires `--destination` pointing to a new directory outside its source folder. Use the top-level `run_revision_checks.py` to execute them together and record unified provenance.

Every mathematical source and certificate input is imported byte-for-byte; see `../SOURCE_PROVENANCE.json`. Generated results remain outside these source folders. The R2 code folder contains only the dependency closure needed for declared validation, not unused exploratory optimizers. Utility modules can contain dormant optimization functions, but the declared entry points do not invoke new optimization searches.

The mathematical model fixes all complete PVMs to `{Z,B0,B1,B2}` and allows all six-valued records, arbitrary correlations, non-independent waiting sets, and all perfectly complete acceptance rules in two layers with no within-layer feed-forward. Individual finite native certificates enumerate exact waiting-set sizes. R2 combines smaller-size certificates explicitly for at-most-width conclusions.

The R3 historical source snapshot and R4 `source_archives/` plus `reference_receipts/` are isolated historical provenance. A successful historical receipt is never reported as a new scientific execution. The top-level runner genuinely reruns the preserved baseline and every declared revision group.

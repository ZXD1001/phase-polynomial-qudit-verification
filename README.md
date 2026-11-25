# Local verification of nonlinear phase-polynomial qudit states

This repository contains the supplementary software for *Local verification of nonlinear phase-polynomial qudit states*. A single Python command runs the finite-dimensional examples, checks the supplied exact certificates, and regenerates the six rows of Table 2. The calculations use a CPU and run offline after the dependencies are installed.

## Reproduce the results

Open a terminal in the repository directory containing `run_revision_checks.py` and `requirements.txt`. Use Python 3.12 or 3.14 and an unused output directory. The commands below create both the virtual environment and the results directory beside the repository.

**macOS or Linux**

```sh
python3.12 -m venv ../cp-venv
../cp-venv/bin/python -m pip install -r requirements.txt
../cp-venv/bin/python run_revision_checks.py --destination ../revision-rerun
```

**Windows PowerShell**

```powershell
py -3.12 -m venv ..\cp-venv
..\cp-venv\Scripts\python.exe -m pip install -r requirements.txt
..\cp-venv\Scripts\python.exe run_revision_checks.py --destination ..\revision-rerun
```

For Python 3.14, substitute `python3.14` or `py -3.14` in the first command. Use normal Python execution with assertions enabled. The full run takes about five minutes on the tested macOS arm64 system. The pinned dependencies are NumPy 2.3.5, SciPy 1.17.0 and SymPy 1.14.0.

The final output line begins with **PASS**. In `../revision-rerun/RUN_RECEIPT.json`, the fields `status`, `scientific_checks` and `import_closure` should each read `PASS`.

## Locate a result

Paths in this table are relative to the result directory you chose.

| Paper result | Output |
| --- | --- |
| Original four- and six-dimensional examples, Eqs. (49) and (50) | `author/original_reproduction/results/reproduction_results.json` |
| All-pass sample counts, including 598 versus 1796 | `author/comparison.csv` |
| Four-site Z/F upper bound and exact 23-setting construction | `r1/lz_n4_exact_upper_result.json`, `r1/lz_fullmenu_exact_half_result.json` |
| Finite waiting-width certificates | `r2/results/`, `r3/results/`, `r4/results/` |
| Necessary and sufficient widths in Table 2 | `tables/claim_values.json`, `tables/width_rows.tex` |

The entries for 100 and 1000 sites evaluate the paper's proved formulas. For Supplementary Table S3, the four-site MUB gap at terminal 4 evaluates to `0.036939006784454276` on the reference system. At 12 decimal places this gives `0.036939006784`. The November 2026 manuscript prints `0.036939006785` in that cell.

The complete [paper-to-code index](MANUSCRIPT_MAPPING.md) lists the calculations, outputs and manuscript sections. The [reproduction guide](REPRODUCTION.md) describes the certificates, numerical tolerances, tested environments and optional checks.

## File integrity

`SHA256SUMS.txt` records the distributed file bytes. The files use their original line endings. If cloning with Git, use `git clone --config core.autocrlf=false <repository-url>` to preserve them. The release ZIP also preserves those bytes.

# Supplementary Software

For **Local verification of nonlinear phase-polynomial qudit states**.

This package accompanies the revised main text and Supplementary Information. It includes complete inputs, executable implementations and recorded numerical results. It contains no experimental data. The scripts require Python and NumPy only; no Zotero library, manuscript files, network service or external optimization solver is needed to execute them.

## Reproduce all calculations

Tested: Python 3.12.14, NumPy 2.3.5, macOS arm64. In a prepared environment:

```sh
python3 run_all.py --destination ../software-rerun
```

The destination must not already exist. The runner copies its source and input files there, performs all six steps, and records stdout, stderr and generated results. Supplied results remain untouched. Dependencies are pinned in `requirements.txt`; install them in your preferred Python environment if necessary. The runner uses the currently selected Python interpreter.

## Contents

| Directory | Calculation |
|---|---|
| original_reproduction | Original finite d=4 and d=6 examples, unmodified historical script, independent full-sector calculation; d=2,3 regressions |
| structured_benchmark | Physical matching tests, exact n=2 rational certificate, floating n=3 CRT optimization |
| lz_benchmark | Uniform Z/F MUB, CMUB, SMUB and SCMUB for n=2,3,4, every terminal, exact conditional support and numerical spectra |
| resources | Independent integer transfer contraction and explicit local control matrices |

`comparison.csv` and `comparison.json` apply one common independent-source, independently randomized all-pass model to the listed gaps. `MANUSCRIPT_MAPPING.md` identifies equations and SI sections. A single-copy test template is not a fixed physical measurement basis.

## Precision and scope

The independent original-example checks use complex128 and absolute tolerance 1e-12; the preserved historical script retains 1e-10. Li-Zhu branch support is computed exactly in the sixth-root cyclotomic integer ring before floating diagonalization. The n=2 CRT optimum 3/11 has a standard-library Fraction certificate; the n=3 primal/dual estimates are floating numerical candidates, not interval-certified bounds. Complete recorded outputs retain residuals and unrounded spectra.

The Li-Zhu reference benchmark is deliberately restricted to n=2,3,4. Its dense matrices are an auditing representation, not a necessary online cost. The separate transfer code operates on one observed record with linear fixed-dimensional arithmetic and uses arbitrary-precision integers.

The historical README is preserved as source material, including any archival suggestions made at that time. It does not report a public deposit. No repository URL, DOI or licence is invented here. This package is supplied as accompanying software.

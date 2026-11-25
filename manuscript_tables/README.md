# Exact formula tables

`generate_tables.py --output NEW_DIRECTORY` computes the frozen formulas with exact fractions. `width_rows.tex` supplies only six table rows for direct LaTeX inclusion, while `width_intervals.tex` includes a complete tabular environment. `claim_values.json` retains fractions and constructive gaps. The top-level revision runner re-generates and compares the JSON and row file byte-for-byte.

The two columns are necessary and sufficient widths, not exact optima. These are evaluations of proven formulas, not independent proofs, new optimization searches or dense-spectrum simulations for n=100/1000.

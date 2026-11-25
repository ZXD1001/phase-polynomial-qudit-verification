#!/usr/bin/env python3
"""Build common-model sample counts from this package's recorded results."""
import csv
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent
epsilon, delta = 0.01, 0.05
rows = []


def append(n, protocol, gap, evidence):
    copies = math.ceil(math.log(delta)/math.log1p(-gap*epsilon)) if gap > 1e-10 else None
    rows.append({'n': n, 'protocol': protocol, 'gap': gap, 'evidence': evidence,
                 'epsilon': epsilon, 'delta': delta, 'N_all_pass': copies})


for n in (2, 3, 4):
    append(n, 'adaptive_two_colors', .5, 'analytic')
    append(n, 'uniform_independent_HT_matching', 1/6, 'analytic')
    append(n, 'uniform_shared_index_matching', 1/6, 'analytic')
    append(n, 'optimal_in_matching_only_class', .25 if n == 2 else 1/6, 'analytic restricted optimum')
    protocols = json.loads((ROOT/'lz_benchmark'/f'results_n{n}.json').read_text())['protocols']
    for name in ('MUB', 'CMUB'):
        tested = [(key, value) for key, value in protocols.items() if key.startswith(name+'_terminal_')]
        key, value = max(tested, key=lambda item: item[1]['spectral_gap_clipped'])
        append(n, 'LZ_best_tested_'+key, value['spectral_gap_clipped'], 'finite-size floating spectrum; uniform basis weights')
    for name in ('SMUB_uniform_terminals', 'SCMUB_uniform_terminals'):
        append(n, 'LZ_'+name, protocols[name]['spectral_gap_clipped'], 'finite-size floating spectrum')
append(2, 'optimal_in_144_CRT_product_Pauli_class', 3/11, 'exact rational certificate')
structured = json.loads((ROOT/'structured_benchmark'/'benchmark_results.json').read_text())
append(3, 'constructed_CRT_product_Pauli_numeric_primal', structured['n3_numeric']['qutrit_primal_gap'], 'floating construction; not interval certified')
assert json.loads((ROOT/'original_reproduction'/'results'/'reproduction_results.json').read_text())['passed']
assert json.loads((ROOT/'lz_benchmark'/'verification_report.json').read_text())['status'] == 'PASS'
assert json.loads((ROOT/'resources'/'resource_checks.json').read_text())['all_checks_passed']
with (ROOT/'comparison.csv').open('w', newline='') as stream:
    writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
(ROOT/'comparison.json').write_text(json.dumps({
    'statistical_model': 'trusted ideal measurements; independent, possibly non-identical preparations; mean tested-state infidelity at least epsilon; all tests pass',
    'formula': 'ceil(log(delta)/log(1-gap*epsilon))',
    'not_applicable_without_extension': 'arbitrarily correlated sources, retained-output-state certification, imperfect-completeness measurements',
    'rows': rows}, indent=2)+'\n')
print(json.dumps({'all_recorded_checks_passed': True, 'comparison_rows': len(rows),
                  'illustrative_epsilon': epsilon, 'illustrative_delta': delta}, indent=2))

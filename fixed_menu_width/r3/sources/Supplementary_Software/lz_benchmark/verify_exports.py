#!/usr/bin/env python3
"""Independent floating-amplitude cross-check of exported cyclotomic branches."""
import json
from pathlib import Path
import numpy as np

folder = Path(__file__).resolve().parent
reports = {}
for n in (2, 3, 4):
    data = np.load(folder/f'conditional_branches_n{n}.npz')
    psi = data['target_state']
    d = 6
    root = np.exp(2j*np.pi/d)
    f = root**np.outer(np.arange(d), np.arange(d))/np.sqrt(d)
    maximum_amplitude_error = 0.
    maximum_probability_error = 0.
    maximum_zero_branch_float_norm = 0.
    for k, (terminal, front, mask) in enumerate(zip(data['terminals_zero_based'], data['front_sites_zero_based'], data['front_masks'])):
        b = np.ones((1, 1), dtype=np.complex128)
        for bit in mask:
            b = np.kron(b, f if bit else np.eye(d))
        psi_matrix = psi.reshape((d,)*n).transpose(list(front)+[terminal]).reshape(-1, d)
        amplitudes = b.conj().T@psi_matrix
        probabilities = data['branch_probability_numerators'][k]/data['branch_probability_denominators'][k]
        exact_amplitudes = data['conditional_vectors'][k]*np.sqrt(probabilities[:, None])
        maximum_amplitude_error = max(maximum_amplitude_error, float(np.max(abs(amplitudes-exact_amplitudes))))
        maximum_probability_error = max(maximum_probability_error, float(np.max(abs((abs(amplitudes)**2).sum(axis=1)-probabilities))))
        if (probabilities==0).any():
            maximum_zero_branch_float_norm = max(maximum_zero_branch_float_norm,float(np.linalg.norm(amplitudes[probabilities==0], axis=1).max()))
    r = json.loads((folder/f'results_n{n}.json').read_text())
    reports[n] = {
        'maximum_conditional_amplitude_absolute_error': maximum_amplitude_error,
        'maximum_branch_probability_absolute_error': maximum_probability_error,
        'maximum_exact_zero_branch_floating_norm': maximum_zero_branch_float_norm,
        'maximum_Q_target_residual_2': max(t['target_eigenvector_residual_2'] for t in r['tests']),
        'maximum_Q_isometry_gram_residual_frobenius': max(t['accepted_columns_orthonormality_residual_frobenius'] for t in r['tests']),
        'maximum_Omega_target_residual_2': max(p['target_eigenvector_residual_2'] for p in r['protocols'].values()),
        'maximum_unit_eigenvalue_error': max(abs(p['largest_eigenvalue']-1) for p in r['protocols'].values()),
        'minimum_Omega_eigenvalue': min(p['smallest_eigenvalue'] for p in r['protocols'].values()),
    }
    assert maximum_amplitude_error < 1e-12
    assert maximum_probability_error < 1e-12
    assert reports[n]['maximum_Q_target_residual_2'] < 1e-12
    assert reports[n]['maximum_Q_isometry_gram_residual_frobenius'] < 1e-11
    assert reports[n]['maximum_Omega_target_residual_2'] < 1e-12
    assert reports[n]['maximum_unit_eigenvalue_error'] < 1e-11
    for p in r['protocols'].values():
        if 'zero_gap_witness' in p:
            assert p['zero_gap_witness']['overlap_with_target_abs'] < 1e-12
            assert p['zero_gap_witness']['unit_eigenvector_residual_2'] < 1e-11
    if n <= 3:
        op = np.load(folder/f'operators_and_spectra_n{n}.npz')
        for terminal in range(n):
            all_q = [op[f'Q_t{terminal}_m'+''.join(map(str,m))] for m in data['front_masks'][data['terminals_zero_based']==terminal]]
            reconstructed = sum(all_q)/len(all_q)
            assert np.linalg.norm(reconstructed-op[f'Omega_MUB_terminal_{terminal+1}'])<1e-12
            zeros, ones = '0'*(n-1), '1'*(n-1)
            reconstructed = (op[f'Q_t{terminal}_m{zeros}']+op[f'Q_t{terminal}_m{ones}'])/2
            assert np.linalg.norm(reconstructed-op[f'Omega_CMUB_terminal_{terminal+1}'])<1e-12
reports['status'] = 'PASS'
(folder/'verification_report.json').write_text(json.dumps(reports, indent=2)+'\n')
print(json.dumps(reports, indent=2))

#!/usr/bin/env python3
"""Li--Zhu Sec.4 uniform Z/F6 protocols on an oriented d=6 path.

Small-system reference calculation, not an online-cost claim.  Exact integer
arithmetic in Z[exp(i*pi/3)] identifies impossible target branches before
double-precision construction/diagonalization of verification operators.
"""
from __future__ import annotations
import argparse
import hashlib
import itertools
import json
import platform
import time
from pathlib import Path
import numpy as np

D = 6
ROOT = np.exp(2j*np.pi/D)
PAIR_A = np.array([1, 0, -1, -1, 0, 1], dtype=np.int64)
PAIR_B = np.array([0, 1, 1, 0, -1, -1], dtype=np.int64)
FOURIER = ROOT ** np.outer(np.arange(D), np.arange(D)) / np.sqrt(D)


def path_phase(n):
    xx = np.indices((D,)*n, dtype=np.int64)
    return sum((2*xx[v]**2+3*xx[v])*xx[v+1] for v in range(n-1)) % D


def exact_fourier_bra(a, b, axis):
    exponent = (-np.outer(np.arange(D), np.arange(D))) % D
    c, e = PAIR_A[exponent], PAIR_B[exponent]
    ac = np.tensordot(c, a, axes=(1, axis))
    be = np.tensordot(e, b, axes=(1, axis))
    ae = np.tensordot(e, a, axes=(1, axis))
    bc = np.tensordot(c, b, axes=(1, axis))
    return np.moveaxis(ac-be, 0, axis), np.moveaxis(ae+bc+be, 0, axis)


def make_test(phase, terminal, mask):
    n = phase.ndim
    front = [v for v in range(n) if v != terminal]
    order = front + [terminal]
    a, b = PAIR_A[phase].copy(), PAIR_B[phase].copy()
    for axis, bit in zip(front, mask):
        if bit:
            a, b = exact_fourier_bra(a, b, axis)
    a = a.transpose(order).reshape(-1, D)
    b = b.transpose(order).reshape(-1, D)
    numerators = (a*a+a*b+b*b).sum(axis=1)
    denominator = D**(n+sum(mask))
    assert int(numerators.sum()) == denominator
    nonzero = numerators > 0
    q = np.zeros(a.shape, dtype=np.complex128)
    q[nonzero] = (a[nonzero]+ROOT*b[nonzero])/np.sqrt(numerators[nonzero, None])
    # Each column is |front basis outcome> tensor |conditional terminal state>.
    front_basis = np.ones((1, 1), dtype=np.complex128)
    for bit in mask:
        front_basis = np.kron(front_basis, FOURIER if bit else np.eye(D))
    v = np.einsum('io,oj->ijo', front_basis, q).reshape((D,)*n + (D**(n-1),))
    v = v.transpose([order.index(j) for j in range(n)] + [n]).reshape(D**n, -1)
    v = v[:, nonzero]
    label = f't{terminal}_m'+''.join(map(str, mask))
    return {
        'label': label, 'terminal_zero_based': terminal,
        'front_sites_zero_based': front, 'mask': list(mask),
        'conditional_vectors': q, 'integer_a': a, 'integer_b': b,
        'branch_probability_numerators': numerators,
        'branch_probability_denominator': denominator,
        'nonzero_branch_count': int(nonzero.sum()),
        'zero_branch_count': int((~nonzero).sum()), 'v': v,
    }


def metrics(operator, psi):
    op = (operator+operator.conj().T)/2
    # Full Hermitian spectrum; numerical example only, with no asymptotic claim.
    eig, vec = np.linalg.eigh(op)
    beta = float(eig[-2])
    result = {
        'spectral_gap_raw': 1-beta,
        'spectral_gap_clipped': max(0., 1-beta),
        'second_eigenvalue': beta,
        'largest_eigenvalue': float(eig[-1]),
        'smallest_eigenvalue': float(eig[0]),
        'eigenvalue_one_multiplicity_at_1e_minus_10': int(np.sum(abs(eig-1)<1e-10)),
        'target_acceptance': float(np.vdot(psi, op@psi).real),
        'target_eigenvector_residual_2': float(np.linalg.norm(op@psi-psi)),
        'hermiticity_residual_frobenius': float(np.linalg.norm(operator-operator.conj().T)),
    }
    witness = None
    # Provide a normalized, target-orthogonal accepted eigenvector if gap vanishes.
    if beta > 1-1e-10:
        candidates = vec[:, abs(eig-1)<1e-10]
        candidates -= np.outer(psi, psi.conj()@candidates)
        norms = np.linalg.norm(candidates, axis=0)
        witness = candidates[:, norms.argmax()]/norms.max()
        result['zero_gap_witness'] = {
            'overlap_with_target_abs': float(abs(np.vdot(psi, witness))),
            'acceptance': float(np.vdot(witness, op@witness).real),
            'unit_eigenvector_residual_2': float(np.linalg.norm(op@witness-witness)),
        }
    return result, eig, witness


def run_n(n, out):
    start = time.perf_counter()
    phase = path_phase(n)
    psi = (ROOT**phase).ravel()/np.sqrt(D**n)
    all_masks = list(itertools.product((0, 1), repeat=n-1))
    tests, test_diagnostics, saved = [], [], {}
    terminal_ops = {}
    for terminal in range(n):
        mub = np.zeros((D**n, D**n), dtype=np.complex128)
        cmub = np.zeros_like(mub)
        for mask in all_masks:
            test = make_test(phase, terminal, mask)
            v = test.pop('v')
            q_operator = v@v.conj().T
            mub += q_operator/len(all_masks)
            if len(set(mask)) == 1:
                cmub += q_operator/2
            gram_residual = np.linalg.norm(v.conj().T@v-np.eye(v.shape[1]))
            diagnostics = {
                key: test[key] for key in ('label', 'terminal_zero_based',
                    'front_sites_zero_based', 'mask', 'nonzero_branch_count',
                    'zero_branch_count', 'branch_probability_denominator')
            }
            diagnostics.update({
                'target_acceptance': float(np.vdot(psi, q_operator@psi).real),
                'target_eigenvector_residual_2': float(np.linalg.norm(q_operator@psi-psi)),
                'accepted_columns_orthonormality_residual_frobenius': float(gram_residual),
                'minimum_nonzero_branch_probability': float(test['branch_probability_numerators'][test['branch_probability_numerators']>0].min()/test['branch_probability_denominator']),
            })
            test_diagnostics.append(diagnostics)
            tests.append(test)
            if n <= 3:
                saved['Q_'+test['label']] = q_operator
        terminal_ops[f'MUB_terminal_{terminal+1}'] = mub
        terminal_ops[f'CMUB_terminal_{terminal+1}'] = cmub
    terminal_ops['SMUB_uniform_terminals'] = sum(terminal_ops[f'MUB_terminal_{k+1}'] for k in range(n))/n
    terminal_ops['SCMUB_uniform_terminals'] = sum(terminal_ops[f'CMUB_terminal_{k+1}'] for k in range(n))/n
    results = {}
    for name, op in terminal_ops.items():
        result, eig, witness = metrics(op, psi)
        is_sym = name.startswith('SMUB') or name.startswith('SCMUB')
        base_count = 2 if 'CMUB' in name else 2**(n-1)
        result['test_template_count'] = base_count*(n if is_sym else 1)
        result['weights'] = 'uniform over terminals and masks' if is_sym else 'uniform over masks'
        results[name] = result
        saved['eigenvalues_'+name] = eig
        if n <= 3:
            saved['Omega_'+name] = op
        if witness is not None:
            saved['zero_gap_witness_'+name] = witness
        print(f'n={n} {name}: gap={result["spectral_gap_raw"]:.12g} multiplicity={result["eigenvalue_one_multiplicity_at_1e_minus_10"]}', flush=True)
    branch_arrays = {
        'target_state': psi, 'phase_mod6': phase,
        'terminals_zero_based': np.array([t['terminal_zero_based'] for t in tests]),
        'front_sites_zero_based': np.array([t['front_sites_zero_based'] for t in tests]),
        'front_masks': np.array([t['mask'] for t in tests]),
        'conditional_vectors': np.array([t['conditional_vectors'] for t in tests]),
        'integer_a': np.array([t['integer_a'] for t in tests]),
        'integer_b': np.array([t['integer_b'] for t in tests]),
        'branch_probability_numerators': np.array([t['branch_probability_numerators'] for t in tests]),
        'branch_probability_denominators': np.array([t['branch_probability_denominator'] for t in tests]),
    }
    np.savez_compressed(out/f'conditional_branches_n{n}.npz', **branch_arrays)
    np.savez_compressed(out/f'operators_and_spectra_n{n}.npz', target_state=psi, **saved)
    result = {'n': n, 'd': D, 'hilbert_dimension': D**n,
              'target_norm_squared': float(np.vdot(psi, psi).real),
              'tests': test_diagnostics, 'protocols': results,
              'elapsed_seconds': time.perf_counter()-start}
    (out/f'results_n{n}.json').write_text(json.dumps(result, indent=2)+'\n')
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--n', nargs='+', type=int, default=[2, 3])
    parser.add_argument('--output', type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    assert all(2 <= n <= 4 for n in args.n), 'This reference calculator is limited to n=2,3,4.'
    args.output.mkdir(parents=True, exist_ok=True)
    meta = {
        'source': 'Li--Zhu Quantum 10, 2011 (2026), Secs.4.1--4.3, Eqs.(30)--(31); App.D.1',
        'family': 'psi(x)=6^(-n/2)*exp(2*pi*i*sum_{v=0}^{n-2}(2*x_v^2+3*x_v)*x_{v+1}/6)',
        'phase_orientation': 'left endpoint is the quadratic variable on every path edge',
        'basis': 'Z=computational; X basis column o has entries exp(+2*pi*i*x*o/6)/sqrt(6)',
        'zero_probability_branches': 'exact cyclotomic integer support; impossible target branches always reject',
        'probability_optimization': False,
        'basis_order': 'C order, last local index fastest; front sites increasing, terminal moved last for conditional arrays',
        'numerical_precision': 'numpy complex128/float64; exact int64 support/probability numerators',
        'n_values': args.n, 'python': platform.python_version(), 'numpy': np.__version__,
        'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'limitations': 'Small-n operator audit; no asymptotic conclusion and no assumption that dense state storage is needed online. Not an implementation of SD/CSD, 3MUB, or optimized variants.'}
    (args.output/'inputs_and_environment.json').write_text(json.dumps(meta, indent=2)+'\n')
    for n in args.n:
        run_n(n, args.output)


if __name__ == '__main__':
    main()

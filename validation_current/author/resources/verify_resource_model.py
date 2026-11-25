#!/usr/bin/env python3
"""Independent transfer contraction and explicit local-gate resource checks.

Requires Python >=3.10 and NumPy. No Zotero, network, SciPy, or private paths.
The transfer algorithm uses Python integers, not fixed-width integer arrays.
Run after ../lz_benchmark/benchmark_lz.py to cross-check its exported branches.
"""
from __future__ import annotations
import itertools
import json
import platform
import time
from pathlib import Path
import numpy as np

D = 6
ROOT = np.exp(2j * np.pi / D)
ROOT_PAIRS = ((1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1))
ONE, ZERO = (1, 0), (0, 0)
G = tuple((2*x*x+3*x) % D for x in range(D))


def add(u, v):
    return u[0]+v[0], u[1]+v[1]


def mul(u, v):
    a, b = u
    c, e = v
    return a*c-b*e, a*e+b*c+b*e


def norm2(u):
    a, b = u
    return a*a+a*b+b*b


def conditional_transfer(n, terminal, mask, outcomes):
    """Unnormalized cyclotomic coefficients and exact Born probability.

    Front sites increase in index, mask 0=Z, 1=positive Fourier ket basis.
    The returned numerator is divided by 6**(n+sum(mask)).
    At fixed d the recurrence uses O(n) ring operations; exact bit costs grow.
    """
    assert 0 <= terminal < n and len(mask) == len(outcomes) == n-1
    front = [j for j in range(n) if j != terminal]
    weights = {}
    for j, bit, outcome in zip(front, mask, outcomes):
        assert bit in (0, 1) and 0 <= outcome < D
        weights[j] = [ROOT_PAIRS[(-x*outcome) % D] if bit else
                      (ONE if x == outcome else ZERO) for x in range(D)]
    left, right = [ONE]*D, [ONE]*D
    multiplies = 0
    max_bits = 1
    for j in range(terminal):
        new = []
        for y in range(D):
            value = ZERO
            for x in range(D):
                value = add(value, mul(mul(left[x], weights[j][x]), ROOT_PAIRS[G[x]*y % D]))
                multiplies += 2
            new.append(value)
        left = new
        max_bits = max(max_bits, *(abs(z).bit_length() for p in left for z in p))
    for j in range(n-1, terminal, -1):
        new = []
        for x in range(D):
            value = ZERO
            for y in range(D):
                value = add(value, mul(ROOT_PAIRS[G[x]*y % D], mul(weights[j][y], right[y])))
                multiplies += 2
            new.append(value)
        right = new
        max_bits = max(max_bits, *(abs(z).bit_length() for p in right for z in p))
    vector = [mul(a, b) for a, b in zip(left, right)]
    multiplies += D
    numerator = sum(map(norm2, vector))
    denominator = D**(n+sum(mask))
    max_bits = max(max_bits, *(abs(z).bit_length() for p in vector for z in p))
    return vector, numerator, denominator, multiplies, max_bits


def normalized(vector, numerator):
    if numerator == 0:
        return np.zeros(D, dtype=complex)
    # Scale before converting Python integers to float. float(numerator) can
    # overflow even when the final normalized six-vector is well-conditioned.
    shift = max(abs(z).bit_length() for pair in vector for z in pair)
    scale = 1 << shift
    return np.array([a/scale+ROOT*(b/scale) for a, b in vector])/np.sqrt(numerator/(scale*scale))


def givens(a, b):
    r = np.sqrt(abs(a)**2+abs(b)**2)
    if r < 1e-15:
        return None
    return np.array([[np.conj(a)/r, np.conj(b)/r], [-b/r, a/r]], dtype=complex)


def embedded(g, low, high):
    result = np.eye(D, dtype=complex)
    result[np.ix_([low, high], [low, high])] = g
    return result


def compile_unitary(u):
    """Return chronological diagonal + two-level SU(2) circuit for a 6x6 U."""
    work = u.copy()
    reductions = []
    for col in range(D-1):
        for row in range(D-1, col, -1):
            if abs(work[row, col]) < 1e-14:
                continue
            g = givens(work[row-1, col], work[row, col])
            work = embedded(g, row-1, row) @ work
            reductions.append((row-1, row, g))
    diagonal = np.diag(work).copy()
    diagonal /= np.abs(diagonal)
    circuit = [(a, b, g.conj().T) for a, b, g in reversed(reductions)]
    reconstructed = np.diag(diagonal)
    for a, b, g in circuit:
        reconstructed = embedded(g, a, b) @ reconstructed
    return diagonal, circuit, float(np.linalg.norm(reconstructed-u, ord=2))


def compile_ray(q):
    """At most five adjacent-level SU(2) gates implement |q><q| acceptance."""
    current = q.copy()
    unitary = np.eye(D, dtype=complex)
    operations = []
    for row in range(D-1, 0, -1):
        if abs(current[row]) < 1e-14:
            continue
        g = givens(current[row-1], current[row])
        full = embedded(g, row-1, row)
        current = full @ current
        unitary = full @ unitary
        operations.append((row-1, row, g))
    accepted = np.outer(unitary[0].conj(), unitary[0])
    return len(operations), float(np.linalg.norm(accepted-np.outer(q, q.conj()), ord=2))


def encode_complex(a):
    a = np.asarray(a)
    return np.stack([a.real, a.imag], axis=-1).tolist()


def main():
    out = Path(__file__).resolve().parent
    report = {'python': platform.python_version(), 'numpy': np.__version__,
              'model': 'ideal continuous two-level SU(2) gates and diagonal phases',
              'cross_checks': [], 'large_n_transfer': []}
    max_ray_gates, max_ray_error = 0, 0.
    for n in (2, 3, 4):
        path = out.parent/'lz_benchmark'/f'conditional_branches_n{n}.npz'
        if not path.exists():
            raise FileNotFoundError(f'Run Li-Zhu benchmark first: {path}')
        data = np.load(path, allow_pickle=False)
        branch_count, max_error = 0, 0.
        for index, terminal in enumerate(data['terminals_zero_based']):
            mask = tuple(map(int, data['front_masks'][index]))
            for flat, outcomes in enumerate(itertools.product(range(D), repeat=n-1)):
                vector, num, den, _, _ = conditional_transfer(n, int(terminal), mask, outcomes)
                assert num == int(data['branch_probability_numerators'][index, flat])
                assert den == int(data['branch_probability_denominators'][index])
                assert [v[0] for v in vector] == data['integer_a'][index, flat].tolist()
                assert [v[1] for v in vector] == data['integer_b'][index, flat].tolist()
                q = normalized(vector, num)
                error = float(np.max(np.abs(q-data['conditional_vectors'][index, flat])))
                max_error = max(max_error, error)
                if num:
                    count, residual = compile_ray(q)
                    max_ray_gates = max(max_ray_gates, count)
                    max_ray_error = max(max_ray_error, residual)
                branch_count += 1
        assert max_error < 1e-12
        report['cross_checks'].append({'n': n, 'branches_checked': branch_count,
            'integer_coefficients_and_probabilities': 'exactly equal',
            'max_normalized_vector_error': max_error})
    for n in (20, 100, 500):
        start = time.perf_counter()
        vector, num, _, operations, bits = conditional_transfer(n, n//2, (1,)*(n-1), (0,)*(n-1))
        q = normalized(vector, num)
        gates, error = compile_ray(q)
        assert num > 0 and abs(np.vdot(q, q)-1) < 1e-12 and error < 1e-12
        report['large_n_transfer'].append({'n': n, 'ring_multiplications': operations,
            'max_coefficient_bits': bits, 'nonzero_branch': num > 0,
            'normalization_error': float(abs(np.vdot(q, q)-1)),
            'ray_SU2_gate_count': gates, 'ray_projector_error_2': error,
            'seconds_in_this_environment': time.perf_counter()-start})
    assert max_ray_gates <= 5 and max_ray_error < 1e-12
    report['conditional_rank_one_measurement'] = {
        'max_SU2_gates_observed': max_ray_gates, 'general_bound': D-1,
        'max_projector_error_2': max_ray_error}
    fourier = ROOT**np.outer(np.arange(D), np.arange(D))/np.sqrt(D)
    circuits = {'convention': 'Initial diagonal phases, then listed gates in chronological order; complex entry=[real,imag]', 'bases': []}
    crt_order = [(3*a+4*b) % D for a in range(2) for b in range(3)]
    f3 = np.exp(2j*np.pi*np.outer(np.arange(3), np.arange(3))/3)/np.sqrt(3)
    h2 = np.array([[1, 1], [1, -1]])/np.sqrt(2)
    menu_checks = []
    for t in range(3):
        phase = np.diag(ROOT**(2*t*np.arange(D)**2))
        basis = phase @ fourier.conj()
        unitary = basis.conj().T
        diagonal, gates, error = compile_unitary(unitary)
        crt_basis = np.kron(h2, np.diag(np.exp(2j*np.pi*t*np.arange(3)**2/3))@f3)
        crt_error = float(np.linalg.norm(basis[np.ix_(crt_order, crt_order)]-crt_basis, ord=2))
        rays = [compile_ray(basis[:, j]) for j in range(D)]
        assert error < 1e-12 and crt_error < 1e-12 and len(gates) <= 15
        assert max(r[1] for r in rays) < 1e-12 and max(r[0] for r in rays) <= 5
        circuits['bases'].append({'t': t, 'initial_diagonal': encode_complex(diagonal),
            'gates': [{'levels': [a, b], 'matrix': encode_complex(g)} for a, b, g in gates]})
        menu_checks.append({'t': t, 'SU2_gate_count': len(gates),
            'unitary_reconstruction_error_2': error, 'CRT_basis_error_2': crt_error,
            'binary_ray_max_SU2_gate_count': max(r[0] for r in rays),
            'binary_ray_max_projector_error_2': max(r[1] for r in rays)})
    report['three_fixed_adaptive_bases'] = menu_checks
    report['all_checks_passed'] = True
    (out/'compiled_three_bases.json').write_text(json.dumps(circuits, indent=2)+'\n')
    (out/'resource_checks.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()

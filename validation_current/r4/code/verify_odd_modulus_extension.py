#!/usr/bin/env python3
"""Exact finite checks for the 2d-dimensional extension, d odd (not necessarily prime).

Only unit differences between chirp coefficients are used. Squared Gauss sums
are reduced modulo the actual cyclotomic polynomial Phi_d, rather than the
prime-only relation 1+w+...+w**(d-1)=0. This is a regression of the analytic
proof, not a numerical optimization of the whole verification class.
"""
from __future__ import annotations

import argparse
from functools import lru_cache
import json
from math import gcd
from pathlib import Path
import time

if not __debug__:
    raise RuntimeError('Run exact checks without -O or -OO.')


def trim(a: list[int]) -> list[int]:
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    return a


def div_monic(a: list[int], b: list[int]) -> tuple[list[int], list[int]]:
    if b[-1] != 1:
        raise ValueError('Only monic divisors are supported.')
    a = trim(a.copy())
    q = [0] * max(1, len(a) - len(b) + 1)
    while len(a) >= len(b) and a != [0]:
        shift, coeff = len(a) - len(b), a[-1]
        q[shift] += coeff
        for j, value in enumerate(b):
            a[shift + j] -= coeff * value
        trim(a)
    return trim(q), trim(a)


@lru_cache(None)
def cyclotomic(d: int) -> tuple[int, ...]:
    if d < 1:
        raise ValueError('Positive modulus required.')
    poly = [-1] + [0] * (d - 1) + [1]
    for divisor in range(1, d):
        if d % divisor == 0:
            poly, rem = div_monic(poly, list(cyclotomic(divisor)))
            assert rem == [0]
    return tuple(poly)


def norm_poly(d: int, a: int, b: int) -> list[int]:
    phase = [(a*x*x+b*x) % d for x in range(d)]
    coeff = [0]*d
    for x in range(d):
        for y in range(d):
            coeff[(phase[x]-phase[y]) % d] += 1
    return div_monic(coeff, list(cyclotomic(d)))[1]


def check(d: int) -> dict:
    if d < 3 or d % 2 != 1:
        raise ValueError('d must be odd and at least three.')
    units = [a for a in range(d) if gcd(a, d) == 1]
    checks = 0
    for a in units:
        for b in range(d):
            assert norm_poly(d,a,b) == [d]
            checks += 1
    for x in range(2*d):
        for y in range(2*d):
            a,b=x%2,x%d
            c,e=y%2,y%d
            assert x == (d*a+(d+1)*b) % (2*d)
            assert ((2*x*x+d*x)*y-d*a*c-2*b*b*e) % (2*d) == 0
    return {'odd_modulus':d, 'native_dimension':2*d,
            'complete_menu_size':d+1,'phi_d':len(units),
            'unit_mismatch_probability':f'{len(units)}/{d}',
            'unit_gauss_norm_cases':checks,'crt_edge_cases':4*d*d,
            'cyclotomic_polynomial_low_first':list(cyclotomic(d)),
            'status':'EXACT_PASS'}


def main() -> None:
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    start=time.monotonic()
    cases=[check(d) for d in (3,5,7,9,11,15,21,25,27,35)]
    # Nonunit differences are NOT universally unbiased. This exact example
    # prevents silently extending the prime-only probability formula.
    assert norm_poly(9,3,0) == [27]
    out={'status':'EXACT_PASS','cases':cases,
         'nonunit_counterexample':{'d':9,'quadratic_coefficient':3,
             'linear_coefficient':0,'unnormalized_squared_sum':27,
             'normalized_overlap_squared':'1/3',
             'mutually_unbiased_value_would_be':'1/9'},
         'scope':'Finite cyclotomic checks only. The analytic width threshold '
                 'extends to every odd d>=3 via unit differences; the sharper '
                 'prime formula must not be assigned to composite d.',
         'seconds':time.monotonic()-start}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(out,indent=2,ensure_ascii=False)+'\n')
    print(json.dumps(out,indent=2,ensure_ascii=False))

if __name__=='__main__':
    main()

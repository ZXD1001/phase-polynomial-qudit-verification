#!/usr/bin/env python3
"""Explicit positive character-error witness for near-optimal adaptive width.

All fractions are exact. The infinite-family proof is in REPORT.md; finite
verification must not be substituted for that proof. Site zero is the left bit.
"""
from __future__ import annotations
from collections import defaultdict
from fractions import Fraction as F
from pathlib import Path
import argparse, json
if not __debug__: raise RuntimeError('Do not run verification with -O/-OO.')

def generate(n:int,r:int)->dict:
    if not isinstance(n,int) or not isinstance(r,int) or n<3 or not 0<=r<n:
        raise ValueError('Require integers n>=3 and 0<=r<n.')
    delta=int(n%2==0);D=17*n-31+delta
    terms=defaultdict(F);zero=(0,)*n
    def bit(*sites):
        v=[0]*n
        for j in sites:v[j]^=1
        return tuple(v)
    for j in range(n-1):terms[zero,bit(j)]+=F(3,D)
    for j in range(n-2):
        terms[bit(j,j+2),zero]+=F(7,D)
        terms[bit(j+1),zero]+=F(7,D)
    if delta:terms[bit(n-1),zero]+=F(1,D)
    assert sum(terms.values())==1
    assert all(w>0 and (any(a) or any(b)) for (a,b),w in terms.items())
    raw=F(1,2)-F(n-1-2*r,D)
    m=[len(range(c,n-1,2)) for c in (0,1)]
    q=[F(1,3)+F(2*min(r,z),3*z) if z else F(1) for z in m]
    v=q[0]*q[1]/sum(q)
    return dict(n=n,adaptive_width=r,menu='Z,B0,B1,B2',
        claim_upper=str(raw),spectral_upper_with_menu_cap=str(min(F(1,2),raw)),
        denominator=D,achievable_lower=str(v),
        terms=[dict(qubit_character=list(a),qutrit_character=list(b),weight=str(w))
               for (a,b),w in sorted(terms.items())],
        class_description='At most r preselected second-layer sites; complete native records; '
        'arbitrary randomized/correlated second-layer menu choices; strongest support acceptance.',
        evidence='Analytic family certificate. claim_upper is the witness bound BEFORE the '
        'separate menu-wide 1/2 spectral cap; it can exceed 1/2 for large widths.')

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('n',type=int);ap.add_argument('r',type=int)
    ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=generate(a.n,a.r)
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2)+'\n')
    print(json.dumps({k:d[k] for k in ('n','adaptive_width','claim_upper','spectral_upper_with_menu_cap','denominator','achievable_lower')}))
if __name__=='__main__':main()

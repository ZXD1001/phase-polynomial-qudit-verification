#!/usr/bin/env python3
"""Generate an O(n)-support exact witness for the fixed-menu width bound.

The validity for every n is proved in REPORT.md. This generator uses only exact
Fractions. Output follows the round-2 finite verifier's schema as well.
Sparse bitstrings printed in the JSON use site 0 at the left.
"""
from __future__ import annotations

if not __debug__:
    raise RuntimeError("Run exact checks without -O or -OO.")
from collections import defaultdict
from fractions import Fraction as F
from pathlib import Path
import argparse, json


def parameters(n: int, r: int) -> dict:
    if not isinstance(n, int) or not isinstance(r, int) or n < 2 or not 0 <= r < n:
        raise ValueError('Require integer n >= 2 and 0 <= r < n.')
    groups = [[v for v in range(n-1) if v % 2 == c] for c in (0, 1)]
    m = [len(x) for x in groups]
    q = [F(1) if z == 0 else F(1,3)+F(2*min(r,z),3*z) for z in m]
    alpha = [q[1]/sum(q), q[0]/sum(q)]
    lower = q[0]*q[1]/sum(q)
    kappa = F(1, 4*n-6)
    upper = F(1,2)-(F(1,2)-lower)/F(4*n-5)
    return dict(n=n, r=r, groups=groups, m=m, q=q, alpha=alpha,
                lower=lower, upper=upper, kappa=kappa, bad_coefficient=4*n-5)


def generate(n: int, r: int) -> dict:
    p=parameters(n,r)
    z=(0,)*n
    terms=defaultdict(F)
    def bit(*positions):
        a=[0]*n
        for v in positions:a[v]^=1
        return tuple(a)
    # Graph part: 2(n-1) edge terms and 2(n-2) triple terms.
    d=4*n-5
    for v in range(n-1):
        terms[bit(v),z]+=F(1,d)
        terms[bit(v+1),z]+=F(1,d)
    for v in range(n-2):
        terms[bit(v,v+2),z]+=F(1,d)
        terms[bit(v+1),z]+=F(1,d)
    # Non-last qutrit errors, with the single last vertex for an empty group.
    for c in (0,1):
        sites=p['groups'][c] or [n-1]
        for v in sites:
            terms[z,bit(v)]+=p['alpha'][c]/len(sites)/d
    assert sum(terms.values())==1
    assert all(w>0 and (any(a) or any(b)) for (a,b),w in terms.items())
    return {
        'status':'ANALYTIC_FAMILY_CERTIFICATE_INPUT',
        'n':n,'adaptive_width':r,'menu':'Z,B0,B1,B2',
        'class':'at most r preselected second-layer sites; arbitrary front product menu; '
                'arbitrary full-transcript-dependent parallel second-layer product menu; '
                'strongest perfectly complete all-outcome support acceptance',
        'claim_upper':str(p['upper']),'achievable_lower':str(p['lower']),
        'minimum_width_for_half':n//2,'m':p['m'],
        'q':[str(x) for x in p['q']],
        'qutrit_color_weights':[str(x) for x in p['alpha']],
        'graph_noncolor_penalty':str(p['kappa']),
        'terms':[{'qubit_character':list(a),'qutrit_character':list(b),'weight':str(w)}
                 for (a,b),w in sorted(terms.items())],
        'evidence':'The infinite family follows from the analytic lemmas in REPORT.md; '
                   'finite tests are independent implementation checks, not the proof.'
    }


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('n',type=int);ap.add_argument('r',type=int)
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args();cert=generate(args.n,args.r)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(cert,indent=2)+'\n')
    print(json.dumps({k:cert[k] for k in ('n','adaptive_width','claim_upper','achievable_lower','minimum_width_for_half')}))

if __name__=='__main__':main()

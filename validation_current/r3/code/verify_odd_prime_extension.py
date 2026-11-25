#!/usr/bin/env python3
"""Exact finite regression of the odd-prime 2p-dimensional corollary.

Cyclotomic squared sums are reduced modulo 1+w+...+w^(p-1)=0.
No floating trigonometric tolerance is used.
"""
from pathlib import Path

if not __debug__:
    raise RuntimeError("Run exact checks without -O or -OO.")
from fractions import Fraction as F
import argparse,json,time


def prime(p):return p>=2 and all(p%k for k in range(2,int(p**0.5)+1))

def check(p):
    if p%2==0 or not prime(p):raise ValueError('p must be an odd prime')
    ctr=0
    for c in range(p):
        for linear in range(p):
            ph=[(c*x*x+linear*x)%p for x in range(p)]
            coeff=[0]*p
            for x in range(p):
                for y in range(p):coeff[(ph[x]-ph[y])%p]+=1
            reduced=[coeff[k]-coeff[-1] for k in range(p-1)]
            expected=p if c else p*p*int(linear==0)
            assert reduced==[expected]+[0]*(p-2)
            ctr+=1
    for x in range(2*p):
        for y in range(2*p):
            a,b=x%2,x%p;c,e=y%2,y%p
            assert x==(p*a+(p+1)*b)%(2*p)
            assert ((2*x*x+p*x)*y-p*a*c-2*b*b*e)%(2*p)==0
    return {'p':p,'native_dimension':2*p,'complete_local_menu_size':p+1,
            'gauss_norm_cases':ctr,'crt_edge_cases':(2*p)**2,'status':'EXACT_PASS'}


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();t=time.monotonic()
    out={'status':'EXACT_PASS','cases':[check(x) for x in (3,5,7,11,13,17,19)],
         'scope':'Finite regression only. Analytic corollary applies to every odd prime, '
                 'not odd composite p or p=2.','seconds':time.monotonic()-t}
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
if __name__=='__main__':main()

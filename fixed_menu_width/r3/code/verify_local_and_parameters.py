#!/usr/bin/env python3
"""Exact Gauss-sum checks and symbolic-family parameter checks.

The finite parameter sweep is a regression, not the proof for all n.
The generator also emits explicit n=7,8 inputs; these are theorem instances,
not claims that all original decision trees at these sizes were enumerated.
"""
from __future__ import annotations

if not __debug__:
    raise RuntimeError("Run exact checks without -O or -OO.")
from pathlib import Path
from fractions import Fraction as F
import argparse,json,time
from family_certificate import parameters,generate
ROOT=((1,0),(0,1),(-1,-1))

def norm2(a,b):return a*a-a*b+b*b

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--inputs',type=Path,required=True)
    args=ap.parse_args();start=time.monotonic();gauss=[]
    for r in range(3):
        for t in range(3):
            for linear in range(3):
                values=[]
                for shift in range(3):
                    raw=[tuple(ROOT[((r-t)*x*x+(linear+shift)*x)%3]) for x in range(3)]
                    a=sum(v[0] for v in raw);b=sum(v[1] for v in raw)
                    q=F(norm2(a,b),9)
                    expected=F(1,3) if r!=t else F(int((linear+shift)%3==0))
                    assert q==expected;values.append(str(q))
                if r!=t:assert len(set(values))==1
                gauss.append({'target_chirp':r,'measurement_chirp':t,'linear_term':linear,
                              'probabilities_by_character_shift':values})
    count=0;monotone=True
    for n in range(2,501):
        lastlo=F(0);lasthi=F(0)
        for r in range(n):
            p=parameters(n,r)
            assert p['lower']<=p['upper']<=F(1,2)
            assert (p['upper']==F(1,2))==(r>=n//2)
            assert (p['lower']==F(1,2))==(r>=n//2)
            assert p['lower']>=lastlo and p['upper']>=lasthi
            assert sum(p['alpha'])==1
            assert p['alpha'][0]*p['q'][0]==p['alpha'][1]*p['q'][1]==p['lower']
            lastlo,lasthi=p['lower'],p['upper'];count+=1
    args.inputs.mkdir(parents=True,exist_ok=True)
    examples=[]
    for n,r in [(2,0),(2,1),(3,0),(3,1),(4,1),(4,2),(5,1),(5,2),(6,2),(6,3),(7,2),(7,3),(8,3),(8,4),(100,49),(100,50)]:
        cert=generate(n,r)
        path=args.inputs/f'family_n{n}_w{r}.json';path.write_text(json.dumps(cert,indent=2)+'\n')
        examples.append({k:cert[k] for k in ('n','adaptive_width','claim_upper','achievable_lower','minimum_width_for_half')})
    out={'status':'EXACT_PASS','gauss_cases':gauss,'gauss_norm_identities':81,
         'finite_parameter_checks':count,'parameter_range':'2<=n<=500; 0<=r<n',
         'examples':examples,'seconds':time.monotonic()-start}
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k!='gauss_cases'},indent=2))
if __name__=='__main__':main()

#!/usr/bin/env python3
"""Evaluate frozen R2/R4 formulas exactly and generate manuscript table inputs.

This is arithmetic transcription QA, not a proof checker or a dense-spectrum
calculation. No optimization over measurements is performed.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as F
import json
from pathlib import Path
import sys


def ceil(x:F)->int:
    return -((-x.numerator)//x.denominator)


def achievable(n:int,r:int)->F:
    if n<2 or not 0<=r<=n:raise ValueError('Require n>=2 and 0<=r<=n.')
    m=[len(range(c,n-1,2)) for c in (0,1)]
    q=[F(1) if v==0 else F(1,3)+F(2*min(r,v),3*v) for v in m]
    return q[0]*q[1]/sum(q)


def necessary(n:int,loss:F)->int:
    if n<3 or not F(0)<=loss<=F(1,2):raise ValueError('Require n>=3, 0<=loss<=1/2.')
    d=17*n-31+int(n%2==0)
    return max(0,ceil((F(n-1)-d*loss)/2))


def sufficient(n:int,loss:F)->int:
    goal=F(1,2)-loss
    return next(r for r in range(n//2+1) if achievable(n,r)>=goal)


def data()->dict:
    expected={(100,F(1,100)):(42,49),(100,F(1,50)):(33,47),(100,F(1,25)):(17,44),
              (1000,F(1,100)):(415,485),(1000,F(1,50)):(330,470),(1000,F(1,25)):(161,440)}
    rows=[]
    for (n,loss),(lo,hi) in expected.items():
        got=(necessary(n,loss),sufficient(n,loss))
        if got!=(lo,hi):raise AssertionError((n,str(loss),got,(lo,hi)))
        rows.append(dict(n=n,ideal_gap=str(F(1,2)-loss),loss=str(loss),necessary_width=lo,sufficient_width=hi,
                         sufficient_strategy_gap=str(achievable(n,hi))))
    finite=[]
    for n,r,low,up in [(4,1,'2/5','29/70'),(5,1,'1/3','26915/68048'),(6,2,'7/16','10049/22500')]:
        if achievable(n,r)!=F(low):raise AssertionError((n,r))
        d=17*n-31+int(n%2==0)
        finite.append(dict(n=n,width=r,constructive_lower=low,round2_exact_upper=up,
                           round4_generic_upper=str(F(1,2)-F(max(0,n-1-2*r),d)),
                           interpretation='Lower and upper bounds; not an exact class optimum.'))
    for n in range(3,251):
        if necessary(n,F(0))!=n//2 or achievable(n,n//2)!=F(1,2):raise AssertionError(n)
    if achievable(2,0)!=F(1,4) or achievable(2,1)!=F(1,2):raise AssertionError('n=2')
    return {'status':'EXACT_ARITHMETIC_PASS','scope':'Evaluation of supplied formulas and table transcription; not a new theorem proof.',
            'source':'Supplementary Information: fixed-menu width proofs; finite certificates in fixed_menu_width/r2' ,
            'near_optimal_width_rows':rows,'finite_width_rows':finite}


def main()->int:
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path);a=p.parse_args()
    result=data()
    if a.output:
        dest=a.output.expanduser().resolve()
        if dest.exists():raise FileExistsError('Output directory must be new.')
        dest.mkdir(parents=True)
        (dest/'claim_values.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8',newline='\n')
        text=['% Evaluated from exact formulas, not dense spectra or exact optimal widths.',
              '\\begin{tabular}{rrrr}','\\hline',
              '$n$ & Required ideal gap & Necessary width & Sufficient width \\\\','\\hline']
        for row in result['near_optimal_width_rows']:
            text.append(f"{row['n']} & {float(F(row['ideal_gap'])):.2f} & {row['necessary_width']} & {row['sufficient_width']} \\\\")
        text+=['\\hline','\\end{tabular}']
        (dest/'width_intervals.tex').write_text('\n'.join(text)+'\n',encoding='utf-8',newline='\n')
        rows_only=[f"{row['n']} & {float(F(row['ideal_gap'])):.2f} & {row['necessary_width']} & {row['sufficient_width']} \\\\" for row in result['near_optimal_width_rows']]
        (dest/'width_rows.tex').write_text('\n'.join(rows_only)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(result,ensure_ascii=False,indent=2));return 0

if __name__=='__main__':
    if not __debug__:raise RuntimeError('Do not use Python -O/-OO.')
    try:raise SystemExit(main())
    except (OSError,ValueError,AssertionError) as exc:
        print(f'ERROR: {exc}',file=sys.stderr);raise SystemExit(1)

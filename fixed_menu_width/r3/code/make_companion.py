#!/usr/bin/env python3
"""Coverage input for exact-size one using the n=6, budget-two witness."""
from pathlib import Path
import argparse,json
if not __debug__:raise RuntimeError('Run without -O/-OO.')

def main():
    p=argparse.ArgumentParser();p.add_argument('--inputs',type=Path,required=True);a=p.parse_args()
    x=json.loads((a.inputs/'family_n6_w2.json').read_text())
    x['adaptive_width']=1;x['witness_parameter_budget']=2
    x['class']='exactly one second-layer site; companion coverage for the at-most-two result using sigma_(6,2)'
    x['companion_note']='The witness and claimed upper are generated for r=2; only the enumerated second-layer size is one.'
    for k in ('achievable_lower','q','qutrit_color_weights'):x.pop(k,None)
    (a.inputs/'family_n6_w1_companion.json').write_text(json.dumps(x,indent=2)+'\n')
    print('PASS: generated exact-size-one companion; not a width-one achievable-gap claim.')
if __name__=='__main__':main()

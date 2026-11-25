#!/usr/bin/env python3
"""Re-evaluate the supplied n=4 Z/F floating feasible weights; no optimization.

This numerical calculation is not an interval-certified lower bound or an exact
class optimum. The separately supplied exact upper certificate is unaffected.
"""
from pathlib import Path
import json
import numpy as np
from scipy.linalg import eigvalsh
from lz_optimize import build
if not __debug__:
    raise RuntimeError('Do not use -O/-OO.')
ROOT=Path(__file__).resolve().parent
c=json.loads((ROOT/'lz_optimized_n4_all.json').read_text())
labels,Q,a,psi,V=build(4)
assert json.loads(json.dumps(labels))==c['labels']
p=np.asarray(c['weights'],dtype=float)
assert p.shape==(32,) and np.all(p>=0)
assert abs(float(p.sum())-1)<1e-12
omega=np.einsum('i,ijk->jk',p,Q)
beta3=float(eigvalsh(V.conj().T@omega@V)[-1])
beta2=float(np.max(p@a[:,1:]))
gap=1-max(beta3,beta2)
residual=abs(gap-c['primal_gap'])
assert residual<1e-12
assert gap<float(197537/500000)
out={'status':'NUMERIC_PASS','evidence_type':'NUMERIC','n':4,'class':'32 terminal/Z-F front settings; unrestricted conditional terminal rays','primal_gap':gap,'reported_candidate_gap':c['primal_gap'],'candidate_gap_residual':residual,'normalization_residual':abs(float(p.sum())-1),'tolerance':1e-12,'scope':'Feasible-weight floating spectrum; not interval certification or an exact class optimum.'}
(ROOT/'n4_primal_recomputed.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))

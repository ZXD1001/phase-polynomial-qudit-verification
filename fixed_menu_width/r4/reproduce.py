#!/usr/bin/env python3
"""Reproduce the round-4 checks in a NEW directory, preserving packaged sources.

The author's full pipeline is NOT re-executed: its received manifest and real
round-3 receipt are checked. This distinction is recorded in the results.
No web access or numerical optimizer is used by this entry point.
"""
from __future__ import annotations
import argparse,hashlib,json,os,platform,shutil,subprocess,sys,time
from pathlib import Path
if not __debug__:raise RuntimeError('Do not use -O/-OO.')
ROOT=Path(__file__).resolve().parent
FOLDERS=('code','inputs','source_archives','reference_receipts','prior_checks')
def hashes(root):
    out={}
    for folder in FOLDERS:
      for p in sorted((root/folder).rglob('*')):
        if not p.is_file() or '__pycache__' in p.parts:continue
        if folder=='prior_checks' and p.name.endswith('_result.json'):continue
        out[str(p.relative_to(root))]=hashlib.sha256(p.read_bytes()).hexdigest()
    return out

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--destination',type=Path,required=True);args=ap.parse_args()
    dest=args.destination.resolve()
    if dest.exists():raise SystemExit('Destination must not already exist.')
    before=hashes(ROOT);dest.mkdir(parents=True)
    for folder in FOLDERS:
        shutil.copytree(ROOT/folder,dest/folder,ignore=shutil.ignore_patterns('__pycache__','*_result.json'))
    (dest/'results').mkdir();(dest/'logs').mkdir()
    steps=[
      ['code/audit_evidence.py','--output','results/evidence_audit.json'],
      ['code/verify_graph_lemmas.py','--max-n','12','--walsh-max-n','6','--output','results/round3_graph_audit.json'],
      ['code/verify_robust_combinatorics.py','--max-n','12','--mask-max-n','8','--output','results/robust_combinatorics.json'],
      ['code/verify_odd_modulus_extension.py','--output','results/round3_odd_modulus_audit.json'],
      ['prior_checks/exact_certificate.py'],
      ['prior_checks/certify_fullmenu_half.py','--verify-only'],
    ]
    for n,rs in [(3,[0,1,2]),(4,[0,1,2,3]),(5,[0,1,2]),(6,[1,2,3])]:
      for r in rs:steps.append(['code/verify_width_certificate.py',f'inputs/robust_n{n}_w{r}.json','--native'])
    for n in (3,4):steps.append(['code/verify_direct_and_local.py','--n',str(n),'--output',f'results/direct_local_n{n}.json'])
    import numpy,scipy,sympy
    report=dict(status='RUNNING',python=sys.version,numpy=numpy.__version__,scipy=scipy.__version__,sympy=sympy.__version__,platform=platform.platform(),steps=[],source_sha256=before)
    env=os.environ.copy();env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1')
    start=time.monotonic()
    for k,cmd in enumerate(steps,1):
        print(f'[{k}/{len(steps)}] {" ".join(cmd)}',flush=True);t=time.monotonic()
        x=subprocess.run([sys.executable,*cmd],cwd=dest,env=env,capture_output=True,text=True)
        (dest/'logs'/f'{k:02d}_stdout.txt').write_text(x.stdout);(dest/'logs'/f'{k:02d}_stderr.txt').write_text(x.stderr)
        report['steps'].append(dict(command=cmd,exit_code=x.returncode,seconds=time.monotonic()-t));report['seconds']=time.monotonic()-start
        if x.returncode:
            report['status']='FAILED';(dest/'run_receipt.json').write_text(json.dumps(report,indent=2)+'\n');print(x.stderr,file=sys.stderr);return x.returncode
        (dest/'run_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    assert before==hashes(ROOT),'Packaged sources changed.'
    report.update(status='PASS',source_bytes_unchanged=True,seconds=time.monotonic()-start,
                  evidence_boundary='Finite checks support implementations; the all-n proof is in REPORT.md. Author full six-stage rerun is historical round-3 evidence, not re-executed in this round.')
    (dest/'run_receipt.json').write_text(json.dumps(report,indent=2)+'\n');print(f'PASS {len(steps)} steps in {report["seconds"]:.2f}s',flush=True);return 0
if __name__=='__main__':raise SystemExit(main())

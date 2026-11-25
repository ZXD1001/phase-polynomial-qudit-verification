#!/usr/bin/env python3
"""Run all round-3 finite checks and the author's six stages in a NEW directory.

No optimizer is required. Exact tests support (but do not replace) the analytic
proof. Existing source/results files are preserved. Never use -O or -OO.
"""
from __future__ import annotations
import argparse,hashlib,json,os,platform,shutil,subprocess,sys,time
from pathlib import Path
if not __debug__:raise RuntimeError('Exact checks must not run under -O/-OO.')
ROOT=Path(__file__).resolve().parent

def hashes(root):
    out={}
    for directory in ('code','inputs','sources/Supplementary_Software'):
        for p in sorted((root/directory).rglob('*')):
            if p.is_file() and '__pycache__' not in p.parts:
                out[str(p.relative_to(root))]=hashlib.sha256(p.read_bytes()).hexdigest()
    return out

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--destination',type=Path,required=True)
    args=ap.parse_args();dest=args.destination.resolve()
    if dest.exists():raise SystemExit('Destination already exists. Supply a NEW directory.')
    before=hashes(ROOT);dest.mkdir(parents=True)
    for folder in ('code','inputs','sources/Supplementary_Software'):
        shutil.copytree(ROOT/folder,dest/folder,ignore=shutil.ignore_patterns('__pycache__'))
    (dest/'results').mkdir();(dest/'logs').mkdir()
    steps=[
      ['code/verify_graph_lemmas.py','--max-n','12','--walsh-max-n','6','--output','results/graph_lemmas.json'],
      ['code/verify_local_and_parameters.py','--output','results/local_and_parameters.json','--inputs','inputs'],
      ['code/make_companion.py','--inputs','inputs'],
      ['code/verify_odd_prime_extension.py','--output','results/odd_prime_extension.json'],
      ['code/verify_odd_modulus_extension.py','--output','results/odd_modulus_extension.json'],
      ['code/verify_width_certificate.py','inputs/family_n4_w1.json','--native'],
      ['code/verify_width_certificate.py','inputs/family_n5_w1.json','--native'],
      ['code/verify_width_certificate.py','inputs/family_n6_w1_companion.json','--native'],
      ['code/verify_width_certificate.py','inputs/family_n6_w2.json','--native'],
      ['code/verify_family_direct_omega6.py'],
      ['code/audit_author.py','--source','sources/Supplementary_Software','--destination','author_rerun','--output','results/author_audit.json'],
    ]
    env=os.environ.copy();env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1')
    import numpy,scipy
    report={'status':'RUNNING','python':sys.version,'numpy':numpy.__version__,'scipy':scipy.__version__,
      'platform':platform.platform(),'steps':[],'source_sha256':before,'note':'Finite checks are not a formal proof for all n; see REPORT.md.'}
    start=time.monotonic()
    for k,cmd in enumerate(steps,1):
        print(f'[{k}/{len(steps)}] {" ".join(cmd)}',flush=True)
        t=time.monotonic();result=subprocess.run([sys.executable,*cmd],cwd=dest,env=env,capture_output=True,text=True)
        (dest/'logs'/f'{k:02d}_stdout.txt').write_text(result.stdout)
        (dest/'logs'/f'{k:02d}_stderr.txt').write_text(result.stderr)
        report['steps'].append({'command':cmd,'exit_code':result.returncode,'seconds':time.monotonic()-t})
        report['seconds']=time.monotonic()-start
        if result.returncode:
            report['status']='FAILED';(dest/'run_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
            print(result.stderr,file=sys.stderr);return result.returncode
    assert before==hashes(ROOT)
    report.update(status='PASS',source_bytes_unchanged=True,seconds=time.monotonic()-start)
    (dest/'run_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    print(f'PASS: {len(steps)} stages (including the author six-stage pipeline); {report["seconds"]:.2f}s',flush=True)
    return 0
if __name__=='__main__':raise SystemExit(main())

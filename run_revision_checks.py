#!/usr/bin/env python3
"""Run the preserved baseline and four revision groups in a NEW directory.

Requires Python without -O/-OO and the versions in requirements.txt. Sources,
certificates and supplied results are never overwritten. All scientific modules
are packaged here; the selected Python environment supplies numerical libraries.
"""
from __future__ import annotations
import argparse
import datetime
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import shlex
import shutil
import subprocess
import sys
import sysconfig

# Keep the supplied package unchanged even when invoked without environment flags.
sys.dont_write_bytecode = True
from portability.import_audit import closure, install_launcher
import time

ROOT=Path(__file__).resolve().parent
if not __debug__:raise RuntimeError('Do not use -O/-OO; exact-check assertions must be enabled.')

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def hashes(root):
    return {str(p.relative_to(root)):sha(p) for p in sorted(root.rglob('*'))
            if p.is_file() and '__pycache__' not in p.parts and not any(x.startswith('.') for x in p.relative_to(root).parts)
            and 'validation_current' not in p.relative_to(root).parts}
def write(path,obj):path.write_text(json.dumps(obj,indent=2,ensure_ascii=False)+'\n',encoding='utf-8',newline='\n')
def copytree(a,b):shutil.copytree(a,b,ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--destination',type=Path,required=True)
    a=p.parse_args();dest=a.destination.expanduser().resolve()
    if dest.exists():p.error('The output directory already exists. Choose a new --destination; existing results are never overwritten.')
    if dest==ROOT or ROOT in dest.parents:p.error('Choose an output directory outside the supplied package, for example ../revision-rerun.')
    expected={'numpy':'2.3.5','scipy':'1.17.0','sympy':'1.14.0'}
    dependencies={};missing=[]
    for name in expected:
        try:dependencies[name]=importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:missing.append(name)
    if missing:
        command=[sys.executable,'-m','pip','install','-r',str(ROOT/'requirements.txt')]
        install_command=('& '+' '.join("'"+arg.replace("'","''")+"'" for arg in command)) if os.name=='nt' else shlex.join(command)
        p.error('Missing dependencies: '+', '.join(missing)+'. Install them with the same Python interpreter:\n'+install_command)
    source_hashes=hashes(ROOT)
    dest.mkdir(parents=True);(dest/'logs').mkdir();(dest/'import_audit').mkdir()
    report={'status':'RUNNING','started_utc':now(),'python':sys.version,'platform':platform.platform(),
            'dependencies':dependencies,
            'source_sha256':source_hashes,'steps':[],
            'scope':'Actual fresh execution of baseline plus revision checks; finite checks do not replace the all-n analytic proofs.',
            'evidence_types':['EXACT','NUMERIC','FORMULA_EVALUATION','PROVENANCE'],
            'scientific_checks':'RUNNING','import_closure':'UNVERIFIED'}
    report['dependency_deviations']={k:{'expected':v,'actual':report['dependencies'][k]} for k,v in expected.items() if report['dependencies'][k]!=v}
    env=os.environ.copy();env.pop('PYTHONHOME',None)
    env.pop('PYTHONPATH',None);env.pop('CP_IMPORT_AUDIT_DIR',None)
    # Apply to every descendant, including byte-preserved historical scripts.
    # UTF-8 avoids Windows locale failures in paths, JSON files and captured logs.
    env.update(PYTHONDONTWRITEBYTECODE='1',PYTHONUTF8='1',PYTHONIOENCODING='utf-8',
               OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
    copytree(ROOT/'revision_benchmarks',dest/'r1')
    specs=[('AUTHOR','MIXED',['run_all.py','--destination',str(dest/'author')],ROOT,6)]
    for name,args,kind in [
        ('exact_certificate.py',[],'EXACT'),
        ('certify_lz_n4_upper.py',['--verify-only'],'EXACT'),
        ('certify_fullmenu_half.py',['--verify-only'],'EXACT'),
        ('check_identities.py',[],'MIXED'),
        ('count_certificate_rays.py',[],'EXACT'),
        ('verify_n4_primal.py',[],'NUMERIC')]:
        specs.append(('R1',kind,[name,*args],dest/'r1',1))
    for nr,count in [(2,13),(3,11),(4,21)]:
        specs.append((f'R{nr}','MIXED',[f'fixed_menu_width/r{nr}/reproduce.py','--destination',str(dest/f'r{nr}')],ROOT,count))
    specs.append(('TABLES','FORMULA_EVALUATION',['manuscript_tables/generate_tables.py','--output',str(dest/'tables')],ROOT,1))
    # The supported pipeline has 11 direct and 60 descendant Python launches.
    # This explicit topology also catches an omitted whole descendant branch.
    expected_roots=[(cwd/cmd[0]).resolve() for _,_,cmd,cwd,_ in specs]
    expected_children={
        ROOT/'run_all.py':6,
        dest/'author/original_reproduction/reproduce.py':1,
        ROOT/'fixed_menu_width/r2/reproduce.py':13,
        ROOT/'fixed_menu_width/r3/reproduce.py':11,
        dest/'r3/code/audit_author.py':1,
        dest/'r3/sources/Supplementary_Software/run_all.py':6,
        dest/'r3/author_rerun/original_reproduction/reproduce.py':1,
        ROOT/'fixed_menu_width/r4/reproduce.py':21,
    }
    install_launcher(dest/'import_audit')
    started=time.monotonic()
    try:
        for i,(group,kind,cmd,cwd,count) in enumerate(specs,1):
            print(f'[{i}/{len(specs)}] {group}: {cmd[0]}',flush=True)
            start=time.monotonic()
            result=subprocess.run([sys.executable,*cmd],cwd=cwd,env=env,capture_output=True,text=True,encoding='utf-8')
            out=dest/'logs'/f'{i:02d}_{group}_stdout.txt';err=dest/'logs'/f'{i:02d}_{group}_stderr.txt'
            out.write_text(result.stdout,encoding='utf-8',newline='\n');err.write_text(result.stderr,encoding='utf-8',newline='\n')
            command=[x.replace(str(dest),'${RUN_ROOT}').replace(str(ROOT),'${PACKAGE_ROOT}') for x in cmd]
            report['steps'].append({'group':group,'evidence_type':kind,'command':['python',*command],
                                    'working_directory':'${PACKAGE_ROOT}' if cwd==ROOT else str(cwd.relative_to(dest)),
                                    'declared_child_stages':count,'exit_code':result.returncode,'seconds':time.monotonic()-start,
                                    'stdout':str(out.relative_to(dest)),'stdout_sha256':sha(out),'stderr':str(err.relative_to(dest)),'stderr_sha256':sha(err)})
            write(dest/'RUN_RECEIPT.json',report)
            if result.returncode:raise RuntimeError(f'{group} failed: {result.stderr[-2000:]}')
        child_receipts={}
        for group,rel,count in [('AUTHOR','author/run_report.json',6),('R2','r2/RUN_RECEIPT.json',13),('R3','r3/run_receipt.json',11),('R4','r4/run_receipt.json',21)]:
            data=json.loads((dest/rel).read_text(encoding='utf-8'));assert data.get('status')=='PASS' or data.get('passed') is True
            assert len(data['steps'])==count and all(x['exit_code']==0 for x in data['steps'])
            child_receipts[group]={'path':rel,'sha256':sha(dest/rel),'steps':count,'status':'PASS'}
        report['child_receipts']=child_receipts
        # Retain raw witness maxima and add a separately named gap-bound cap.
        bounds=[]
        from fractions import Fraction as F
        for path in sorted((dest/'r4/results').glob('robust_*_native_verified.json')):
            data=json.loads(path.read_text(encoding='utf-8'));raw=F(data['actual_maximum_rejection'])
            bounds.append({'n':data['n'],'exact_second_layer_size':data['adaptive_width'],
                           'raw_maximum_rejection':str(raw),'gap_upper_after_menu_half_cap':str(min(F(1,2),raw)),
                           'input_sha256':data['input_sha256'],'result':str(path.relative_to(dest)),
                           'coverage':'exact second-layer size only; at-most bounds require all smaller sizes or a separately proved class inclusion'})
        write(dest/'raw_witness_and_gap_bounds.json',bounds)
        for name in ['claim_values.json','width_rows.tex']:
            assert (dest/'tables'/name).read_bytes()==(ROOT/'manuscript_tables'/name).read_bytes(),name
        report['manuscript_table_bytes_match']=True
        # Scientific imports must resolve inside this package/run; third-party
        # libraries may resolve inside the chosen Python environment.
        runtime_roots={Path(sys.base_prefix).resolve(),Path(sys.prefix).resolve()}
        runtime_roots.update(Path(v).resolve() for k,v in sysconfig.get_paths().items() if k in ['stdlib','platstdlib','purelib','platlib'])
        report['scientific_checks']='PASS'
        audit=closure(dest/'import_audit',ROOT,dest,runtime_roots,expected_roots,expected_children)
        write(dest/'IMPORT_CLOSURE.json',audit)
        report['import_closure']=audit['status']
        if audit['status']!='PASS':raise RuntimeError('Import provenance failed: '+repr(audit['errors']))
        assert source_hashes==hashes(ROOT),'Supplied package changed during execution.'
        report.update(status='PASS',finished_utc=now(),seconds=time.monotonic()-started,source_bytes_unchanged=True,
                      import_closure='PASS',groups={'AUTHOR':6,'R1':6,'R2':13,'R3':11,'R4':21,'TABLES':1},
                      repetitions='R3 includes a repeated six-stage historical-author pipeline; repeated checks are not independent mathematical evidence.')
    except BaseException as exc:
        if report['scientific_checks']=='RUNNING':report['scientific_checks']='FAIL'
        report.update(status='INTERRUPTED' if isinstance(exc,(KeyboardInterrupt,SystemExit)) else 'FAIL',error=repr(exc),finished_utc=now(),seconds=time.monotonic()-started)
        write(dest/'RUN_RECEIPT.json',report);raise
    # Hash generated artifacts before writing the outer receipt/manifest to avoid self-reference.
    report['result_sha256']={str(p.relative_to(dest)):sha(p) for p in sorted(dest.rglob('*')) if p.is_file() and p.name!='RUN_RECEIPT.json' and '__pycache__' not in p.parts}
    write(dest/'RUN_RECEIPT.json',report)
    print(f'PASS: baseline 6; R1 6; R2 13; R3 11; R4 21; tables 1. Total {report["seconds"]:.2f}s',flush=True)
    return 0
if __name__=='__main__':raise SystemExit(main())

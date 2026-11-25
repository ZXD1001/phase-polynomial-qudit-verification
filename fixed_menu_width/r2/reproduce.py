#!/usr/bin/env python3
"""Reproduce the certified round-2 results in a NEW destination directory.

This entry point does not run exploratory optimizers or the unfinished n=7
proposal. It preserves the packaged source/input bytes and never writes to
an author's manuscript directory. Run without Python -O/-OO.
"""
from __future__ import annotations
import argparse,hashlib,json,os,platform,shutil,subprocess,sys,time
from pathlib import Path
from fractions import Fraction

ROOT=Path(__file__).resolve().parent

def digest(p: Path)->str:
    return hashlib.sha256(p.read_bytes()).hexdigest()

def main()->None:
    if not __debug__:
        raise RuntimeError('Run without Python -O or -OO.')
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--destination',type=Path,required=True)
    args=ap.parse_args();dest=args.destination.expanduser().resolve()
    if dest.exists():raise FileExistsError(f'Destination must be new: {dest}')
    if ROOT==dest or ROOT in dest.parents:
        raise ValueError('Use a sibling/outside destination, not a directory inside this package.')
    cases=[('witness_n4_w1','29/70'),('witness_n5_w1','26915/68048'),('witness_n6_w2','10049/22500'),('witness_n6_w1_companion','1319/3000')]
    files=list((ROOT/'code').glob('*.py'))+[ROOT/'inputs'/f'{name}.json' for name,_ in cases]
    files += [ROOT/'prior_checks'/name for name in ['exact_certificate.py','certify_fullmenu_half.py','lz_fullmenu_exact_half_certificate.json']]
    before={str(p.relative_to(ROOT)):digest(p) for p in files}
    dest.mkdir(parents=True)
    for d in ['code','inputs','results','logs','prior_checks']:(dest/d).mkdir()
    for p in files:shutil.copy2(p,dest/p.relative_to(ROOT))
    env=os.environ.copy();env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
    cmds=[]
    for native in [True,False]:
        for name,_ in cases:
            cmd=[sys.executable,'code/verify_width_certificate.py',f'inputs/{name}.json']
            if native:cmd.append('--native')
            cmds.append(cmd)
    cmds.extend([[sys.executable,'code/verify_direct_omega6_n4.py'],[sys.executable,'code/verify_achievable_tradeoff.py'],
      [sys.executable,'code/verify_structural_identities.py'],[sys.executable,'prior_checks/exact_certificate.py'],
      [sys.executable,'prior_checks/certify_fullmenu_half.py','--verify-only']])
    receipt={'status':'RUNNING','started_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'python':sys.version,
      'platform':platform.platform(),'source_hashes':before,'steps':[],'scope':'Round-2 exact finite certificates, analytic-construction checks and two round-1 exact certificate reruns; NOT the original authors complete software.'}
    try:
        for i,cmd in enumerate(cmds,1):
            # A slower computer must be allowed to finish the exact check.
            # The reviewer can still interrupt the run with Ctrl+C.
            t=time.time();p=subprocess.run(cmd,cwd=dest,env=env,text=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
            log=dest/'logs'/f'{i:02d}.log';log.write_text(p.stdout+'\nSTDERR\n'+p.stderr)
            receipt['steps'].append({'command':cmd[1:],'exit_code':p.returncode,'seconds':time.time()-t,'log':str(log.relative_to(dest))})
            print(f'{i}/{len(cmds)} exit={p.returncode} seconds={time.time()-t:.3f} {cmd[1:]}',flush=True)
            if p.returncode:raise RuntimeError(f'Stage {i} failed; inspect {log}')
        checks=[]
        for name,expected in cases:
            a=json.loads((dest/'results'/f'{name}_native_verified.json').read_text())
            b=json.loads((dest/'results'/f'{name}_reduced_verified.json').read_text())
            if a['actual_maximum_rejection']!=expected:raise AssertionError((name,expected,a['actual_maximum_rejection']))
            aa=a['per_template_maximum_numerators'];bb=b['per_template_maximum_numerators']
            same=len(aa)==len(bb) and all(Fraction(x,a['common_denominator'])==Fraction(y,b['common_denominator']) for x,y in zip(aa,bb))
            if not same:raise AssertionError(f'Native/reduced template mismatch: {name}')
            checks.append({'case':name,'expected_exact_upper':expected,'all_template_values_identical':True,'templates':len(aa)})
        receipt['crosschecks']=checks
        receipt['at_most_width_coverage']={'n4_width1':'all size1; size0 embeds','n5_width1':'all size1; size0 embeds','n6_width2':'all size1 and size2 checked separately; size0 embeds','n6_combined_upper':'10049/22500'}
        if before!={str(p.relative_to(ROOT)):digest(p) for p in files}:raise AssertionError('Packaged source/input bytes changed')
        receipt['source_inputs_preserved']=True;receipt['status']='PASS'
    except Exception as exc:
        receipt['status']='FAIL';receipt['error']=repr(exc)
        raise
    finally:
        receipt['finished_utc']=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())
        (dest/'RUN_RECEIPT.json').write_text(json.dumps(receipt,indent=2))
    print('PASS: exact values, native/reduced agreement, direct omega6 check, construction checks, and prior certificate reruns.',flush=True)

if __name__=='__main__':main()

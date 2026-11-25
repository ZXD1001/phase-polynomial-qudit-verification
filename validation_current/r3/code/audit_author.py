#!/usr/bin/env python3
"""Verify supplied source manifest, then run the author's unmodified six stages."""
from __future__ import annotations
from pathlib import Path
import argparse,hashlib,json,subprocess,sys,time
if not __debug__:raise RuntimeError('Run without -O/-OO.')

def check(root):
    result=[]
    for line in (root/'SHA256SUMS.txt').read_text().splitlines():
        if not line.strip():continue
        expected,name=line.split(maxsplit=1);name=name.lstrip('*');path=(root/name).resolve()
        if not path.is_relative_to(root.resolve()):raise ValueError('Unsafe manifest path')
        actual=hashlib.sha256(path.read_bytes()).hexdigest()
        if actual!=expected:raise ValueError(f'Hash mismatch: {name}')
        result.append({'file':name,'sha256':actual})
    return result

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source',type=Path,required=True)
    ap.add_argument('--destination',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    a=ap.parse_args();tic=time.monotonic();before=check(a.source)
    cmd=[sys.executable,str(a.source/'run_all.py'),'--destination',str(a.destination)]
    subprocess.run(cmd,check=True)
    after=check(a.source);assert before==after
    author=json.loads((a.destination/'run_report.json').read_text())
    assert author['passed'] and author['source_bytes_unchanged']
    assert len(author['steps'])==6 and all(x['exit_code']==0 for x in author['steps'])
    result={'status':'PASS','manifest_entries':len(before),'all_source_manifest_entries_unchanged':True,
            'author_six_stage_rerun_passed':True,'author_run_report':str(a.destination/'run_report.json'),
            'seconds':time.monotonic()-tic,'manifest':before}
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='manifest'},indent=2))
if __name__=='__main__':main()

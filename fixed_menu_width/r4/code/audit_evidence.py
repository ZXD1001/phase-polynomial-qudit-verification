#!/usr/bin/env python3
"""Check the received author archive and prior receipts without claiming a new full rerun."""
from pathlib import Path
import hashlib,json,zipfile,argparse
ROOT=Path(__file__).resolve().parents[1]
if not __debug__:raise RuntimeError('Run without -O/-OO.')
def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    zpath=ROOT/'source_archives/Supplementary_Software.zip';rows=[]
    with zipfile.ZipFile(zpath) as z:
        assert z.testzip() is None
        manifest=next(n for n in z.namelist() if n.endswith('/SHA256SUMS.txt') and '__MACOSX' not in n)
        base=manifest[:-len('SHA256SUMS.txt')]
        for line in z.read(manifest).decode().splitlines():
            if not line.strip():continue
            digest,name=line.split(maxsplit=1);name=name.strip().lstrip('*')
            actual=hashlib.sha256(z.read(base+name)).hexdigest();assert actual==digest,(name,digest,actual)
            rows.append(dict(file=name,sha256=actual))
    prior=json.loads((ROOT/'reference_receipts/round3_validation_receipt.json').read_text())
    author=json.loads((ROOT/'reference_receipts/author_round3_run_report.json').read_text())
    assert prior['status']=='PASS' and all(s['exit_code']==0 for s in prior['steps'])
    assert author['passed'] and len(author['steps'])==6 and all(s['exit_code']==0 for s in author['steps'])
    with zipfile.ZipFile(zpath) as z:
        for name,digest in author['source_sha256'].items():assert hashlib.sha256(z.read(base+name)).hexdigest()==digest
    out=dict(status='PASS',author_archive_sha256=hashlib.sha256(zpath.read_bytes()).hexdigest(),
        author_manifest_entries=len(rows),manifest=rows,prior_round3_steps=len(prior['steps']),
        prior_author_six_stage_receipt_matches_received_source=True,
        author_full_pipeline_rerun_this_round=False,
        note='Archive bytes/manifests and prior receipts checked now. The six-stage run itself occurred in round 3, not this round.')
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k!='manifest'},indent=2))
if __name__=='__main__':main()

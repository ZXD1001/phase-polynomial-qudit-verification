"""Optional import provenance observer used by run_revision_checks.py.

No scientific functions or import resolution are patched. Each subprocess
records loaded module filenames at normal termination for closure auditing.
"""
import atexit
import json
import os
from pathlib import Path
import sys

def record():
    dest=os.environ.get('CP_IMPORT_AUDIT_DIR')
    if not dest:return
    rows={}
    for name,module in list(sys.modules.items()):
        f=getattr(module,'__file__',None)
        if f:
            rows[name]=str(Path(f).resolve())
    Path(dest).mkdir(parents=True,exist_ok=True)
    (Path(dest)/f'{os.getpid()}.json').write_text(json.dumps({'modules':rows},indent=2)+'\n')
atexit.register(record)

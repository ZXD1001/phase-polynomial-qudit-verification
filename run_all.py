#!/usr/bin/env python3
"""Reproduce all calculations in a NEW directory, preserving delivered data."""
from __future__ import annotations
import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path
import platform
import stat
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
SOURCES = [
    'original_reproduction/reproduce.py', 'original_reproduction/inputs.json',
    'original_reproduction/provenance.json',
    'original_reproduction/historical_original/validate_phase_polynomial_blocks.py',
    'lz_benchmark/benchmark_lz.py', 'lz_benchmark/verify_exports.py',
    'structured_benchmark/run_benchmark.py', 'structured_benchmark/pauli_search.py',
    'structured_benchmark/exact_n2_certificate.py',
    'resources/verify_resource_model.py', 'summarize_results.py',
]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    stamp = datetime.datetime.now().strftime('%Y%m%d-%H%M%S')
    parser.add_argument('--destination', type=Path, default=ROOT.parent/f'reproduction-rerun-{stamp}')
    args = parser.parse_args()
    dest = args.destination.resolve()
    if dest.exists():
        raise SystemExit('Destination must not already exist; recorded results will not be overwritten.')
    dest.mkdir(parents=True)
    before = {name: digest(ROOT/name) for name in SOURCES}
    for name in SOURCES:
        target = dest/name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes((ROOT/name).read_bytes())
    commands = [
        ['original_reproduction/reproduce.py'],
        ['lz_benchmark/benchmark_lz.py', '--n', '2', '3', '4'],
        ['lz_benchmark/verify_exports.py'],
        ['structured_benchmark/run_benchmark.py'],
        ['resources/verify_resource_model.py'],
        ['summarize_results.py'],
    ]
    report = {'python': sys.version, 'platform': platform.platform(),
              'source_sha256': before, 'steps': [], 'destination': str(dest)}
    (dest/'logs').mkdir()
    env = os.environ.copy()
    env['PYTHONDONTWRITEBYTECODE'] = '1'
    for index, command in enumerate(commands, 1):
        print(f'[{index}/{len(commands)}] {command[0]}', flush=True)
        invocation = [sys.executable, *[str(dest/command[0]), *command[1:]]]
        run = subprocess.run(invocation, cwd=dest, capture_output=True, text=True, env=env)
        (dest/'logs'/f'{index:02d}_stdout.txt').write_text(run.stdout)
        (dest/'logs'/f'{index:02d}_stderr.txt').write_text(run.stderr)
        report['steps'].append({'command': command, 'exit_code': run.returncode})
        if run.returncode:
            report['passed'] = False
            (dest/'run_report.json').write_text(json.dumps(report, indent=2)+'\n')
            print(run.stderr, file=sys.stderr)
            return run.returncode
    assert before == {name: digest(ROOT/name) for name in SOURCES}
    report['source_bytes_unchanged'] = True
    report['passed'] = True
    (dest/'run_report.json').write_text(json.dumps(report, indent=2)+'\n')
    # New generated deliverables must also be Finder-visible on macOS.
    if sys.platform == 'darwin':
        subprocess.run(['chflags', '-R', 'nohidden', str(dest)], check=True)
        assert all(not (p.stat().st_flags & stat.UF_HIDDEN) for p in [dest, *dest.rglob('*')])
    print(f'PASS: all calculations reproduced in {dest}', flush=True)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())

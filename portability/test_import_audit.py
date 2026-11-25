"""Small regression tests for observer startup and fail-closed process coverage.

Run directly with the same Python interpreter used for the scientific pipeline.
These tests use disposable scripts and never execute or alter scientific data.
"""
from __future__ import annotations
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import sysconfig
import tempfile
import unittest

sys.dont_write_bytecode = True
PACKAGE = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('audit_under_test', PACKAGE/'portability/import_audit.py')
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)


class ImportAuditTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix='cp observer 审稿测试 ')
        self.run = Path(self.temporary.name).resolve()
        self.directory = self.run/'import_audit'
        self.directory.mkdir()
        (self.run/'leaf.py').write_text('import fractions\nassert fractions.Fraction(2, 4) == fractions.Fraction(1, 2)\n', encoding='utf-8')
        (self.run/'grandchild.py').write_text('import decimal\nassert decimal.Decimal("0.5") * 2 == 1\n', encoding='utf-8')
        (self.run/'branch.py').write_text('import subprocess, sys\nsubprocess.run([sys.executable, "grandchild.py"], check=True)\n', encoding='utf-8')
        (self.run/'parent.py').write_text('import subprocess, sys\nfor child in ("branch.py", "leaf.py"):\n    subprocess.run([sys.executable, child], check=True)\n', encoding='utf-8')
        self.env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1',
                        PYTHONUTF8='1', PYTHONIOENCODING='utf-8')
        self.env.pop('CP_IMPORT_AUDIT_DIR', None)
        self.env.pop('PYTHONPATH', None)
        self.runtime = {Path(sys.base_prefix), Path(sys.prefix)}
        self.runtime.update(Path(v) for k, v in sysconfig.get_paths().items() if k in ('stdlib', 'platstdlib', 'purelib', 'platlib'))

    def tearDown(self):
        self.temporary.cleanup()

    def execute(self, target='parent.py', bypass=False, environment=None):
        # The separate driver avoids persisting Python audit hooks between tests.
        driver = self.run/'driver.py'
        driver.write_text('import sys, subprocess\n'
                          f'sys.path.insert(0, {str(PACKAGE/"portability")!r})\n'
                          'from import_audit import install_launcher, ORIGINAL_POPEN\n'
                          f'install_launcher({str(self.directory)!r})\n'
                          + (f'proc=ORIGINAL_POPEN([sys.executable, {target!r}]); assert proc.wait()==0\n' if bypass else
                             f'subprocess.run([sys.executable, {target!r}], check=True)\n'), encoding='utf-8')
        result = subprocess.run([sys.executable, str(driver)], cwd=self.run,
                                env=environment or self.env, capture_output=True,
                                text=True, encoding='utf-8')
        self.assertEqual(result.returncode, 0, result.stderr)

    def result(self, target='parent.py', children=None):
        if children is None:
            children = {self.run/'parent.py': 2, self.run/'branch.py': 1}
        return audit.closure(self.directory, PACKAGE, self.run, self.runtime,
                             [self.run/target], children)

    def test_normal_child_and_grandchild_coverage(self):
        self.execute()
        result = self.result()
        self.assertEqual(result['status'], 'PASS', result)
        self.assertEqual(result['launched_process_count'], 4)
        self.assertEqual(result['observed_process_count'], 4)
        self.assertGreater(result['observed_module_file_count'], 0)

    def test_empty_observation_directory_never_passes(self):
        result = self.result()
        self.assertEqual(result['status'], 'FAIL')
        self.assertEqual(result['observed_process_count'], 0)

    def test_observer_absent_bypass_is_detected(self):
        self.execute('leaf.py', bypass=True)
        result = self.result('leaf.py', {})
        self.assertEqual(result['status'], 'FAIL')
        self.assertEqual(result['observed_process_count'], 0)
        self.assertEqual(len(result['unwrapped_launch_records']), 1)

    def test_one_missing_child_completion_never_passes(self):
        self.execute()
        missing = next(p for p in self.directory.glob('*.observed.json')
                       if Path(json.loads(p.read_text(encoding='utf-8'))['entrypoint']).name == 'grandchild.py')
        missing.unlink()
        result = self.result()
        self.assertEqual(result['status'], 'FAIL')
        self.assertEqual(result['observed_process_count'], 3)
        self.assertEqual(len(result['missing_invocations']), 1)

    def test_missing_entire_child_registration_never_passes(self):
        self.execute()
        launch = next(p for p in self.directory.glob('*.launch.json')
                      if Path(json.loads(p.read_text(encoding='utf-8'))['entrypoint']).name == 'grandchild.py')
        invocation = json.loads(launch.read_text(encoding='utf-8'))['invocation']
        launch.unlink()
        (self.directory/(invocation+'.observed.json')).unlink()
        result = self.result()
        self.assertEqual(result['status'], 'FAIL')
        self.assertTrue(any('Child coverage mismatch' in error for error in result['errors']))

    def test_unregistered_extra_completion_never_passes(self):
        self.execute()
        observation = json.loads(next(self.directory.glob('*.observed.json')).read_text(encoding='utf-8'))
        observation['invocation'] = 'unregistered-test-invocation'
        (self.directory/'unregistered-test-invocation.observed.json').write_text(json.dumps(observation), encoding='utf-8')
        result = self.result()
        self.assertEqual(result['status'], 'FAIL')
        self.assertEqual(len(result['unregistered_invocations']), 1)

    def test_empty_module_snapshots_never_pass(self):
        self.execute()
        for path in self.directory.glob('*.observed.json'):
            observation = json.loads(path.read_text(encoding='utf-8'))
            observation['modules'] = {}
            path.write_text(json.dumps(observation), encoding='utf-8')
        result = self.result()
        self.assertEqual(result['status'], 'FAIL')
        self.assertTrue(any('Empty or invalid observed modules' in error for error in result['errors']))

    def test_preloaded_sitecustomize_cannot_displace_explicit_observer(self):
        conflicting = self.run/'startup'
        conflicting.mkdir()
        (conflicting/'sitecustomize.py').write_text('"""Unrelated preloaded sitecustomize; contains no CP observer."""\n', encoding='utf-8')
        environment = dict(self.env, PYTHONPATH=str(conflicting))
        self.execute(environment=environment)
        result = self.result()
        self.assertEqual(result['status'], 'PASS', result)
        self.assertEqual(result['observed_process_count'], 4)
        for path in self.directory.glob('*.observed.json'):
            self.assertEqual(json.loads(path.read_text(encoding='utf-8'))['modules']['sitecustomize'], str(conflicting/'sitecustomize.py'))

    def test_unexpected_external_sitecustomize_is_not_whitelisted(self):
        with tempfile.TemporaryDirectory(prefix='cp-external-startup-') as outside:
            conflicting = Path(outside).resolve()
            (conflicting/'sitecustomize.py').write_text('"""External startup module, deliberately outside permitted roots."""\n', encoding='utf-8')
            self.execute(environment=dict(self.env, PYTHONPATH=str(conflicting)))
            result = self.result()
            self.assertEqual(result['status'], 'FAIL')
            self.assertEqual(result['observed_process_count'], 4)
            self.assertIn(str(conflicting/'sitecustomize.py'), result['unexpected_module_files'])

    def test_abrupt_child_exit_has_no_completion_and_fails(self):
        (self.run/'abrupt.py').write_text('import os\nos._exit(0)\n', encoding='utf-8')
        self.execute('abrupt.py')
        result = self.result('abrupt.py', {})
        self.assertEqual(result['status'], 'FAIL')
        self.assertEqual(result['launched_process_count'], 1)
        self.assertEqual(result['observed_process_count'], 0)


if __name__ == '__main__':
    unittest.main(verbosity=2)

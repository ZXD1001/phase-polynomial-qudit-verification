"""Small regressions for reviewer startup and platform-independent table bytes.

These tests use temporary output only. They do not run scientific simulations
or require NumPy, SciPy or SymPy to be installed.
"""
from __future__ import annotations
import contextlib
import io
import os
from pathlib import Path
import runpy
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

sys.dont_write_bytecode = True
PACKAGE = Path(__file__).resolve().parents[1]


class ReviewerPortabilityTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix='cp reviewer 审稿测试 ')
        self.run = Path(self.temporary.name).resolve()
        self.env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1',
                        PYTHONUTF8='1', PYTHONIOENCODING='utf-8')
        for name in ('PYTHONPATH', 'PYTHONHOME', 'CP_IMPORT_AUDIT_DIR'):
            self.env.pop(name, None)

    def tearDown(self):
        self.temporary.cleanup()

    def preflight_without_site_packages(self, destination):
        # -S really removes site-package discovery; no dependency API is mocked.
        return subprocess.run(
            [sys.executable, '-S', str(PACKAGE/'run_revision_checks.py'),
             '--destination', str(destination)],
            cwd=self.run, env=self.env, capture_output=True,
            text=True, encoding='utf-8')

    def test_missing_dependencies_explain_install_without_creating_output(self):
        destination = self.run/'new results'
        result = self.preflight_without_site_packages(destination)
        self.assertEqual(result.returncode, 2, result.stderr)
        self.assertIn('Missing dependencies: numpy, scipy, sympy', result.stderr)
        self.assertIn('same Python interpreter', result.stderr)
        self.assertIn(sys.executable, result.stderr)
        self.assertIn('pip', result.stderr)
        self.assertIn('requirements.txt', result.stderr)
        self.assertNotIn('Traceback', result.stderr)
        self.assertFalse(destination.exists())

    def test_existing_output_is_preserved_with_clear_error(self):
        destination = self.run/'existing results'
        destination.mkdir()
        sentinel = destination/'keep.txt'
        sentinel.write_bytes(b'Existing reviewer results must survive.\n')
        before = sentinel.read_bytes()
        result = self.preflight_without_site_packages(destination)
        self.assertEqual(result.returncode, 2, result.stderr)
        self.assertIn('output directory already exists', result.stderr)
        self.assertIn('never overwritten', result.stderr)
        self.assertNotIn('Traceback', result.stderr)
        self.assertEqual(list(destination.iterdir()), [sentinel])
        self.assertEqual(sentinel.read_bytes(), before)

    def test_actual_table_generator_matches_reference_under_windows_newlines(self):
        generator = PACKAGE/'manuscript_tables/generate_tables.py'
        namespace = runpy.run_path(str(generator))
        destination = self.run/'generated tables'
        original_open = Path.open

        def windows_text_open(path, mode='r', buffering=-1, encoding=None,
                              errors=None, newline=None):
            # Simulate the real platform difference at the text I/O boundary.
            # Explicit LF and binary I/O retain their normal semantics.
            if 'b' not in mode and any(flag in mode for flag in 'wax') and newline is None:
                newline = '\r\n'
            return original_open(path, mode, buffering, encoding, errors, newline)

        with mock.patch.object(Path, 'open', windows_text_open):
            probe = self.run/'default newline probe.txt'
            probe.write_text('first\nsecond\n', encoding='utf-8')
            self.assertEqual(probe.read_bytes(), b'first\r\nsecond\r\n')
            with mock.patch.object(sys, 'argv', [str(generator), '--output', str(destination)]):
                with contextlib.redirect_stdout(io.StringIO()):
                    self.assertEqual(namespace['main'](), 0)

        for name in ('claim_values.json', 'width_rows.tex'):
            self.assertEqual((destination/name).read_bytes(),
                             (PACKAGE/'manuscript_tables'/name).read_bytes(), name)
        self.assertNotIn(b'\r\n', (destination/'width_intervals.tex').read_bytes())


if __name__ == '__main__':
    unittest.main(verbosity=2)

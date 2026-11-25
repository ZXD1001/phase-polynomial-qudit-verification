#!/usr/bin/env python3
"""Reproduce manuscript Eq. (38)/(39) and label d=2,3 as regressions.

The preserved historical script uses cyclic generator averages. This separate
implementation uses direct D-conjugation of tensor Fourier projectors and a
full branch measurement constructed from finite differences of the complete
polynomial. It neither imports nor edits the historical script.
"""
from __future__ import annotations

import argparse
import contextlib
import datetime
import hashlib
import io
import itertools
import json
import os
from pathlib import Path
import platform
import re
import subprocess
import sys
import warnings

import numpy as np

ROOT = Path(__file__).resolve().parent


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, data):
    Path(path).write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n")


def opnorm(matrix):
    return float(np.linalg.norm(matrix, 2))


def tensor(matrices):
    result = np.ones((1, 1), dtype=complex)
    for matrix in matrices:
        result = np.kron(result, matrix)
    return result


def environment():
    config = io.StringIO()
    with warnings.catch_warnings(record=True) as observed_warnings:
        with contextlib.redirect_stdout(config):
            np.show_config()
    return {"python_version": platform.python_version(),
            "python_executable": sys.executable, "numpy_version": np.__version__,
            "platform": platform.platform(), "machine": platform.machine(),
            "dtype": "complex128", "numpy_configuration": config.getvalue(),
            "numpy_configuration_warnings": [str(w.message) for w in observed_warnings],
            "thread_environment": {k: os.environ.get(k) for k in
                ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                 "VECLIB_MAXIMUM_THREADS")}}


def independent_case(case, tolerance):
    d, n = case["dimension"], case["parties"]
    blocks = case["blocks"]
    if sorted(v for block in blocks for v in block) != list(range(n)):
        raise ValueError("The supplied blocks must form a nonempty color partition")
    terms = [(term["coefficient"], {int(v): s for v, s in term["powers"].items()})
             for term in case["terms"]]
    for block in blocks:
        if not block or any(len(set(block) & set(powers)) > 1
                            for coefficient, powers in terms if coefficient % d):
            raise ValueError("A supplied block is not independent")

    def polynomial(values):
        total = 0
        for coefficient, powers in terms:
            term = coefficient % d
            for vertex, exponent in powers.items():
                term = term * pow(int(values[vertex]), exponent, d) % d
            total = (total + term) % d
        return total

    states = list(itertools.product(range(d), repeat=n))
    index = {state: i for i, state in enumerate(states)}
    size = d**n
    roots = np.exp(2j*np.pi*np.arange(d)/d)
    diagonal = roots[[polynomial(state) for state in states]]
    target = diagonal/np.sqrt(size)
    x = np.arange(d)
    fourier = roots[(-x[:, None]*x[None, :]) % d]/np.sqrt(d)
    pvm_factors = [np.outer(fourier[:, j], fourier[:, j].conj()) for j in range(d)]
    identity = np.eye(size, dtype=complex)
    accepted_blocks, block_results = [], []

    for block in blocks:
        complement = [v for v in range(n) if v not in block]
        labels = list(itertools.product(range(d), repeat=len(block)))
        records = list(itertools.product(range(d), repeat=len(complement)))
        measurement = np.zeros((size, size), dtype=complex)
        column = 0
        factorization_checks = 0
        max_local_gram_error = 0.0
        for record in records:
            base = [0]*n
            for vertex, value in zip(complement, record, strict=True):
                base[vertex] = value
            base_phase = polynomial(base)
            local_phases = {}
            for vertex in block:
                phase = []
                for value in range(d):
                    varied = base.copy()
                    varied[vertex] = value
                    phase.append((polynomial(varied)-base_phase) % d)
                local_phases[vertex] = phase
                local_basis = roots[np.array(phase)[:, None]]*fourier
                max_local_gram_error = max(max_local_gram_error,
                    opnorm(local_basis.conj().T@local_basis-np.eye(d)))
            for terminal_values in labels:
                values = base.copy()
                for v, value in zip(block, terminal_values, strict=True):
                    values[v] = value
                reconstructed = (base_phase + sum(local_phases[v][values[v]] for v in block)) % d
                if reconstructed != polynomial(values):
                    raise ValueError("Exact conditional factorization failed")
                factorization_checks += 1
            for sector in labels:
                for terminal_values in labels:
                    values = base.copy()
                    for v, value in zip(block, terminal_values, strict=True):
                        values[v] = value
                    phase = sum(local_phases[v][values[v]]-j*values[v]
                                for v, j in zip(block, sector, strict=True)) % d
                    measurement[index[tuple(values)], column] = roots[phase]/np.sqrt(d**len(block))
                column += 1

        completeness = np.zeros((size, size), dtype=complex)
        sector_residuals = []
        for label_index, sector in enumerate(labels):
            selected = measurement[:, label_index::len(labels)]
            effect = selected@selected.conj().T
            completeness += effect
            sector_map = dict(zip(block, sector, strict=True))
            before_D = tensor([pvm_factors[sector_map[v]] if v in sector_map else np.eye(d)
                               for v in range(n)])
            ideal = diagonal[:, None]*before_D*diagonal.conj()[None, :]
            sector_residuals.append(opnorm(effect-ideal))
            if all(j == 0 for j in sector):
                accepted = effect
                accepted_ideal = ideal
        accepted_blocks.append(accepted)
        errors = {
            "local_basis_gram_operator_norm": max_local_gram_error,
            "full_branch_gram_operator_norm": opnorm(measurement.conj().T@measurement-identity),
            "full_measurement_completeness_operator_norm": opnorm(completeness-identity),
            "all_sector_effect_vs_direct_D_projector_max_operator_norm": max(sector_residuals),
            "accept_effect_vs_direct_D_projector_operator_norm": opnorm(accepted-accepted_ideal),
            "accept_effect_idempotence_operator_norm": opnorm(accepted@accepted-accepted),
            "accept_effect_hermiticity_operator_norm": opnorm(accepted-accepted.conj().T),
            "target_acceptance_vector_norm": float(np.linalg.norm(accepted@target-target))}
        block_results.append({"block": block, "complement": complement,
            "complement_records": len(records), "joint_sector_outcomes": len(labels),
            "complete_rank_one_branches": size, "exact_modular_factorization_cases": factorization_checks,
            "accepted_projector_expected_rank": d**len(complement),
            "accepted_effect_trace": float(np.trace(accepted).real),
            "errors": errors, "passed": all(v < tolerance for v in errors.values())})

    verification = sum(accepted_blocks)/len(blocks)
    eigenvalues = np.linalg.eigvalsh(verification)
    expected_eigenvalues = sorted(sum(all(state[v] == 0 for v in block) for block in blocks)/len(blocks)
                                  for state in states)
    spectral_error = float(np.max(np.abs(eigenvalues-expected_eigenvalues)))
    observed_gap = float(1-eigenvalues[-2])
    spectrum = [{"value": value, "expected_multiplicity": expected_eigenvalues.count(value),
                 "observed_multiplicity_within_tolerance": int(np.sum(np.abs(eigenvalues-value)<tolerance))}
                for value in sorted(set(expected_eigenvalues))]
    global_errors = {"complete_spectrum_max_absolute_error": spectral_error,
        "gap_absolute_error": abs(observed_gap-case["expected_gap"]),
        "target_acceptance_vector_norm": float(np.linalg.norm(verification@target-target)),
        "verification_hermiticity_operator_norm": opnorm(verification-verification.conj().T)}
    return {"case_id": case["id"], "role": case["role"], "dimension": d, "parties": n,
        "hilbert_dimension": size, "blocks": block_results,
        "maximum_accept_effect_residual": max(b["errors"]["accept_effect_vs_direct_D_projector_operator_norm"] for b in block_results),
        "maximum_all_sector_residual": max(b["errors"]["all_sector_effect_vs_direct_D_projector_max_operator_norm"] for b in block_results),
        "expected_gap": case["expected_gap"], "observed_gap": observed_gap,
        "spectrum": spectrum, "global_errors": global_errors,
        "passed": all(b["passed"] for b in block_results) and all(v < tolerance for v in global_errors.values())}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=ROOT/"inputs.json")
    parser.add_argument("--output-dir", type=Path, default=ROOT/"results")
    args = parser.parse_args()
    started = datetime.datetime.now(datetime.timezone.utc).isoformat()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    spec = json.loads(args.input.read_text())
    tolerance = spec["absolute_tolerance"]
    if tolerance <= 0:
        raise ValueError("Tolerance must be positive")
    historical = ROOT/"historical_original"/"validate_phase_polynomial_blocks.py"
    provenance = json.loads((ROOT/"provenance.json").read_text())
    expected_hash = next(item["sha256"] for item in provenance["files"]
                         if item.get("preserved_path") == "historical_original/validate_phase_polynomial_blocks.py")
    if sha256(historical) != expected_hash:
        raise ValueError("The preserved historical script does not match provenance.json")
    old_run = subprocess.run([sys.executable, str(historical)], cwd=historical.parent,
                             capture_output=True, text=True, check=False)
    (args.output_dir/"historical_stdout.txt").write_text(old_run.stdout)
    (args.output_dir/"historical_stderr.txt").write_text(old_run.stderr)
    parsed = [{"dimension": int(d), "operator_norm_residual_printed": float(r),
               "gap_printed": float(g)} for d, r, g in
              re.findall(r"d=(\d+): residual=([\d.eE+-]+), gap=([\d.eE+-]+)", old_run.stdout)]
    env = environment()
    write_json(args.output_dir/"environment.json", env)
    results = [independent_case(case, tolerance) for case in spec["cases"]]
    passed = old_run.returncode == 0 and len(parsed) == 4 and all(case["passed"] for case in results)
    report = {"schema_version": 1, "run_started_utc": started,
        "run_completed_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "scope": "Original Eq. (38)/(39) and d=2,3 regression matrices; not all manuscript theorems or new six-level path candidate.",
        "input_sha256": sha256(args.input), "runner_sha256": sha256(__file__),
        "absolute_tolerance": tolerance, "dtype": "complex128",
        "historical_run": {"script_sha256": sha256(historical), "exit_code": old_run.returncode,
            "historical_assertion_tolerance": spec["historical_absolute_tolerance"],
            "printed_results": parsed,
            "note": "The unmodified script prints residuals to three digits after the decimal in scientific notation."},
        "independent_results": results, "passed": passed}
    write_json(args.output_dir/"reproduction_results.json", report)
    for result in results:
        print(f"{result['case_id']}: pass={result['passed']}, "
              f"accept_residual={result['maximum_accept_effect_residual']:.6e}, "
              f"all_sector_residual={result['maximum_all_sector_residual']:.6e}, "
              f"gap={result['observed_gap']:.15f}")
    print(f"historical_exit_code={old_run.returncode}; overall_pass={passed}")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())

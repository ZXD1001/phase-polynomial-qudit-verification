"""Explicit subprocess provenance for the supported revision pipeline.

This instruments process launch/recording only, not scientific imports or math.
All supported Python commands are script invocations via sys.executable. An
unsupported Python invocation is rejected rather than silently left unobserved.
"""
from __future__ import annotations
import collections
import json
import os
from pathlib import Path
import subprocess
import sys
import uuid

BOOTSTRAP = Path(__file__).with_name("observed_python.py").resolve()
ORIGINAL_POPEN = subprocess.Popen


def write_json(path, value):
    path = Path(path)
    temporary = path.with_suffix(".tmp")
    temporary.write_text(json.dumps(value, indent=2) + "\n")
    temporary.replace(path)


def install_launcher(directory, parent=None):
    """Register every child before launching it through the explicit observer."""
    directory = Path(directory).resolve()
    directory.mkdir(parents=True, exist_ok=True)
    interpreter = Path(sys.executable).resolve()

    def python_command(command):
        return isinstance(command, (list, tuple)) and command and Path(command[0]).resolve() == interpreter

    def audit_hook(event, args):
        if event != "subprocess.Popen":
            return
        executable, command, cwd, env = args
        if Path(executable).resolve() != interpreter:
            return
        # The audit event independently catches a direct use of an unwrapped
        # Popen object, even if the launcher below has been bypassed.
        if (not isinstance(command, (list, tuple)) or len(command) < 5
                or Path(command[1]).resolve() != BOOTSTRAP
                or command[2] != str(directory)
                or not (directory / (str(command[3]) + ".launch.json")).is_file()):
            write_json(directory / (uuid.uuid4().hex + ".unobserved.json"),
                       {"parent": parent, "command": list(command) if isinstance(command, (list, tuple)) else str(command),
                        "reason": "Python subprocess bypassed the explicit observer"})
    sys.addaudithook(audit_hook)

    class ObservedPopen(ORIGINAL_POPEN):
        def __init__(self, args, *positional, **kwargs):
            if python_command(args):
                if (positional or kwargs.get("shell") or kwargs.get("executable")
                        or len(args) < 2 or str(args[1]).startswith("-")):
                    raise RuntimeError("Unsupported Python invocation in the audited pipeline")
                invocation = uuid.uuid4().hex
                cwd = Path(kwargs.get("cwd") or os.getcwd()).resolve()
                entrypoint = (cwd / str(args[1])).resolve()
                if not entrypoint.is_file() or entrypoint.suffix != ".py":
                    raise RuntimeError("Audited Python command must name an existing .py script")
                launch = {"invocation": invocation, "parent": parent, "entrypoint": str(entrypoint),
                          "argv": [str(x) for x in args[1:]], "cwd": str(cwd)}
                launch_path = directory / (invocation + ".launch.json")
                write_json(launch_path, launch)
                command = [sys.executable, str(BOOTSTRAP), str(directory), invocation, *args[1:]]
                super().__init__(command, **kwargs)
                launch["pid"] = self.pid
                write_json(launch_path, launch)
            else:
                super().__init__(args, *positional, **kwargs)
    subprocess.Popen = ObservedPopen


def closure(directory, package_root, run_root, runtime_roots,
            expected_roots, expected_children):
    """Fail closed on absent/incomplete records, topology mismatch or bad paths."""
    directory, package_root, run_root = map(lambda p: Path(p).resolve(), (directory, package_root, run_root))
    errors, launches, observations = [], {}, {}
    for suffix, destination in [("launch", launches), ("observed", observations)]:
        for path in sorted(directory.glob("*." + suffix + ".json")):
            try:
                data = json.loads(path.read_text())
                key = data["invocation"]
                if key in destination or path.name != key + "." + suffix + ".json":
                    raise ValueError("duplicate or mismatched invocation identifier")
                destination[key] = data
            except (ValueError, KeyError, TypeError) as exc:
                errors.append(f"Invalid {suffix} record {path.name}: {exc}")
    if not launches:
        errors.append("No Python launch records were registered")
    if not observations:
        errors.append("No explicit observer completion records were produced")
    missing = sorted(set(launches) - set(observations))
    extra = sorted(set(observations) - set(launches))
    if missing:
        errors.append(f"Missing observer records for {len(missing)} registered Python invocations")
    if extra:
        errors.append(f"Unregistered observer records: {len(extra)}")
    unobserved = sorted(p.name for p in directory.glob("*.unobserved.json"))
    if unobserved:
        errors.append(f"Unwrapped Python subprocess launches: {len(unobserved)}")
    roots = collections.Counter(x["entrypoint"] for x in launches.values() if x.get("parent") is None)
    expected_roots = collections.Counter(str(Path(p).resolve()) for p in expected_roots)
    if roots != expected_roots:
        errors.append("Top-level Python invocation paths/counts do not match the declared pipeline")
    expected_children = {str(Path(p).resolve()): n for p, n in expected_children.items()}
    children = collections.Counter(x.get("parent") for x in launches.values() if x.get("parent") is not None)
    for invocation, launch in launches.items():
        parent = launch.get("parent")
        if parent is not None and parent not in launches:
            errors.append(f"Unknown parent for {invocation}")
        wanted = expected_children.get(launch["entrypoint"], 0)
        if children[invocation] != wanted:
            errors.append(f"Child coverage mismatch for {launch['entrypoint']}: expected {wanted}, observed {children[invocation]}")
    absent_parents = set(expected_children) - {x["entrypoint"] for x in launches.values()}
    if absent_parents:
        errors.append("Declared parent entrypoints were not executed: " + repr(sorted(absent_parents)))
    module_files = set()
    for invocation, observation in observations.items():
        launch = launches.get(invocation)
        if not launch:
            continue
        if any(observation.get(k) != launch.get(k) for k in ("entrypoint", "parent", "pid")):
            errors.append(f"Launch/observer identity mismatch: {invocation}")
        if observation.get("status") != "COMPLETE":
            errors.append(f"Observer did not complete: {invocation}")
        modules = observation.get("modules")
        if not isinstance(modules, dict) or not modules:
            errors.append(f"Empty or invalid observed modules: {invocation}")
        else:
            module_files.update(modules.values())
        module_files.add(observation["entrypoint"])
    unexpected, scientific = [], []
    runtime_roots = [Path(p).resolve() for p in runtime_roots]
    for filename in sorted(module_files):
        file = Path(filename).resolve()
        if file.is_relative_to(package_root):
            scientific.append({"location": "PACKAGE", "path": str(file.relative_to(package_root))})
        elif file.is_relative_to(run_root):
            scientific.append({"location": "RUN", "path": str(file.relative_to(run_root))})
        elif not any(file.is_relative_to(root) for root in runtime_roots):
            unexpected.append(str(file))
    if not module_files:
        errors.append("No module filenames were observed")
    if unexpected:
        errors.append("Observed module files outside the package, run and selected Python environment")
    return {"status": "FAIL" if errors else "PASS", "errors": errors,
            "launched_process_count": len(launches), "observed_process_count": len(observations),
            "missing_invocations": missing, "unregistered_invocations": extra,
            "unwrapped_launch_records": unobserved,
            "scientific_module_files": scientific, "unexpected_module_files": unexpected,
            "observed_module_file_count": len(module_files),
            "expected_top_level_process_count": sum(expected_roots.values()),
            "expected_total_process_count": sum(expected_roots.values()) + sum(expected_children.values()),
            "coverage": "Explicit observer for each registered Python process, checked against declared parent/child counts; all other observed scripts must be leaves.",
            "external_dependencies": "Python standard library and installed packages only; raw paths and process identities are under import_audit/."}

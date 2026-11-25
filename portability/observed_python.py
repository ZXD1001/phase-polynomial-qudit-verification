"""Run one script under an explicit observer; internal to run_revision_checks."""
from __future__ import annotations
import atexit
import json
import os
from pathlib import Path
import runpy
import sys

from import_audit import install_launcher, write_json


def main():
    if not __debug__:
        raise RuntimeError("Run without -O/-OO")
    directory = Path(sys.argv[1]).resolve()
    invocation = sys.argv[2]
    target_args = sys.argv[3:]
    entrypoint = Path(target_args[0]).resolve()
    launch = json.loads((directory / (invocation + ".launch.json")).read_text())
    if launch["entrypoint"] != str(entrypoint) or launch["argv"] != target_args:
        raise RuntimeError("Observer invocation does not match its launch registration")
    def record():
        modules = {name: str(Path(module.__file__).resolve())
                   for name, module in list(sys.modules.items())
                   if getattr(module, "__file__", None)}
        write_json(directory / (invocation + ".observed.json"),
                   {"invocation": invocation, "parent": launch["parent"], "pid": os.getpid(),
                    "entrypoint": str(entrypoint), "status": "COMPLETE", "modules": modules})
    atexit.register(record)
    install_launcher(directory, invocation)
    sys.argv = target_args
    # Match the script-directory search entry used by direct Python execution.
    # This does not add the observer directory to scientific import resolution.
    sys.path[0] = str(entrypoint.parent)
    runpy.run_path(str(entrypoint), run_name="__main__")


if __name__ == "__main__":
    main()

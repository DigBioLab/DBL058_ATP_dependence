#!/usr/bin/env python3
"""
Thin wrapper around jupyter_utils.make_sub_script.

Two problems with the shared helper, neither worth patching in DBL_secrets
(other campaigns depend on its current behaviour):

  1. It takes no `logs_dir`, but every step in this pipeline passes one, so
     every call raises TypeError. The pipeline as committed cannot submit a
     single job against the current jupyter_utils.
  2. It hardcodes `#BSUB -o logs/<job>_%I.stdout` — a path relative to whatever
     directory bsub happens to run in — and creates `logs/` relative to the
     *calling process*. For steps that run as LSF jobs the CWD is not the
     campaign directory, so logs land somewhere arbitrary or not at all.

This wrapper accepts `logs_dir`, delegates, then rewrites the -o line to an
absolute path and adds a matching -e line.
"""
import os
import sys


def make_sub_script(cmds, n_task, logs_dir=None, **kwargs):
    """make_sub_script with a working, absolute logs_dir. Returns the .sh path."""
    sys.path.insert(0, os.environ.get(
        "ENZDES_SCRIPTS", "/work3/s212001/DBL_secrets/enzdes_pipeline-main"))
    import jupyter_utils

    jupyter_utils.make_sub_script(cmds=cmds, n_task=n_task, **kwargs)
    sub_script = cmds.replace(".cmds", ".sh")

    if logs_dir:
        os.makedirs(logs_dir, exist_ok=True)
        job_name = kwargs.get("job_name", "slurm")
        out = []
        for line in open(sub_script):
            if line.startswith("#BSUB -o "):
                out.append(f"#BSUB -o {logs_dir}/{job_name}_%I.stdout\n")
                out.append(f"#BSUB -e {logs_dir}/{job_name}_%I.stderr\n")
            else:
                out.append(line)
        with open(sub_script, "w") as f:
            f.writelines(out)

    return sub_script

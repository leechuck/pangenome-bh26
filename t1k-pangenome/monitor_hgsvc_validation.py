"""Automatically release the frozen HGSVC pipeline as prerequisites complete."""
import subprocess
import sys
import time
from pathlib import Path

HERE=Path(__file__).resolve().parent
OUT=HERE/'validation-next'


def run(script):
    subprocess.run([sys.executable,str(HERE/script)],check=True,timeout=300)


while True:
    if not (OUT/'FROZEN_VALIDATION.json').exists():
        run('freeze_hgsvc_validation.py')
    if (OUT/'FROZEN_VALIDATION.json').exists():
        if not (OUT/'EXECUTION_LAUNCH.json').exists():
            run('launch_hgsvc_validation.py')
        from advance_hgsvc_validation import main
        done,failed=main()
        if failed:raise RuntimeError('Terminal Slurm failure requires diagnosis; no automatic resubmission')
        if done:
            print('All HGSVC inference stages complete; independent provenance review and scoring still required',flush=True)
            break
    time.sleep(60)

#!/usr/bin/env python3
"""Read current Slurm states for the submitted mixed-Asian pilot."""
import json,subprocess
from pathlib import Path
from remote import remote
R=Path(__file__).resolve().parent
jobs=json.loads((R/'source/jobs.json').read_text())['jobs']
remote('sacct -j '+','.join(map(str,jobs.values()))+' --format=JobID%24,JobName%22,State%24,ExitCode,Elapsed')

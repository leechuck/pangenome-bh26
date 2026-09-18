#!/usr/bin/env python3
"""Stream existing public benchmark assets between clusters without relocating prefixes."""
import json
import shlex
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
config = json.loads((HERE/'ibex.json').read_text())
# Omit whole-genome CRAMs, results and stale graph/index products.
assets = ['mm/envs/asian50-spechla', 'mm/envs/typing', 'mm/envs/pggb053',
          'cactus/cactus-bin-v3.3.0', 'spechla-pg/source', 'hla-typer/source/genes',
          'hla-typer/t1kdb', 'hla-typer/reads']
source = 'tar -C /home/leechuck/hla -cf - ' + ' '.join(map(shlex.quote, assets))
hop = 'ssh -o BatchMode=yes -o IdentitiesOnly=yes -i ~/.ssh/id_ed25519 a001 ' + shlex.quote(source)
dest = config['workdir'] + '/hla'
subprocess.run(['ssh', '-o', 'BatchMode=yes', config['host'], 'mkdir -p ' + shlex.quote(dest)], check=True)
with (HERE/'transfer-source.log').open('w') as src_log, (HERE/'transfer-destination.log').open('w') as dst_log:
    src = subprocess.Popen(['ssh', '-o', 'BatchMode=yes', 'ddbj', hop], stdout=subprocess.PIPE, stderr=src_log)
    dst = subprocess.Popen(['ssh', '-o', 'BatchMode=yes', config['host'], 'tar -C '+shlex.quote(dest)+' -xf -'], stdin=src.stdout, stderr=dst_log)
    src.stdout.close()
    dst_rc = dst.wait()
    src_rc = src.wait()
if src_rc or dst_rc:
    raise SystemExit(f'Transfer failed: source={src_rc}, destination={dst_rc}; safe to rerun')
(HERE/'transfer-complete.json').write_text(json.dumps(dict(assets=assets, workdir=config['workdir'], source_exit=src_rc, destination_exit=dst_rc), indent=2)+'\n')

"""Record HG00611's reproduced upstream failure without scoring partial output."""
import json
import shlex
from pathlib import Path
from launch_ibex import remote

SCRIPT = r'''
from pathlib import Path
import hashlib,json,time
r=Path('/ibex/scratch/projects/c2014/rob/dogohla-benchmark/hla/dogohla-series20260918/gourraud')
attempts=[r/'attempts/arpack-HG00611-hprc/A',r/'arms/hprc/runs/HG00611/A']
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
records=[json.loads((p/'manifest.json').read_text()) for p in attempts]
assert all(m['status']=='failed' and m['configuration']['donor']=='HG00611' for m in records)
assert records[0]['configuration']==records[1]['configuration']
assert records[0]['finished'] < records[1]['started']
signature='ARPACK error -9: Starting vector is zero.'
assert all(signature in (p/'phase.HLA_A.log').read_text() for p in attempts)
audit=dict(donor='HG00611',arm='hprc',reason=signature,identical_configuration=True,
           attempts=[dict(path=str(p),manifest_sha256=sha(p/'manifest.json'),
                          log_sha256=sha(p/'phase.HLA_A.log')) for p in attempts],
           retry_job='52042583',scoring_policy='Zero credit for every eligible planned genotype; retain denominator')
for variant in ('DogoHLA','DogoHLA-no-graph'):
 p=r/'arms/hprc/final/runs/HG00611'/variant
 if p.exists():
  old=json.loads((p/'manifest.json').read_text())
  assert old['status']=='failed' and old['configuration'].get('upstream_failure')==audit
 else:
  p.mkdir(parents=True)
  (p/'manifest.json').write_text(json.dumps(dict(status='failed',finished=time.time(),
     configuration=dict(donor='HG00611',arm=variant,upstream_failure=audit),
     error='Blocked by reproducible upstream phasing failure'),indent=2)+'\n')
 assert not (p/'COMPLETE').exists()
 (p/'TERMINAL_FAILURE.json').write_text(json.dumps(dict(audit,manifest_sha256=sha(p/'manifest.json')),indent=2)+'\n')
p=attempts[-1]
(p/'TERMINAL_FAILURE.json').write_text(json.dumps(dict(audit,manifest_sha256=sha(p/'manifest.json')),indent=2)+'\n')
print(json.dumps(audit,indent=2))
'''

if __name__ == '__main__':
    result = remote('python3 -c ' + shlex.quote(SCRIPT))
    Path(__file__).with_name('arpack-HG00611-final.json').write_text(result+'\n')
    print(result)

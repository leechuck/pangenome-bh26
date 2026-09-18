"""Record a reproduced upstream ARPACK failure without scoring partial output."""
import argparse
import json
import re
import shlex
from pathlib import Path
from launch_ibex import remote

SCRIPT = r'''
from pathlib import Path
import hashlib,json,time
r=Path('/ibex/scratch/projects/c2014/rob/dogohla-benchmark/hla/dogohla-series20260918/gourraud')
attempts=[r/f'attempts/arpack-{donor}-hprc/A',r/f'arms/hprc/runs/{donor}/A']
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
records=[json.loads((p/'manifest.json').read_text()) for p in attempts]
assert all(m['status']=='failed' and m['configuration']['donor']==donor for m in records)
assert records[0]['configuration']==records[1]['configuration']
assert records[0]['finished'] < records[1]['started']
signature='ARPACK error -9: Starting vector is zero.'
log=f'phase.HLA_{gene}.log'
assert all(signature in (p/log).read_text() for p in attempts)
audit=dict(donor=donor,arm='hprc',reason=signature,identical_configuration=True,
           attempts=[dict(path=str(p),manifest_sha256=sha(p/'manifest.json'),
                          log_sha256=sha(p/log)) for p in attempts],
           retry_job=retry_job,scoring_policy='Zero credit for every eligible planned genotype; retain denominator')
for variant in ('DogoHLA','DogoHLA-no-graph'):
 p=r/f'arms/hprc/final/runs/{donor}'/variant
 if p.exists():
  old=json.loads((p/'manifest.json').read_text())
  assert old['status']=='failed' and old['configuration'].get('upstream_failure')==audit
 else:
  p.mkdir(parents=True)
  (p/'manifest.json').write_text(json.dumps(dict(status='failed',finished=time.time(),
     configuration=dict(donor=donor,arm=variant,upstream_failure=audit),
     error='Blocked by reproducible upstream phasing failure'),indent=2)+'\n')
 assert not (p/'COMPLETE').exists()
 (p/'TERMINAL_FAILURE.json').write_text(json.dumps(dict(audit,manifest_sha256=sha(p/'manifest.json')),indent=2)+'\n')
p=attempts[-1]
(p/'TERMINAL_FAILURE.json').write_text(json.dumps(dict(audit,manifest_sha256=sha(p/'manifest.json')),indent=2)+'\n')
print(json.dumps(audit,indent=2))
'''

if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--donor', required=True)
    p.add_argument('--gene', choices=['A','B','C','DPA1','DPB1','DQA1','DQB1','DRB1'], required=True)
    p.add_argument('--retry-job', required=True)
    a = p.parse_args()
    if not re.fullmatch(r'(HG|NA)[0-9]+', a.donor) or not a.retry_job.isdigit():
        p.error('Invalid donor or job identifier')
    setup = f'donor={a.donor!r};gene={a.gene!r};retry_job={a.retry_job!r}\n'
    result = remote('python3 -c ' + shlex.quote(setup+SCRIPT))
    Path(__file__).with_name(f'arpack-{a.donor}-final.json').write_text(result+'\n')
    print(result)

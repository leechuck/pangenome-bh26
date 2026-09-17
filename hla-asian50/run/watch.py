#!/usr/bin/env python3
"""Monitor live submitted jobs every 50 seconds; stop on failure or completion."""
import datetime,json,time
from pathlib import Path
from remote import remote,ROOT
R=Path(__file__).resolve().parent
while True:
 try:
  r=remote('python3 '+ROOT+'/monitor_remote.py',capture_output=True,timeout=45)
  s=json.loads(r.stdout);(R/'source/latest_status.json').write_text(json.dumps(s,indent=2)+'\n')
  print(s['utc'], 'indices',s['index_logs_complete'], '/10', ' | '.join(t+': '+str(s['jobs'][t]['states'])+' files='+str(s['outputs'][t]['donors_with_files']) for t in s['jobs']),flush=True)
  if any(v['failed'] for v in s['jobs'].values()) or s.get('panel_index',{}).get('state') in ['FAILED','TIMEOUT','OUT_OF_MEMORY','CANCELLED']:
   print('FAILURE_REQUIRES_INSPECTION',flush=True);raise SystemExit(2)
  if s['all_complete']:print('ALL_SUBMITTED_RUNS_COMPLETE',flush=True);break
 except (TimeoutError,ValueError) as e:print('MONITOR_RETRY',repr(e),flush=True)
 time.sleep(50)

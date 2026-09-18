from pathlib import Path
import json,subprocess
root=Path('/home/leechuck/hla/t1k-pangenome')
vg='/home/leechuck/hla/cactus/cactus-bin-v3.3.0/bin/vg'
for panel,gene in [('hprc','DRB1'),('hprc_asian','DQB1')]:
 gfa=root/'graphs/build-v1'/panel/gene/'graph.gfa'
 lengths={int(f[1]):len(f[2]) for line in gfa.open() if (f:=line.rstrip().split('\t'))[0]=='S'}
 print('GRAPH',panel,gene,'nodes',len(lengths),'max_id',max(lengths),'max_length',max(lengths.values()),flush=True)
 source=root/'development/mapping-v2/HG00658'/panel/'unsampled'/gene/'alignments.jsonl'
 bad=0
 for line in source.open():
  a=json.loads(line)
  for m in a.get('path',{}).get('mapping',[]):
   p=m.get('position',{});node=int(p.get('node_id',0))
   if node not in lengths:
    print('MISSING',json.dumps(dict(node=node,position=p,mapping_keys=list(m),alignment_keys=list(a),score=a.get('score'),mapping_quality=a.get('mapping_quality'))),flush=True)
    bad+=1;break
  if bad==3:break
subprocess.run([vg,'convert','--help'])
subprocess.run([vg,'giraffe','--help'])

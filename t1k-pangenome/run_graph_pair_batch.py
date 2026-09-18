"""Run unchanged pair refinement across the predefined development panel."""
import argparse
import json
from pathlib import Path
from evidence_io import sha
from run_graph_pair_refinement import run


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--index',type=int,required=True)
    p.add_argument('--config',type=Path,required=True)
    p.add_argument('--internal-pairs',action='store_true')
    p.add_argument('--include-pilot',action='store_true')
    p.add_argument('--pair-specific',action='store_true')
    a=p.parse_args();plan=json.loads(a.config.read_text())
    if sha(Path('/home/leechuck/hla/dogohla-series20260918/cohort.tsv'))!=plan['cohort_sha256']:
        raise ValueError('Development cohort changed')
    donors=([plan['existing_pilot']] if a.include_pilot else [])+plan['donors']
    if not 0<=a.index<4*len(donors):raise ValueError('Invalid development task')
    donor=donors[a.index//4]['donor']
    panel,selection=(('hprc','sampled'),('hprc','unsampled'),('hprc_asian','sampled'),('hprc_asian','unsampled'))[a.index%4]
    root=Path('/home/leechuck/hla/t1k-pangenome')
    if a.pair_specific and not a.internal_pairs:raise ValueError('Pair-specific batch requires internal pairs')
    version='v3-pairwise' if a.pair_specific else 'v2-internal' if a.internal_pairs else 'v1'
    run(root/'development/linear-v1/ipd_genome'/donor,root/'development/join-v3'/donor/panel/selection,
        root/'references/observed-v2'/panel,root/'development/library-calibration-v1'/donor/'COMPLETE.json',
        root/'development'/('graph-pair-refinement-'+version)/donor/panel/selection,
        internal_pairs=a.internal_pairs,pair_specific=a.pair_specific)

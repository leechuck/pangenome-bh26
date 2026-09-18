"""Verify the complete HGSVC execution freeze before dispatching any stage."""
import argparse
import json
from pathlib import Path
from evidence_io import sha

ROOT = Path('/home/leechuck/hla/t1k-pangenome')


def verify(assets, digest):
    frozen = assets/'FROZEN_VALIDATION.json'
    if sha(frozen)!=digest:raise ValueError('Execution freeze changed')
    plan = json.loads(frozen.read_text())
    for group in ('code_sha256','metadata_sha256'):
        for name,expected in plan[group].items():
            if sha(assets/name)!=expected:raise ValueError('Frozen input changed: '+name)
    for name,expected in plan['asset_sha256'].items():
        if sha(ROOT/name)!=expected:raise ValueError('Reference asset changed: '+name)
    return plan


def main(stage,index,assets,digest):
    verify(assets,digest)
    output = ROOT/'hgsvc-validation-v1/validation-runs'
    if stage=='baseline':
        from run_hgsvc_baseline import run
        run(index,assets,output/'frozen-t1k-v1')
    elif stage=='genomic':
        from run_hgsvc_genomic import run
        run(index,assets,output/'genomic-ipd-v1')
    elif stage=='anchor':
        from run_hgsvc_anchor import run
        run(index,assets)
    else:
        from run_hgsvc_graph import run
        run(stage,index,assets,digest)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--stage',required=True);p.add_argument('--index',type=int,required=True)
    p.add_argument('--assets',type=Path,required=True);p.add_argument('--freeze-sha256',required=True)
    a=p.parse_args();main(a.stage,a.index,a.assets,a.freeze_sha256)

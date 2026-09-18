"""Fetch development refinement calls and decisions, never reserved outcomes."""
import subprocess
from score_graph_pair_refinement import HERE,run


if __name__=='__main__':
    snapshot=HERE/'work/graph-pair-refinement-snapshot';snapshot.mkdir(parents=True,exist_ok=True)
    source=('hohndor@ilogin.ibex.kaust.edu.sa:'
            '/ibex/scratch/projects/c2014/rob/dogohla-benchmark/hla/t1k-pangenome/development/graph-pair-refinement-v1/')
    subprocess.run(['rsync','-az','--include=*/','--include=manifest.json','--include=COMPLETE.json',
                    '--include=calls.json','--include=decisions.json','--exclude=*',source,str(snapshot)+'/'],check=True)
    run()

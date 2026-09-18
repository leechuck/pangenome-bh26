"""Fetch small development-control outputs and score with hash verification."""
from pathlib import Path
import subprocess
from score_linear_controls import main


if __name__ == '__main__':
    here = Path(__file__).resolve().parent
    snapshot = here/'work/linear-snapshot'
    snapshot.mkdir(parents=True,exist_ok=True)
    source = ('hohndor@ilogin.ibex.kaust.edu.sa:'
              '/ibex/scratch/projects/c2014/rob/dogohla-benchmark/hla/'
              't1k-pangenome/development/linear-v1/')
    subprocess.run(['rsync','-az','--include=*/','--include=manifest.json',
                    '--include=COMPLETE.json','--include=calls.json',
                    '--include=*_genotype.tsv','--exclude=*',source,str(snapshot)+'/'],check=True)
    # Completed manifests must match completion records and every output hash.
    # A concurrently changing/partial transfer fails rather than being scored.
    main(snapshot,here/'development/linear-analysis')

#!/usr/bin/env python3
"""Read-only input availability audit; no reads, genotypes or labels inferred."""
import csv,json,os,subprocess
from pathlib import Path
root=Path(__file__).resolve().parent
rows=list(csv.DictReader(open(root/'source/donors.tsv'),delimiter='\t'))
records=[]
for r in rows:
    d=r['donor']
    full=Path('/usr/local/shared_data/public-human-genomes/GRCh38/1000Genomes/CRAM')/d/(d+'.cram')
    mhc=Path('/home/asianhla/data/upload/1000G_MHC/cram')/(d+'.mhc.cram')
    fresh=Path('/home/leechuck/hla/codex-targeted/fresh_reads')/(d+'.mhc.cram')
    t1k=Path('/home/asianhla/data/upload/1000G_MHC/t1k')/d/(d+'_genotype.tsv')
    records.append(dict(donor=d,full_cram=str(full),full_readable=int(os.access(full,os.R_OK)),mhc_cram=str(mhc if mhc.is_file() else fresh),mhc_readable=int(mhc.is_file() or fresh.is_file()),existing_t1k=int(t1k.is_file()),assembly1=int(Path('/home/leechuck/hla/mhc_all',d+'_1.mhc.fa').is_file()),assembly2=int(Path('/home/leechuck/hla/mhc_all',d+'_2.mhc.fa').is_file())))
with open(root/'source/remote_availability.tsv','w') as f:
    w=csv.DictWriter(f,fieldnames=list(records[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(records)
print(json.dumps({k:sum(r[k] for r in records) for k in ['full_readable','mhc_readable','existing_t1k','assembly1','assembly2']},indent=2))

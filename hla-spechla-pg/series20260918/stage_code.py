#!/usr/bin/env python3
"""Stage an isolated IPD-365 experiment; never edit frozen v0.1 code or references."""
import csv,hashlib,io,json,subprocess,tarfile
from pathlib import Path
HERE=Path(__file__).resolve().parent
REPO=HERE.parents[1]
OUT=HERE/'transfer-cache/code-stage'
ROOT='/home/leechuck/hla/dogohla-series20260918'
FILES=['dogohla.py','experiment.py','phase_linkage.py','build_refs.py','graph_diagnostic.py','structural_overlay.py','spechla_pg.sh','build_graph_pggb.sh','designate.py','native_long_control.py']

def main():
 (OUT/'code').mkdir(parents=True,exist_ok=True);(OUT/'source').mkdir(exist_ok=True)
 manifest={}
 for name in FILES:
  original=(HERE.parent/name).read_text();text=original
  # These changes affect only the separately staged experiment.
  if name=='experiment.py':
   text=text.replace("ENV = Path('/home/leechuck/hla/mm/envs/asian50-spechla')",f"ENV = Path('{ROOT}/env')")
  if name=='build_refs.py':
   text=text.replace("SPECHLA_DB=Path('/home/leechuck/hla/mm/envs/asian50-spechla/share/spechla/db')",f"SPECHLA_DB=Path('{ROOT}/env/share/spechla/db')")
   text=text.replace("PANEL=Path('/home/leechuck/hla/hla-typer/source/genes')","PANEL=Path(os.environ['DOGOHLA_PANEL'])")
  if name=='spechla_pg.sh':text=text.replace('ENV=/home/leechuck/hla/mm/envs/asian50-spechla',f'ENV={ROOT}/env').replace('-p Unknown','-p nonuse')
  if name in ('dogohla.py','experiment.py','native_long_control.py','structural_overlay.py'):
   text=text.replace("'Unknown'","'nonuse'")
  if name=='dogohla.py':text=text.replace('DogoHLA 0.1.0','DogoHLA 0.2.0-experimental')
  if name.endswith('.py'):compile(text,name,'exec')
  (OUT/'code'/name).write_text(text)
  manifest[name]={'source_sha256':hashlib.sha256(original.encode()).hexdigest(),'staged_sha256':hashlib.sha256(text.encode()).hexdigest()}
 # Database preparation imports this staged code from its sibling directory.
 (OUT/'series').mkdir(exist_ok=True)
 for name in ('prepare_series.py','configure_series.py','run_t1k4.py','prepare_arm.sh','prepare_ibex.sh'):(OUT/'series'/name).write_bytes((HERE/name).read_bytes())
 for name in ('excluded_paths.tsv','validation_groups.tsv','donors_folds.tsv'):(OUT/'source'/name).write_bytes((HERE/name).read_bytes())
 (OUT/'source/haplotype_donors.tsv').write_bytes((REPO/'hla-typer/source/haplotype_donors.tsv').read_bytes())
 rows=list(csv.DictReader((HERE/'cohort.tsv').open(),delimiter='\t'))
 urls=dict(line.split('\t')[:2] for line in (REPO/'hla-targeted/source/C4Investigator/resources/1000Genomes_resources/tgp_full_30x.tsv').read_text().splitlines())
 for r in rows:r['cram']='https://s3.amazonaws.com/1000genomes/'+urls[r['donor']]
 for name,keys in [('cohort.tsv',list(rows[0])),('recruit.tsv',['donor','cram'])]:
  with (OUT/name).open('w') as f:
   w=csv.DictWriter(f,fieldnames=keys,delimiter='\t',lineterminator='\n',extrasaction='ignore');w.writeheader();w.writerows(rows)
 # Same recruitment regions/read-pair policy as the existing benchmark.
 text=(REPO/'hla-typer/stageA/recruit.sbatch').read_text()
 text='\n'.join(l for l in text.splitlines() if not l.startswith('#SBATCH'))+'\n'
 text=text.replace('asian50-genotyping','asian50-spechla').replace('samtools view -@4','samtools view -@3').replace('samtools collate -@4','samtools collate -@1').replace('samtools fastq -@4','samtools fastq -@1').replace('samtools flagstat -@4','samtools flagstat -@3')
 (OUT/'code/recruit.sh').write_text(text)
 manifest['changes']=['IPD 3.65 candidate and naming database','Same nonuse frequency policy for all updated methods','Distinct panel references controlled by DOGOHLA_PANEL','Recruitment policy unchanged, thread requests bounded']
 (OUT/'code/MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n')
 archive=io.BytesIO()
 with tarfile.open(fileobj=archive,mode='w:gz') as tar:
  for p in OUT.rglob('*'):
   if p.is_file():tar.add(p,arcname=str(p.relative_to(OUT)))
 (HERE/'transfer-cache/code.tar.gz').write_bytes(archive.getvalue())
 print('Staged',len(archive.getvalue()),'bytes')

if __name__=='__main__':main()

#!/usr/bin/env python3
"""Build a separate IPD 3.65 SpecHLA environment and size-controlled panel inputs."""
import argparse,collections,csv,hashlib,json,os,random,shutil,subprocess,sys,io
from pathlib import Path
from Bio import SeqIO
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'code'))
from build_refs import GENES,BIN_GENES,read_fasta,write_fasta
from experiment import sha,json_write
from native_long_control import patch_empty_assembly
REAL=Path('/home/leechuck/hla/mm/envs/asian50-spechla')
OLD=Path('/home/leechuck/hla/spechla-pg')
PANEL=Path('/home/leechuck/hla/hla-typer/source/genes')
V2='hla_gen.format.filter.extend.DRB.no26789.v2.fasta'

def run(cmd):
 print('COMMAND',list(map(str,cmd)),flush=True);subprocess.run(list(map(str,cmd)),check=True)

def build(root,threads):
 if (root/'PREPARED.json').exists():raise FileExistsError('Series already prepared')
 env=root/'env';env.mkdir();(env/'bin').symlink_to(REAL/'bin');(env/'lib').symlink_to(REAL/'lib')
 installed=env/'share/spechla';installed.mkdir(parents=True)
 shutil.copytree(REAL/'share/spechla/script',installed/'script')
 assembly=installed/'script/run.assembly.realign.sh'
 assembly.write_text(patch_empty_assembly(assembly.read_text()))
 db=installed/'db';shutil.copytree(REAL/'share/spechla/db',db,symlinks=False)
 os.environ.update(PATH=str(REAL/'bin')+':'+os.environ['PATH'],CONDA_PREFIX=str(REAL),LD_LIBRARY_PATH=str(REAL/'lib'))
 data=Path('/home/leechuck/hla/hla-typer/t1kdb/hla.dat')
 genomic={g:[] for g in BIN_GENES};cds={g:[] for g in GENES};exons=[];empty_records=[]
 # The same pinned flat file used to build T1K's IPD 3.65 database.
 # IPD uses a six-field legacy ID header. Adapt only that header to the
 # seven-field EMBL layout expected by Biopython; sequences/features are untouched.
 normalized=root/'ipd365.biopython.embl'
 with data.open() as src, normalized.open('w') as dst:
  for line in src:
   if line.startswith('ID   '):
    fields=line.rstrip().split(';')
    if len(fields)==6 and fields[2].strip()=='standard':
     line=';'.join(fields[:2]+[' linear',fields[3],' STD']+fields[4:])+'\n'
   dst.write(line)
 for record in SeqIO.parse(normalized,'embl'):
  annotated=[a for f in record.features for a in f.qualifiers.get('allele',[])]
  allele=(annotated[0] if annotated else record.description.split(',')[0]).removeprefix('HLA-');gene=allele.split('*')[0]
  if gene not in BIN_GENES:continue
  sequence=str(record.seq).upper()
  if not sequence:
   empty_records.append(allele);continue
  if set(sequence)-set('ACGTN'):raise ValueError('Unexpected nucleotide alphabet')
  genomic[gene].append((allele,sequence))
  if gene not in GENES:continue
  coding=[f for f in record.features if f.type=='CDS']
  if coding:cds[gene].append((allele,str(coding[0].extract(record.seq)).upper()))
  wanted={'2','3'} if gene in ('A','B','C') else {'2'}
  for f in record.features:
   number=(f.qualifiers.get('number') or [''])[0]
   if f.type=='exon' and number in wanted:
    exons.append((f'HLA-{allele}|Exon|{number}',str(f.extract(record.seq)).upper()))
 assert all(genomic[g] for g in BIN_GENES) and all(cds[g] for g in GENES) and exons
 # Replace candidate databases; retain coordinate references and their existing indexes.
 for directory in (db/'HLA/whole',db/'HLA/exon'):
  for p in directory.iterdir():
   if p.is_file():p.unlink()
 for gene in GENES:
  # Official genomic FASTAs are the same pinned phase reference used in 0.1.0.
  seqs=[(n.split()[1],s.upper()) for n,s in read_fasta(OLD/'source/imgt_gen'/f'{gene}_gen.fasta')]
  write_fasta(db/'HLA/whole'/f'HLA_{gene}.fasta',seqs)
  write_fasta(db/'HLA/exon'/f'HLA_{gene}.fasta',cds[gene])
  for folder,prefix in [('whole',f'HLA_{gene}'),('exon',f'HLA_{gene}.fasta')]:
   fa=db/'HLA'/folder/(f'HLA_{gene}.fasta')
   run(['makeblastdb','-in',fa,'-dbtype','nucl','-out',db/'HLA'/folder/prefix]);run(['samtools','faidx',fa])
 write_fasta(db/'HLA/whole/HLA_DRB1.exon.fasta',cds['DRB1'])
 run(['makeblastdb','-in',db/'HLA/whole/HLA_DRB1.exon.fasta','-dbtype','nucl','-out',db/'HLA/whole/HLA_DRB1.exon'])
 for p in (db/'HLA').glob('hla_exons.fasta*'):p.unlink()
 write_fasta(db/'HLA/hla_exons.fasta',exons)
 run(['makeblastdb','-in',db/'HLA/hla_exons.fasta','-dbtype','nucl','-out',db/'HLA/hla_exons.fasta'])
 for name in ('hla_nom_g.txt','Allelelist.txt'):shutil.copy2(root/'source'/name,db/'HLA'/name)
 assert '3.65' in (db/'HLA/hla_nom_g.txt').read_text()[:1000]
 for p in (db/'ref').glob(V2+'*'):p.unlink()
 write_fasta(db/'ref'/V2,[r for g in BIN_GENES for r in genomic[g]])
 run(['bowtie2-build','--threads',threads,'-q',db/'ref'/V2,db/'ref'/V2])
 excluded={r['hap_id'] for r in csv.DictReader(open(root/'source/excluded_paths.tsv'),delimiter='\t')}
 meta=list(csv.DictReader(open(root/'source/haplotype_donors.tsv'),delimiter='\t'))
 units=collections.defaultdict(list)
 for r in meta:
  if r['hap_id'] not in excluded and r['cohort']!='REF':units[(r['donor_id'],'HPRC' if r['cohort'].startswith('HPRC') else r['cohort'])].append(r['hap_id'])
 hprc=sorted(u for u in units if u[1]=='HPRC');other=sorted(u for u in units if u[1]!='HPRC')
 target=sum(len(units[u]) for u in hprc)
 chosen=list(other);n=sum(len(units[u]) for u in chosen);pool=hprc[:];random.Random(20260918).shuffle(pool)
 if n>target:raise ValueError('Asian-enriched set exceeds matched target')
 for u in pool:
  if n==target:break
  if n+len(units[u])<=target:chosen.append(u);n+=len(units[u])
 assert n==target
 ref={r['hap_id'] for r in meta if r['cohort']=='REF'}
 arms={'full':list(units),'hprc':hprc,'asian_matched':chosen};audit={}
 for arm,us in arms.items():
  keep=ref|{h for u in us for h in units[u]};assert not keep&excluded
  d=root/'panels'/arm;d.mkdir(parents=True);audit[arm]={'haplotypes':len(keep),'genes':{}}
  (d/'members.txt').write_text('\n'.join(sorted(keep))+'\n')
  for gene in BIN_GENES:
   records=[(name,seq) for name,seq in read_fasta(PANEL/f'HLA-{gene}.fa') if '#'.join(name.split('#')[:2]) in keep]
   write_fasta(d/f'HLA-{gene}.fa',records);audit[arm]['genes'][gene]=len(records)
 json_write(root/'PREPARED.json',dict(ipd='3.65.0',ignored_empty_records=empty_records,flatfile_sha256=sha(data),exclusions_sha256=sha(root/'source/excluded_paths.tsv'),panels=audit,
   naming_policy='nonuse: no legacy frequency filter; same policy across updated SpecHLA and DogoHLA',
   reference_coordinates='Original SpecHLA linear references retained; allele genomic/CDS/ARS and binning databases updated',
   database_files={str(p.relative_to(db)):sha(p) for p in db.rglob('*') if p.is_file() and p.suffix in ('.fasta','.fa','.txt')}))
 print('SERIES_PREPARED',flush=True)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('root',type=Path);p.add_argument('--threads',type=int,default=4);a=p.parse_args();build(a.root,a.threads)

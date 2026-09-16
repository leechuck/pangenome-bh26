#!/usr/bin/env python3
"""Build development-only C4 gene alignments and module-boundary probe contexts."""
import csv,gzip,hashlib,json,re
from pathlib import Path
from Bio import SeqIO
ROOT=Path(__file__).resolve().parent
def read(p):return list(csv.DictReader(open(p),delimiter='\t'))
def write(p,rows):
    with open(p,'w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
def main():
    heldout={r['family'] for r in read(ROOT/'source/validation_donors.tsv')}
    panel=[r for r in read(ROOT.parent/'hla-structural/source/catalogue.tsv') if r['locus']=='RCCX' and r['family'] not in heldout]
    with gzip.open(ROOT.parent/'hla-structural/source/locus_sequences.fa.gz','rt') as f:seq={r.id:str(r.seq) for r in SeqIO.parse(f,'fasta')}
    genes=[]
    with open(ROOT/'source/training_C4_genes.fa','w') as out:
        for r in panel:
            hap=r['hap_id'];entry=hap.replace('#','_');s=seq[r['sequence_id']];offset=int(r['start0']);n=0
            with gzip.open(ROOT.parent/f'hla-analysis/source/gtf_all/{entry}.gtf.gz','rt') as f:
                for line in f:
                    if line.startswith('#'):continue
                    x=line.rstrip().split('\t')
                    if len(x)<9 or x[2]!='gene':continue
                    m=re.search(r'gene_name "(C4[AB][LS])"',x[8])
                    if not m:continue
                    assert x[0]==r['contig'];a=int(x[3])-1;b=int(x[4]);g=s[a-offset:b-offset]
                    assert len(g)==b-a and set(g)<=set('ACGT')
                    if x[6]=='-':g=g.translate(str.maketrans('ACGT','TGCA'))[::-1]
                    n+=1;name=entry+'__C4_'+str(n)
                    genes.append(dict(gene_id=name,hap_id=hap,donor=r['donor'],family=r['family'],annotated_form=m[1],start0=a,end0=b,strand=x[6]))
                    out.write('>'+name+'\n'+g+'\n')
    write(ROOT/'source/training_C4_genes.tsv',genes)
    reference=list(SeqIO.parse(ROOT/'source/C4Investigator/resources/c4only_onelines_oneDel_bShort.fasta','fasta'));a=str(reference[0].seq);b=str(reference[1].seq)
    blocks=list(re.finditer('[ACGT]+',a));assert len(blocks)==41
    sites=[blocks[25].start()+i-1 for i in [129,132,140,143,145]]
    assert ''.join(a[i] for i in sites)=='CGTGC' and ''.join(b[i] for i in sites)=='TCACT'
    gap=max(re.finditer(r'\.+',b),key=lambda x:len(x[0]));assert len(gap[0])==6367
    (ROOT/'source/C4_reference.fa').write_text('>C4A_long_reference\n'+a.upper()+'\n')
    meta=dict(reference_source='https://github.com/Hollenbach-lab/C4Investigator',diagnostic_sites0=sites,A_motif='CGTGC',B_motif='TCACT',HERV_start0=gap.start(),HERV_end0=gap.end(),training_haplotypes=len(panel),training_genes=len(genes),heldout_families=sorted(heldout))
    (ROOT/'source/target_definitions.json').write_text(json.dumps(meta,indent=2)+'\n')
    # Context around annotated tandem-module starts, not asserted SV breakpoints.
    byhap={r['hap_id']:r for r in panel};boundaries=[]
    with open(ROOT/'source/module_boundary_contexts.fa','w') as out:
        for x in read(ROOT.parent/'hla-analysis/results/rccx_homology_loci.tsv'):
            if x['hap_id'] not in byhap or x['family']!='WHR1' or x['eligible']!='1':continue
            r=byhap[x['hap_id']];s=seq[r['sequence_id']];pos=int(x['start0'])-int(r['start0'])
            if pos<256 or pos+256>len(s):continue
            name='BOUNDARY_'+str(len(boundaries));context=s[pos-256:pos+256]
            boundaries.append(dict(id=name,hap_id=r['hap_id'],family=r['family'],prototype=x['prototype'],start0=int(x['start0']),interpretation='annotated_module_start_context_not_proven_breakpoint'))
            out.write('>'+name+'\n'+context+'\n')
    write(ROOT/'source/module_boundary_contexts.tsv',boundaries)
    print(json.dumps(meta,indent=2));print('Boundary contexts',len(boundaries))
if __name__=='__main__':main()

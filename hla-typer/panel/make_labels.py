#!/usr/bin/env python3
"""Per-haplotype, per-gene allele labels and donor identity for panel construction (local, metadata only).

Output source/haplotype_labels.tsv.gz: one row per catalogue gene copy with the IPD-IMGT/HLA 3.65.0
exact-match lists from hla-analysis/results/sequence_catalogue.tsv, plus the manifest donor_id and
cohort, and the Locityper haplotype name (sample dots -> '_', '#' -> '.').
"""
import csv,gzip
from pathlib import Path
P=Path(__file__).resolve().parent.parent.parent
csv.field_size_limit(10**9)

def locityper_name(hap_id):
    sample,hap=hap_id.split('#');return sample.replace('.','_')+'.'+hap

def main():
    man={r['hap_id']:r for r in csv.DictReader(open(P/'hla-audit/2026-09-16/panel_manifest.tsv'),delimiter='\t')}
    out=gzip.open(Path(__file__).resolve().parent.parent/'source/haplotype_labels.tsv.gz','wt')
    cols=['hap_id','haplotype','donor_id','cohort','gene','cds_complete','assessment','qc_flags','exact_genomic_alleles','exact_cds_alleles','exact_protein_alleles','gene_sha256','cds_sha256']
    w=csv.DictWriter(out,fieldnames=cols,delimiter='\t',lineterminator='\n');w.writeheader()
    for r in csv.DictReader(open(P/'hla-analysis/results/sequence_catalogue.tsv'),delimiter='\t'):
        m=man[r['hap_id']]
        w.writerow(dict(hap_id=r['hap_id'],haplotype=locityper_name(r['hap_id']),donor_id=m['donor_id'],cohort=m['cohort'],gene=r['gene'],
                        **{k:r[k] for k in cols[5:]}))
    out.close()
    with open(Path(__file__).resolve().parent.parent/'source/haplotype_donors.tsv','w') as f:
        f.write('hap_id\thaplotype\tdonor_id\tcohort\n')
        for h,m in sorted(man.items()):f.write(f"{h}\t{locityper_name(h)}\t{m['donor_id']}\t{m['cohort']}\n")
if __name__=='__main__':main()

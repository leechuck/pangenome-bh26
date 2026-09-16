#!/usr/bin/env python3
"""Descriptive phase counterexamples from frozen validation predictions."""
import csv,json,re
from collections import Counter
from pathlib import Path
ROOT=Path(__file__).resolve().parent

def read(p):return list(csv.DictReader(open(p),delimiter='\t'))
def main():
    catalogue={r['structure_id']:r['structural_signature'] for r in read(ROOT.parent/'hla-structural/source/catalogue.tsv') if r['locus']=='RCCX'}
    cases=read(ROOT/'results/validation_case_analysis.tsv');rows=[]
    def forms(pair):return [re.findall(r'C4([AB][LS])[+-]',catalogue[s]) for s in pair.split(';')]
    def marginals(fs):
        flat=sum(fs,[])
        return dict(A=sum(f[0]=='A' for f in flat),B=sum(f[0]=='B' for f in flat),L=sum(f[1]=='L' for f in flat),S=sum(f[1]=='S' for f in flat))
    for r in cases:
        if r['targeted_correct']=='1':continue
        t=forms(r['truth_pair']);p=forms(r['predicted_pair'])
        rows.append(dict(donor=r['donor'],truth_C4_paths=' / '.join('-'.join(x) for x in t),imputed_C4_paths=' / '.join('-'.join(x) for x in p),truth_marginals=json.dumps(marginals(t),sort_keys=True),imputed_marginals=json.dumps(marginals(p),sort_keys=True),same_marginals=int(marginals(t)==marginals(p)),category=r['error_category']))
    with open(ROOT/'results/validation_structural_error_details.tsv','w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
    lines=['# What the structural errors mean','', 'The frozen targeted model recovers 21/24 full RCCX signature pairs, versus 22/24 with the old sketch. These three mismatches also change the inferred C4 gene-form arrangement; they are not solely CYP21/TNX prototype-label differences.','', 'AL = long C4A; AS = short C4A; BL = long C4B; BS = short C4B. A slash separates the two assembly haplotypes; chromosome ordering within a diploid pair is arbitrary.','', '| Donor | Assembly C4 paths | Imputed C4 paths | Same A/B and L/S marginals? |','|---|---|---|---|']
    for r in rows:lines.append(f"| {r['donor']} | {r['truth_C4_paths']} | {r['imputed_C4_paths']} | {'Yes' if r['same_marginals'] else 'No'} |")
    lines += ['', '**HG00673:** both configurations have A=3, B=1, L=3 and S=1. Correct marginal dosage cannot distinguish them. The true pair remains in the candidate set, but the sequence-profile ranking chooses the wrong combination.','', '**HG03540:** the true trimodular reference signature is absent after excluding the held-out families. The replacement has the same A/B and L/S marginal counts but different gene forms and order. The caller does not reliably recognize an unknown structure.','', '**HG04204:** the targeted A/B estimate rounds to A=3, B=1, while the assembly has A=2, B=2. This erroneous constraint excludes the correct pair, which the old sketch recovered. This is the single loss versus the old model; there are no gains in these 24 donors.','', 'Every validation donor has multiple dosage-compatible signature pairs (3–87). Thus none obtains a uniquely resolved full structure from these dosage constraints alone, even though ranking often predicts the assembly label correctly.','', 'A useful extended typing record would keep classical HLA types, C4 total/A/B/L/S dosage, the full set of compatible RCCX structures, and the source of phase evidence as separate fields. Direct long-molecule evidence linking diagnostic sites and module boundaries is the appropriate next test of arrangement; a unique ranking optimum is insufficient.','', 'No caller thresholds or reference paths were changed after observing these errors. These are descriptive analyses of the frozen validation run.','']
    (ROOT/'STRUCTURAL_ERRORS.md').write_text('\n'.join(lines));print('Wrote phase counterexamples for',len(rows),'errors')
if __name__=='__main__':main()

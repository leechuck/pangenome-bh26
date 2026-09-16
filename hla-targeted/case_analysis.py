#!/usr/bin/env python3
"""Post-evaluation error attribution; not used to modify the frozen caller."""
import csv,json,re
from collections import Counter,defaultdict
from pathlib import Path
ROOT=Path(__file__).resolve().parent

def read(p):return list(csv.DictReader(open(p),delimiter='\t'))
def write(p,rows):
    with open(p,'w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
def main():
    panel=defaultdict(list)
    for r in read(ROOT.parent/'hla-structural/source/catalogue.tsv'):
        if r['locus']=='RCCX':panel[r['donor']].append(r)
    pred={r['donor']:r for r in read(ROOT/'results/validation_dosage_predictions.tsv')}
    paths={(r['donor'],r['method']):r for r in read(ROOT/'results/validation_path_predictions.tsv')}
    reference={r['structure_id'] for r in json.loads((ROOT/'source/training_rows.json').read_text())}
    audit=defaultdict(list)
    for r in read(ROOT/'results/validation_C4_diagnostic_audit.tsv'):audit[r['donor']].append(r)
    result=[]
    for donor,dose in sorted(pred.items()):
        rs=panel[donor];truth=';'.join(sorted(r['structure_id'] for r in rs));forms=[g for r in rs for g in re.findall(r'C4([AB][LS])[+-]',r['structural_signature'])]
        counts={'C4':len(forms),'A':sum(f[0]=='A' for f in forms),'B':sum(f[0]=='B' for f in forms),'L':sum(f[1]=='L' for f in forms),'S':sum(f[1]=='S' for f in forms)}
        p=paths[donor,'targeted_dosage_paths'];old=paths[donor,'old_sketch_CN'];abl=paths[donor,'targeted_without_module_context']
        compatible=json.loads(p['dosage_compatible_pairs']);represented=all(r['structure_id'] in reference for r in rs)
        wrong=[k for k,v in counts.items() if dose[k+'_call']!='' and int(dose[k+'_call'])!=v];absent=[k for k in counts if dose[k+'_call']=='']
        correct=p['ranked_pair']==truth;oldcorrect=old['ranked_pair']==truth
        if correct:reason='correct_ranked_imputation'
        elif not represented:reason='truth_structure_absent_from_reference'
        elif wrong:reason='incorrect_marginal_dosage'
        elif truth.split(';') in compatible:reason='ranking_error_with_truth_compatible'
        elif absent:reason='insufficient_dosage_evidence'
        else:reason='truth_excluded_other_reason'
        result.append(dict(donor=donor,truth_pair=truth,predicted_pair=p['ranked_pair'],targeted_correct=int(correct),old_correct=int(oldcorrect),without_context_correct=int(abl['ranked_pair']==truth),change='gain' if correct and not oldcorrect else 'loss' if oldcorrect and not correct else 'unchanged',truth_represented=int(represented),truth_compatible=int(truth.split(';') in compatible),compatible_pairs=len(compatible),wrong_marginals=';'.join(wrong),uncalled_marginals=';'.join(absent),error_category=reason,assembly_annotation_disagreements=sum(r['status']!='match' for r in audit[donor]),noncanonical_fragments=dose['NONCANONICAL_fragments']))
    write(ROOT/'results/validation_case_analysis.tsv',result)
    summary=dict(n=len(result),changes=dict(Counter(r['change'] for r in result)),categories=dict(Counter(r['error_category'] for r in result)),reference_missing=sum(not r['truth_represented'] for r in result),compatible_set_sizes=dict(Counter(str(r['compatible_pairs']) for r in result)),description='Descriptive post-evaluation attribution, preserving the frozen model. Categories are ordered; multiple failure causes can coexist.')
    (ROOT/'results/validation_case_summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
if __name__=='__main__':main()

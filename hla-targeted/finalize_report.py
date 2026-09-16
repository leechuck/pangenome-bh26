#!/usr/bin/env python3
"""Integrate completed secondary controls into the report; changes no predictions."""
import csv,json,re
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent

def read(p):return list(csv.DictReader(open(p),delimiter='\t'))
def main():
    assert json.loads((ROOT/'results/completion_checks.json').read_text())['status']=='passed'
    specific=json.loads((ROOT/'results/validation_MHC_specificity.json').read_text());assert specific['haplotypes']==48 and specific['haplotype_group_combinations_with_extra_matches']==0
    metrics=read(ROOT/'results/validation_metrics.tsv');cal=read(ROOT/'results/validation_C4Investigator_calibrated_total.tsv');assert len(cal)==24
    def get(method,feature):return next(r for r in metrics if r['method']==method and r['feature']==feature)
    target=get('targeted_markers','total');rank=get('targeted_dosage_paths','full_signature_pair');old=get('old_sketch_CN','full_signature_pair');native=get('C4Investigator','total');correct=sum(int(r['correct']) for r in cal);called=sum(r['call']!='' for r in cal)
    report=(ROOT/'REPORT.md').read_text()
    report=report.replace('All 24 selected families were excluded from assay discovery and candidate reference paths.', 'All 24 selected families were excluded from new targeted-probe discovery, training statistics and candidate reference paths. The inherited global marker vocabulary is distinguished from its training-only eligibility filters in [the marker provenance audit](MARKER_PROVENANCE.md).')
    report=report.replace('Large sequence and read-measurement inputs remain excluded from Git; remote paths and regeneration procedures are in `README.md` and the Slurm scripts.', 'Raw sequence inputs and reads remain outside Git. Compact derived profiles and read measurements are bundled in `source/inference_inputs.tar.gz`, with a verified byte-identical replay recorded in `results/replay_check.json`. The separate held-out MHC background archive supports the specificity control; full sequence-level checks require the original sequence inputs. Paths and regeneration procedures are in `README.md` and the Slurm scripts.')
    report=report.replace('Read-level support and assembly-label discordance must be interpreted separately.', 'NA19700 has 16 deduplicated fragments supporting the noncanonical motif in the development data; the other two candidates were not read-tested here. See [the candidate record](NONCANONICAL_CANDIDATES.md). Read support does not establish novelty, function or physical phase.')
    report=report.replace('a ranking gain does not establish physical phase or a graph-specific advantage', 'ranked imputation does not establish physical phase or a graph-specific advantage')
    case=json.loads((ROOT/'results/validation_case_summary.json').read_text());changes=case['changes']
    report=re.sub(r'Validation superpopulation counts: \{[^}]+\}\.', 'The validation set contains nine AFR, seven AMR, four EAS and four SAS donors.', report)
    disposition=f"Compared with the old sketch, {changes.get('unchanged',0)} donors are unchanged, {changes.get('gain',0)} improve and {changes.get('loss',0)} worsen. The three errors comprise one incorrect marginal dosage, one wrong ranking despite a compatible truth, and one absent reference structure."
    report=re.sub(r'Paired structural changes: \{[^}]+\}\. Post-evaluation error categories: \{[^}]+\}\.',disposition,report)
    # Idempotent integration after the frozen report generator has run.
    begin='<!-- completed-controls:start -->';end='<!-- completed-controls:end -->'
    if begin in report:report=report[:report.index(begin)]+report[report.index(end)+len(end):]
    lead=f'''{begin}
**Main finding:** targeted diagnostic constraints did not improve held-out structural-pair accuracy: {rank['correct']}/24 versus {old['correct']}/24 with the old sketch. Total C4 dosage was {target['correct']}/24, tied with calibrated ordinary reference depth. Module-start context probes added no accuracy in this set. These are new-donor results within the original assembly panel, not external-cohort or clinical validation.

The completed native C4Investigator total-copy result is {native['correct']}/24; its separately declared pilot-calibrated total-copy result is **{correct}/24**. Native component comparisons can inherit total-depth bias and are not evidence of general superiority of the targeted method. See [the calibration supplement](COMPARATOR_CALIBRATION.md) and [the configuration check](COMPARATOR_METHOD_CHECK.md).

The three structural errors illustrate incorrect A/B dosage, unresolved linkage of A/B with long/short forms, and an absent reference structure. Every donor retains multiple compatible structures. See [the actual phase counterexamples](STRUCTURAL_ERRORS.md). All 91 validation genes contain their full diagnostic probe sets, and the whole-MHC control finds no extra diagnostic-probe matches outside those genes in any of the 48 validation haplotypes. This excludes those specific explanations for the A/B error, not all sampling or mapping effects.
{end}
'''
    title,body=report.split('\n',1);report=title+'\n\n'+lead+'\n'+body.lstrip()
    table_row=f'| C4Investigator, pilot-calibrated total | total | {correct}/24 | {called} |'
    lines=[s for s in report.splitlines() if not s.startswith('| C4Investigator, pilot-calibrated total |')]
    for i,s in enumerate(lines):
        if s.startswith('| C4Investigator | total |'):lines.insert(i+1,table_row);break
    (ROOT/'REPORT.md').write_text('\n'.join(lines)+'\n')
    fig,axes=plt.subplots(1,2,figsize=(12,5),layout='constrained');features=['total','A','B','L','S'];x=np.arange(len(features))
    for offset,method,label,color in [(-.26,'targeted_markers','Targeted, calibrated','#2166ac'),(0,'C4Investigator','C4Investigator, native','#b35806')]:axes[0].bar(x+offset,[int(get(method,f)['correct']) for f in features],.24,label=label,color=color)
    axes[0].bar([.26],[correct],.24,label='C4Investigator, calibrated total only',color='#1b7837')
    axes[0].set_xticks(x,features);axes[0].set_ylim(0,25);axes[0].set_ylabel('Correct calls / 24 attempted');axes[0].set_title('Marginal dosage');axes[0].legend(fontsize=7,loc='upper left',bbox_to_anchor=(0,-.12),frameon=False)
    methods=['old_sketch_CN','targeted_dosage_paths','targeted_without_module_context'];vals=[int(get(m,'full_signature_pair')['correct']) for m in methods]
    axes[1].bar(range(3),vals,color=['#888888','#2166ac','#67a9cf']);axes[1].set_xticks(range(3),['Old sketch','Targeted','Without module\ncontexts']);axes[1].set_ylim(0,25);axes[1].set_title('Structural-pair imputation');axes[1].set_ylabel('Exact signature pairs / 24')
    for i,v in enumerate(vals):axes[1].text(i,v+.2,str(v),ha='center',fontsize=9)
    fig.suptitle('Held-out donors within the assembly panel; structural ranking is not physical phase',fontsize=11)
    for ext in ['png','pdf']:fig.savefig(ROOT/f'results/validation_summary.{ext}',dpi=180)
    plt.close(fig)
    audit='''# Completion audit

The following requirements were checked against the recorded artifacts. This audit establishes completion of this experiment, not clinical accuracy or a global novelty claim.

| Requirement | Evidence | Disposition |
|---|---|---|
| A/B, long/short and module-context typing | Frozen probe definitions, dosage model and validation predictions | Complete |
| Development-only use of prior 106 donors | Selection manifest, calibration inputs and family-exclusion checks | Complete |
| New read validation | 24 new donors; selected families absent from 543 reference paths; 24 complete measurement sets | Complete within the existing panel; not external-cohort validation |
| Preserve ambiguity and distinguish phase | Complete compatible-pair sets, imputation labels, structural error counterexamples | Complete; no unique dosage-only structure in this set |
| Dedicated comparator | Native outputs for 18 development and 24 validation donors; independently specified total-dosage calibration | Complete |
| Module-context contribution | Frozen ablation evaluated on all 24 donors | Complete; no accuracy gain |
| Sequence-label and probe controls | 91-gene diagnostic audit, complete probe coverage, 48-haplotype whole-MHC specificity audit | Complete within MHC; whole-genome specificity unproven |
| Reproducibility | Frozen file hashes, training-supported inherited-feature audit, software manifests, 934-file input archive, byte-identical replay of 24 dosage and 72 path rows | Complete for derived-input inference; raw-read and sequence reconstruction require original inputs |
| Reporting | Main report, calibration/configuration supplements, error/candidate records and standalone figures | Complete |

`results/completion_checks.json` verifies the primary invariants, while `results/replay_check.json` and `results/validation_MHC_specificity.json` record the additional controls. The named reports explain the limits of each check. UK Biobank phenotype validation remains a later study, as agreed at the outset; no disease association or validated novel named allele is claimed here.
'''
    (ROOT/'COMPLETION_AUDIT.md').write_text(audit)
    print('Integrated completed controls, regenerated figures, and wrote completion audit')
if __name__=='__main__':main()

import copy
import unittest
from hgsvc_handoff import LOCATIONS,verify_reports


class HandoffTest(unittest.TestCase):
    def fixture(self):
        names=['run_hgsvc_baseline.py','run_linear_control.py','run_hgsvc_genomic.py','run_graph_pair_native_locus.py','graph_pair_refinement_native_locus.py','graph_pair_refinement_alignment_only.py','graph_pair_refinement.py','run_hgsvc_anchor.py','anchor_t1k_coarse.py']
        plan=dict(code_sha256={n:n for n in names},asset_sha256={'t1k-references/genome-v2/COMPLETE.json':'genref','references/observed-v2/COMPLETE.json':'graphref'},metadata_sha256={'cohort.tsv':'cohort'},parameters=dict(temperature=10,locus_margin=10,minimum_fragments=20,minimum_gap=10,noise=.01,internal_pairs=True,pair_specific=True))
        reads={'r1.fq.gz':'read1','r2.fq.gz':'read2'}
        reports=[]
        for method in LOCATIONS:
            r=dict(donor='D',reads_sha256=reads)
            if method=='T1K':r.update(plan_sha256='bp',recipe={'tool':'t1k'},driver_sha256='run_hgsvc_baseline.py')
            elif method=='raw_genomic':r.update(validation_plan_sha256='gp',tool_sha256='tool',driver_sha256='run_linear_control.py',validation_driver_sha256='run_hgsvc_genomic.py',reference_manifest_sha256='genref',cohort_sha256='cohort')
            elif method.startswith('raw_'):
                r.update(baseline_sha256='raw_genomic',parameters=dict(plan['parameters'],alignment_only=True,native_locus_guard=True),panel=method.removeprefix('raw_'),driver_sha256='run_graph_pair_native_locus.py',model_sha256='graph_pair_refinement_native_locus.py',alignment_model_sha256='graph_pair_refinement_alignment_only.py',base_model_sha256='graph_pair_refinement.py',reference_manifest_sha256='graphref')
            else:r.update(original_sha256='T1K',candidate_sha256='raw_genomic' if method=='ipd_genome' else 'raw_'+method.removeprefix('graph_'),driver_sha256='run_hgsvc_anchor.py',anchor_sha256='anchor_t1k_coarse.py',method=method)
            reports.append(dict(donor='D',method=method,status='complete',record=r,manifest_sha256=method))
        return reports,(['D'],plan,{'D':reads},{'recipe':{'tool':'t1k'}},{'tool_sha256':'tool'},'bp','gp')

    def test_complete_chain_and_incomplete_grid(self):
        reports,args=self.fixture()
        self.assertTrue(verify_reports(reports,*args))
        with self.assertRaises(ValueError):verify_reports(reports[:-1],*args)
        reports[-1]['status']='running'
        self.assertFalse(verify_reports(reports,*args))

    def test_wrong_parent_reads_and_missing_guard_rejected(self):
        reports,args=self.fixture()
        for method,key,value in [('graph_hprc_asian','candidate_sha256','wrong'),('T1K','reads_sha256',{'1':'wrong','2':'read2'}),('raw_hprc','parameters',args[1]['parameters'])]:
            bad=copy.deepcopy(reports)
            next(r for r in bad if r['method']==method)['record'][key]=value
            with self.assertRaises(ValueError):verify_reports(bad,*args)


if __name__=='__main__':unittest.main()

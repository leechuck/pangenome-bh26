import tempfile
import unittest
import csv
import json
from pathlib import Path
from build_panel import FLANK, GENES, build, digest, fasta, group_paths, validate_sequence


class ReferenceTests(unittest.TestCase):
    def test_build_excludes_test_aliases_before_using_sequences_or_labels(self):
        with tempfile.TemporaryDirectory() as directory:
            p = Path(directory)
            def tsv(name, records):
                file = p/name
                with file.open('w') as stream:
                    w = csv.DictWriter(stream, fieldnames=list(records[0]), delimiter='\t')
                    w.writeheader(); w.writerows(records)
                return file
            meta = tsv('metadata.tsv', [dict(hap_id='train#1',cohort='HPRC-Rest'),
                                         dict(hap_id='alias#2',cohort='APR')])
            tsv('excluded_paths.tsv',[dict(hap_id='alias#2')])
            tsv('cohort.tsv',[dict(donor='held_out')])
            (p/'RESERVATION.json').write_text(json.dumps(dict(
                cohort_sha256=digest((p/'cohort.tsv').read_bytes()),
                excluded_paths_sha256=digest((p/'excluded_paths.tsv').read_bytes()),
                source_sha256=dict(metadata=digest(meta.read_bytes())))))
            annotations=[]; exons=[]
            sequence='A'*FLANK+'ACGTACGT'+'T'*FLANK
            for gene in GENES:
                name='train#1#HLA-'+gene
                (p/f'HLA-{gene}.fa').write_text('>'+name+'\n'+sequence+'\n>alias#2#HLA-'+gene+'\nINVALID_HELD_OUT_SEQUENCE\n')
                annotations.append(dict(name=name,hap_id='train#1',gene='HLA-'+gene,
                                        gene_sha256=digest(b'ACGTACGT'),cds_sha256=digest(b'ACGTACGT'),
                                        exact_genomic_alleles=gene+'*01:01:01:01',exact_cds_alleles=gene+'*01:01:01:01'))
                annotations.append(dict(annotations[-1],name='alias#2#HLA-'+gene,hap_id='alias#2',gene_sha256='bad'))
                exons.append(dict(name=name,start0=FLANK,end0=FLANK+8))
            cat=tsv('catalogue.tsv',annotations); ex=tsv('exons.tsv',exons)
            result=build(p,meta,cat,ex,p,p/'out')
            self.assertEqual(result['status'],'complete')
            self.assertEqual(len(result['output_sha256']),32)
            side=json.loads((p/'out/hprc_asian/HLA-A.json').read_text())
            self.assertEqual(side[0]['sources'],['train#1#HLA-A'])
            (p/'cohort.tsv').write_text('changed')
            with self.assertRaises(ValueError):
                build(p,meta,cat,ex,p,p/'changed')

    def test_excluded_sequence_is_not_exported(self):
        with tempfile.TemporaryDirectory() as directory:
            p = Path(directory)/'in.fa'
            p.write_text('>test#1#HLA-A\nSECRET\n>train#1#HLA-A\nACGT\n')
            self.assertEqual(list(fasta(p, {'train#1#HLA-A'})), [('train#1#HLA-A','ACGT')])

    def test_observed_gene_and_exon_hashes_are_checked(self):
        seq = 'A'*FLANK + 'ACGTACGT' + 'T'*FLANK
        ann = dict(gene_sha256=digest(b'ACGTACGT'), cds_sha256=digest(b'CGTA'))
        iv = [(FLANK+1,FLANK+5)]
        self.assertEqual(validate_sequence('x',seq,ann,iv),iv)
        with self.assertRaises(ValueError):
            validate_sequence('x',seq,ann,[(FLANK+2,FLANK+6)])
        with self.assertRaises(ValueError):
            validate_sequence('x',seq[:-FLANK]+'G'+'T'*(FLANK-1),
                              dict(ann,gene_sha256='bad'),iv)

    def test_duplicates_do_not_add_read_evidence_or_upgrade_cds_labels(self):
        a = dict(name='a',sequence='ACGT',exons=[(0,4)],genomic_labels=[],cds_labels=['A*01:01:01:01'])
        b = dict(a,name='b')
        g = group_paths([a,b])
        self.assertEqual(len(g),1)
        self.assertEqual(g[0]['sources'],['a','b'])
        self.assertEqual(g[0]['genomic_labels'],[])
        with self.assertRaises(ValueError):
            group_paths([a,dict(b,exons=[(1,4)])])


if __name__ == '__main__':
    unittest.main()

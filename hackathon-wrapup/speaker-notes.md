# Talk track (~3 minutes)

## 1 — Whose HLA is in your reference? (~50 seconds)

“HLA is a good place to ask whether our reference represents the people we want to study. This week we brought five projects together: assemblies and graph paths, through MHC extraction, allele annotation and graph construction.

The result is 754 haplotypes, including 390 Asian or Arab haplotypes. We now have a shared graph and annotation resource, a CWL workflow, and MHC reads from 2,504 samples to work with.

That gives us something we can test: what happens when those extra haplotypes enter the reference panel?”

If asked: 752 biological haplotypes plus GRCh38 and CHM13. Projects are APR, HPRC release 2, JaSaPaGe, K-PanRef and CPC. The drawing is a schematic of paths meeting and diverging. It is not measured graph topology. Data location: `/home/asianhla/data/upload/HLA/` on NIG. Code: https://github.com/leechuck/pangenome-bh26.

## 2 — More haplotypes. Better SV calls. (~55 seconds)

“We took the same reads from 40 donors and used the same caller with two panels: HPRC-only and the expanded panel. The clearest gain was in genotypes carrying structural variants.

In the East Asian donors, exact recovery went from about 72 to 76 percent. In South Asian donors, it went from 75 to 77 percent. SNV gains were much smaller.

This is encouraging pilot evidence. We changed panel size and ancestry together, and measured against assembly-derived truth in the existing graph. Independent validation is the next step.”

If asked: PanGenie; 20 EAS and 20 SAS donors. SV-bearing means a truth allele has a length difference of at least 50 bp from reference. Success requires the entire unordered diploid allele-sequence pair to match. This is genotype recovery at graph bubbles. All eligible truth genotypes form the denominator, including absent predictions. The intervals resample donors. Test families were excluded from training panels; their assemblies remain in graph topology. The HPRC-only panel is drawn from that same topology.

## 3 — Better sequences. Same names. (~65 seconds)

“The surprise came from a small change to an existing tool. We added pangenome sequences to the database SpecHLA uses when collecting reads. Exact gene reconstructions went from 41 to 51 out of 128. The two-field allele calls stayed the same.

So an allele-name benchmark can miss a useful sequence improvement. The sequence itself matters if we want to study variation beyond the familiar name.

In parallel, our Locityper panel experiment still trails T1K on classical allele names, while recovering 111 of 113 DRB3/4/5 type-and-copy-number genotypes. These are separate development experiments.

The next step is locked validation on a much larger set. If you have independent HLA truth or want to test the workflow on another cohort, let's compare notes.”

If asked: SpecHLA sequence experiment has eight development donors, eight genes and two haplotypes per gene. Exact reconstruction is the report's edlib infix comparison to assembled gene bodies; masked N bases count as mismatches. The 40-donor panel comparison is Stage B before exon grafting. T1K 314/319 vs panel 295/319 is the harmonised classical-gene denominator. T1K lacks explicit copy-number reporting, so avoid presenting the DRB3/4/5 contrast as a general typing ranking. The 228/946 validation sets remain future work; successful C1 jobs are not yet scored evidence.

# What the structural errors mean

The frozen targeted model recovers 21/24 full RCCX signature pairs, versus 22/24 with the old sketch. These three mismatches also change the inferred C4 gene-form arrangement; they are not solely CYP21/TNX prototype-label differences.

AL = long C4A; AS = short C4A; BL = long C4B; BS = short C4B. A slash separates the two assembly haplotypes; chromosome ordering within a diploid pair is arbitrary.

| Donor | Assembly C4 paths | Imputed C4 paths | Same A/B and L/S marginals? |
|---|---|---|---|
| HG00673 | AL-AS / AL-BL | AL-BS / AL-AL | Yes |
| HG03540 | AL-BL / AL-BS-BS | BL-AS-BS / AL-BL | Yes |
| HG04204 | AL-BS / AL-BL | AL-BS / AL-AL | No |

**HG00673:** both configurations have A=3, B=1, L=3 and S=1. Correct marginal dosage cannot distinguish them. The true pair remains in the candidate set, but the sequence-profile ranking chooses the wrong combination.

**HG03540:** the true trimodular reference signature is absent after excluding the held-out families. The replacement has the same A/B and L/S marginal counts but different gene forms and order. The caller does not reliably recognize an unknown structure.

**HG04204:** the targeted A/B estimate rounds to A=3, B=1, while the assembly has A=2, B=2. This erroneous constraint excludes the correct pair, which the old sketch recovered. This is the single loss versus the old model; there are no gains in these 24 donors.

Every validation donor has multiple dosage-compatible signature pairs (3–87). Thus none obtains a uniquely resolved full structure from these dosage constraints alone, even though ranking often predicts the assembly label correctly.

A useful extended typing record would keep classical HLA types, C4 total/A/B/L/S dosage, the full set of compatible RCCX structures, and the source of phase evidence as separate fields. Direct long-molecule evidence linking diagnostic sites and module boundaries is the appropriate next test of arrangement; a unique ranking optimum is insufficient.

No caller thresholds or reference paths were changed after observing these errors. These are descriptive analyses of the frozen validation run.

## Post-freeze probe-coverage audit

All 91 validation C4 genes contain all 15 A/B diagnostic probes for their assembly-assigned class. All 60 long genes contain 16 probes at each insertion junction; all 31 short genes contain 16 deletion-junction probes. Thus no assembled diagnostic window is absent from the frozen vocabulary in this set.

HG04204 has 39 A-supporting and 19 B-supporting fragments despite complete probe coverage of both A and both B genes. Probe dropout does not explain this discrepancy. Stochastic sampling or upstream mapping/recruitment effects remain possible; this audit does not establish their cause. See `results/validation_probe_coverage.tsv`.

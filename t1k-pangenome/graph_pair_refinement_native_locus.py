"""Exploratory native-locus veto before four-field graph refinement.

Reject a graph locus when positive native T1K assignments support only other
loci, including non-classical HLA genes missing from the eight-locus graph.
Mixed-locus native evidence and graph-only fragments are not resolved by this
guard. No native and graph scores are treated as comparable probabilities.
"""
import math
from graph_pair_refinement_alignment_only import coarse,internal_fragment,refine_gene
from graph_pair_refinement_alignment_only import assign_fragment as alignment_assign


def native_locus_allowed(record, gene):
    supported=set()
    for allele,assignments in record['native_candidates'].items():
        for assignment in assignments:
            weight=float(assignment['adjusted_weight'])
            if not math.isfinite(weight) or weight<0:
                raise ValueError('Invalid native assignment weight')
            if weight>0:
                supported.add(allele.removeprefix('HLA-').split('*')[0])
    return not supported or gene in supported


def assign_fragment(record, paths, temperature=10, margin=10):
    assigned=alignment_assign(record,paths,temperature,margin)
    if assigned is not None and not native_locus_allowed(record,assigned[0]):
        return None
    return assigned

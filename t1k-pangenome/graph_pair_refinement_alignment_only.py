"""Exploratory alignment-only emissions; frozen v1 inference remains unchanged.

Do not let full observed-path length make identical local alignment evidence
discriminate between alleles. This is an ablation candidate, not calibrated
genotype likelihood or an independently validated improvement.
"""
from graph_pair_refinement import coarse, internal_fragment, refine_gene
from graph_pair_refinement import assign_fragment as _assign_fragment


class _EqualLengthPaths:
    def __init__(self, paths):
        self.paths=paths

    def __contains__(self, path):
        return path in self.paths

    def __getitem__(self, path):
        meta=self.paths[path]
        if meta['effective_length']<=0:
            raise ValueError('Nonpositive effective path length')
        return dict(meta,effective_length=1.)


def assign_fragment(record, paths, temperature=10, margin=10):
    # Adapt only accessed paths, rather than copying the whole panel per read.
    return _assign_fragment(record,_EqualLengthPaths(paths),temperature,margin)

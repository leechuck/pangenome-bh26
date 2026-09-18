"""Experimental joint diploid fitting from relative fragment alignment scores.

One matrix row is one fragment across all competing genes. This optimizes a
score-based mixture objective, not calibrated genotype probabilities. A caller
must supply complete candidate/fallback evidence and preserve tied solutions.
"""
import itertools
import numpy as np


def emission_matrix(fragments, paths, effective_lengths, temperature=10.0):
    if temperature <= 0 or any(effective_lengths[p] <= 0 for p in paths):
        raise ValueError('Positive score temperature and effective lengths required')
    index = {p:i for i,p in enumerate(paths)}
    if len(index)!=len(paths):
        raise ValueError('Duplicate path candidates')
    matrix = np.zeros((len(fragments),len(paths)),dtype=float)
    for row,fragment in enumerate(fragments):
        scores = {p:max(v['alignment_score'] for v in placements)
                  for p,placements in fragment['candidates'].items() if placements}
        if any(p not in index for p in scores):
            raise ValueError('Evidence references an unknown candidate')
        if not scores:
            continue
        maximum = max(scores.values())
        for path,score in scores.items():
            matrix[row,index[path]] = np.exp((score-maximum)/temperature)/effective_lengths[path]
        # A common row factor cancels from genotype comparisons and improves scaling.
        matrix[row] /= matrix[row].max()
    return matrix


def fit(emissions, genes, noise=1e-6, maximum_iterations=100, tolerance=1e-8):
    x = np.asarray(emissions,dtype=float)
    if x.ndim!=2 or x.shape[1]!=len(genes) or not np.isfinite(x).all() or (x<0).any() or (x>1).any():
        raise ValueError('Expected finite fragment-by-candidate scores in [0,1]')
    if not 0 < noise < 1 or maximum_iterations < 1:
        raise ValueError('Invalid fitting parameters')
    groups = {g:[i for i,gene in enumerate(genes) if gene==g and np.any(x[:,i]>0)]
              for g in sorted(set(genes))}
    groups = {g:indices for g,indices in groups.items() if indices}
    if not groups:
        return dict(converged=True,iterations=0,genes={},unassigned_fragments=len(x),objective_trace=[])
    names = list(groups)
    # Initialization is uniform across genes, never based on panel donor counts.
    weights = np.full(len(names),1.0/len(names))
    chosen = {}
    components = []
    for gene in names:
        candidates = groups[gene]
        first = max(candidates,key=lambda i:float(x[:,i].sum()))
        chosen[gene]=(first,first)
        components.append(x[:,first].copy())
    components = np.asarray(components)
    trace=[];converged=False
    for iteration in range(maximum_iterations):
        total = noise+(1-noise)*(weights@components)
        for j,gene in enumerate(names):
            base = np.maximum(noise,total-(1-noise)*weights[j]*components[j])
            best_score=-np.inf;best_pair=None;best_component=None
            for a,b in itertools.combinations_with_replacement(groups[gene],2):
                component=(x[:,a]+x[:,b])*0.5
                score=float(np.log(base+(1-noise)*weights[j]*component).sum())
                if score > best_score+tolerance:
                    best_score,best_pair,best_component=score,(a,b),component
            chosen[gene]=best_pair;components[j]=best_component
            total=base+(1-noise)*weights[j]*best_component
        responsibilities=(1-noise)*weights[:,None]*components/total[None,:]
        mass=responsibilities.sum(axis=1)
        weights=mass/mass.sum()
        score=float(np.log(noise+(1-noise)*(weights@components)).sum())
        if trace and score < trace[-1]-tolerance:
            raise ArithmeticError('Joint objective decreased')
        trace.append(score)
        if len(trace)>1 and abs(trace[-1]-trace[-2])<=tolerance:
            converged=True;break
    total=noise+(1-noise)*(weights@components)
    results={}
    # Recompute tied alternatives under the final gene weights and other genes.
    for j,gene in enumerate(names):
        base=np.maximum(noise,total-(1-noise)*weights[j]*components[j])
        scored=[]
        for a,b in itertools.combinations_with_replacement(groups[gene],2):
            score=float(np.log(base+(1-noise)*weights[j]*(x[:,a]+x[:,b])*0.5).sum())
            scored.append((score,(a,b)))
        best=max(s for s,_ in scored)
        ties=[list(pair) for s,pair in scored if best-s<=tolerance]
        inferior=[s for s,_ in scored if best-s>tolerance]
        results[gene]=dict(pairs=ties,gene_fraction=float(weights[j]),
                           score_gap=best-max(inferior) if inferior else None)
    return dict(converged=converged,iterations=len(trace),genes=results,
                unassigned_fragments=int(np.sum(~np.any(x>0,axis=1))),objective_trace=trace,
                limitations='Coordinate-ascent local solution; score gaps are not calibrated confidence. '
                            'Candidate completeness and all-IPD fallback are caller responsibilities.')

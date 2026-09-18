"""Development-only graph refinement within native T1K two-field calls.

One uniquely assigned fragment contributes once. Observed path multiplicity is
not an allele prior. Scores are exploratory, not calibrated probabilities.
"""
import collections
import itertools
import math
import numpy as np


def coarse(allele):
    return ':'.join(allele.removeprefix('HLA-').split(':')[:2])


def internal_fragment(record, paths, gene):
    """Require both mates inside the exon envelope at every winning-locus placement.

    Reject the whole fragment on ambiguous boundary placement, rather than
    deleting individual candidate placements and creating artificial absences.
    Locus competition must already have used all unfiltered evidence.
    """
    seen=False
    for source in record['graph_sources'].values():
        for path,placements in source['candidates'].items():
            meta=paths[path]
            if meta['gene']!=gene:continue
            for placement in placements:
                seen=True
                span=meta.get('gene_span')
                if span is None:return False
                for mate in ('first','second'):
                    start,end=placement[mate+'_start'],placement[mate+'_end']
                    if not span[0]<=start<end<=span[1]:return False
    return seen


def refine_gene(native_pair, candidates, rows, minimum_fragments=20, minimum_gap=10,
                noise=.01, tolerance=1e-8):
    """Candidate entries carry label (or None) and possible two-field families."""
    if len(native_pair)!=2 or any(not a or len(a.split(':'))!=4 for a in native_pair):
        return dict(changed=False,reason='Native four-field pair unresolved')
    if not 0<noise<1 or minimum_fragments<1 or minimum_gap<0:
        raise ValueError('Invalid refinement thresholds')
    labels={r['label']:i for i,r in enumerate(candidates) if r['label'] is not None}
    if len(labels)!=sum(r['label'] is not None for r in candidates):
        raise ValueError('Allele candidates must be deduplicated')
    if any(a not in labels for a in native_pair):
        return dict(changed=False,reason='Native allele absent from observed graph; retain IPD fallback')
    x=np.asarray(rows,dtype=float)
    if x.size==0 or len(x)<minimum_fragments:
        return dict(changed=False,reason='Insufficient uniquely assigned graph fragments')
    if x.ndim!=2 or x.shape[1]!=len(candidates) or not np.isfinite(x).all() or (x<0).any() or (x>1).any():
        raise ValueError('Invalid fragment emissions')
    family=[coarse(a) for a in native_pair]
    def compatible(a,b):
        fa,fb=candidates[a]['families'],candidates[b]['families']
        return ((family[0] in fa and family[1] in fb) or
                (family[1] in fa and family[0] in fb))
    def score(pair):return float(np.log(noise+(1-noise)*(x[:,pair[0]]+x[:,pair[1]])*.5).sum())
    scored=[(score((a,b)),(a,b)) for a,b in itertools.combinations_with_replacement(range(len(candidates)),2) if compatible(a,b)]
    if not scored:raise ValueError('Native pair is incompatible with its own coarse labels')
    best=max(s for s,_ in scored);ties=[p for s,p in scored if best-s<=tolerance]
    native_score=score(tuple(labels[a] for a in native_pair))
    informative=int(np.sum(np.ptp(x,axis=1)>tolerance))
    result=dict(changed=False,fragments=len(x),informative_fragments=informative,
                score_gain_over_native=best-native_score,tied_pairs=len(ties))
    if len(ties)!=1:return dict(result,reason='Graph genotype is tied')
    pair=ties[0];proposal=[candidates[i]['label'] for i in pair]
    if any(a is None for a in proposal):return dict(result,reason='Best graph pair includes an unlabelled path')
    if sorted(proposal)==sorted(native_pair):return dict(result,reason='Graph agrees with native pair')
    if informative<minimum_fragments or best-native_score<minimum_gap:
        return dict(result,reason='Insufficient discriminating graph evidence')
    inferior=[s for s,p in scored if p!=pair]
    gap=best-max(inferior) if inferior else math.inf
    if gap<minimum_gap:return dict(result,reason='Competing graph pair remains close',score_gap=gap)
    pair_x=x[:,list(pair)];totals=pair_x.sum(axis=1)
    support=np.divide(pair_x,totals[:,None],out=np.zeros_like(pair_x),where=totals[:,None]>0).sum(axis=0)
    if np.any(support<5):return dict(result,reason='Insufficient support for a proposed allele')
    return dict(result,changed=True,reason='Graph-supported four-field refinement',pair=proposal,
                score_gap=gap,effective_allele_support=support.tolist())


def assign_fragment(record,paths,temperature=10,margin=10):
    """Return one gene's relative candidate scores, or reject competing loci."""
    if temperature<=0 or margin<0:raise ValueError('Invalid fragment scoring parameters')
    scores={}
    for source in record['graph_sources'].values():
        for path,placements in source['candidates'].items():
            if path not in paths:raise ValueError('Unknown graph path')
            if placements:
                score=max(p['alignment_score'] for p in placements)
                if not math.isfinite(score):raise ValueError('Nonfinite alignment score')
                scores[path]=max(scores.get(path,-math.inf),score)
    if not scores:return None
    by_gene=collections.defaultdict(list)
    for path,score in scores.items():by_gene[paths[path]['gene']].append(score)
    ranked=sorted(((max(v),g) for g,v in by_gene.items()),reverse=True)
    if len(ranked)>1 and ranked[0][0]-ranked[1][0]<margin:return None
    top,gene=ranked[0]
    values={}
    for path,score in scores.items():
        meta=paths[path]
        if meta['gene']!=gene:continue
        if meta['effective_length']<=0:raise ValueError('Nonpositive effective path length')
        value=math.exp((score-top)/temperature)/meta['effective_length']
        for candidate in meta['candidates']:
            values[candidate]=max(values.get(candidate,0),value)
    maximum=max(values.values(),default=0)
    return (gene,{a:v/maximum for a,v in values.items()}) if maximum else None

"""Refine four-field calls only within original T1K's resolved two-field pair.

This truth-free composition is an experimental candidate. It preserves original
T1K's coarse predictions exactly, including ambiguity and missing-call handling.
"""
import copy
from graph_pair_refinement import coarse


def resolved_pair(call, fields):
    slots=call.get('resolutions',{}).get(str(fields),[])
    if len(slots)!=2:
        return None
    if any(s.get('unresolved',True) or len(s.get('alleles',[]))!=1 for s in slots):
        return None
    pair=[s['alleles'][0] for s in slots]
    return pair if all(len(a.split(':'))==fields for a in pair) else None


def anchor(original, candidate):
    updated=copy.deepcopy(original);decisions={}
    for gene,call in original.items():
        native=resolved_pair(call,2)
        proposed=resolved_pair(candidate.get(gene,{}),4)
        if native is None:
            decisions[gene]='Original two-field pair unresolved; retain original call'
            continue
        if proposed is None:
            decisions[gene]='Candidate four-field pair unresolved; retain original call'
            continue
        families=[coarse(a) for a in proposed]
        if sorted(families)!=sorted(native):
            decisions[gene]='Candidate conflicts with original two-field pair; retain original call'
            continue
        slots=copy.deepcopy(candidate[gene]['resolutions']['4'])
        if families!=native:
            slots.reverse();proposed.reverse()
        updated[gene]['resolutions']['4']=slots
        updated[gene]['copy_count']=1 if proposed[0]==proposed[1] else 2
        decisions[gene]='Compatible four-field refinement within original T1K pair'
    assert set(updated)==set(original)
    for gene,call in original.items():
        assert updated[gene].get('resolutions',{}).get('2')==call.get('resolutions',{}).get('2')
    return updated,decisions

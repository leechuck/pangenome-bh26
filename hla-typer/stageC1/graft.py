#!/usr/bin/env python3
"""Graft a coding sequence onto a panel backbone haplotype.

The backbone supplies everything outside the CDS (UTRs, introns, flanks); the new CDS replaces the
backbone's CDS base by base along a global alignment. Insertions relative to the backbone go into
the exon that holds the preceding aligned base (the first exon for leading insertions). Intervals
are 0-based half-open in genomic (locus) orientation; the CDS is in transcript orientation.
"""
import edlib
COMP=str.maketrans('ACGTN','TGCAN')

def rc(s):return s.translate(COMP)[::-1]

def cds_of(seq,intervals,strand):
    s=''.join(seq[a:b] for a,b in intervals)
    return rc(s) if strand=='-' else s

def column_ops(old,new):
    """Global alignment old->new as a list of (old_index or None, new_base or '')."""
    aln=edlib.align(new,old,mode='NW',task='path')
    ops=[];i=j=0;num=''
    for ch in aln['cigar']:
        if ch.isdigit():num+=ch;continue
        n=int(num);num=''
        for _ in range(n):
            if ch in '=X':ops.append((i,new[j]));i+=1;j+=1
            elif ch=='I':ops.append((None,new[j]));j+=1       # base in new (query) only
            elif ch=='D':ops.append((i,''));i+=1              # base in old (target) only
    assert i==len(old) and j==len(new)
    return ops

def graft(seq,intervals,strand,new_cds):
    """Return (grafted locus sequence, new CDS intervals). intervals: genomic-order list of (start,end)."""
    intervals=sorted(intervals)
    if strand=='-':
        # work in transcript orientation, then map back
        L=len(seq);t_int=[(L-b,L-a) for a,b in reversed(intervals)]
        s,iv=graft(rc(seq),t_int,'+',new_cds)
        M=len(s);return rc(s),sorted((M-b,M-a) for a,b in iv)
    old=cds_of(seq,intervals,'+')
    # old CDS index -> (exon number, locus position)
    where=[(k,p) for k,(a,b) in enumerate(intervals) for p in range(a,b)]
    per_exon=[[] for _ in intervals];last=0
    for oi,base in column_ops(old,new_cds):
        if oi is not None:last=where[oi][0]
        if base:per_exon[last].append(base)
    out=[];prev=0;new_iv=[];pos=0
    for k,(a,b) in enumerate(intervals):
        out.append(seq[prev:a]);pos+=a-prev
        ex=''.join(per_exon[k]);out.append(ex);new_iv.append((pos,pos+len(ex)));pos+=len(ex);prev=b
    out.append(seq[prev:])
    new_seq=''.join(out)
    assert cds_of(new_seq,new_iv,'+')==new_cds
    return new_seq,[x for x in new_iv if x[1]>x[0]]

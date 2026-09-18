"""Exact genomic truth from the documented HGSVC3 sense-oriented allele archive.

No predicted calls or archive-provided allele labels are used. Import is inert.
"""
import re
from collections import defaultdict

GENES = ('A','B','C','DPA1','DPB1','DQA1','DQB1','DRB1')
IDENTIFIER = re.compile(r'^(HG\d{5}|NA\d{5})\.([12])_HLA-([A-Z0-9]+)(?:_|$)')
ALLELE = re.compile(r'^([A-Z0-9]+)\*((?:\d+:)*\d+)([A-Z]?)$')


def fasta(stream):
    name = None
    parts = []
    for line in stream:
        if line.startswith('>'):
            if name is not None:
                yield name, ''.join(parts).upper()
            name, parts = line[1:].strip(), []
        elif line.strip():
            if name is None:raise ValueError('Sequence before FASTA header')
            parts.append(line.strip())
    if name is not None:
        yield name, ''.join(parts).upper()


def reference_index(records, gene):
    index = defaultdict(set)
    for header, sequence in records:
        fields = header.split()
        if len(fields) < 2:raise ValueError('Missing IPD allele label')
        label = fields[1]
        parsed = ALLELE.fullmatch(label)
        if parsed is None:raise ValueError('Invalid IPD allele name')
        if parsed[1] != gene:continue
        if sequence and set(sequence) <= set('ACGT'):
            index[sequence].add(label)
    return {s:sorted(names) for s,names in index.items()}


def derive(records, donors, references):
    """Return exactly two QC records for every donor/gene, including missing loci."""
    donors = tuple(donors)
    if not donors or len(set(donors)) != len(donors):raise ValueError('Empty or duplicate donors')
    observed = defaultdict(list)
    for header, sequence in records:
        match = IDENTIFIER.match(header.split()[0])
        if match and match[1] in donors and match[3] in GENES:
            observed[match[1],match[3],match[2]].append(sequence)
    result = {}
    for donor in donors:
        for gene in GENES:
            slots = []
            for hap in ('1','2'):
                sequences = observed[donor,gene,hap]
                row = dict(haplotype=hap, status='missing', exact_genomic_alleles=[], core_length=None)
                if len(sequences) > 1:
                    row['status'] = 'duplicate_haplotype_gene'
                elif sequences:
                    sequence = sequences[0]
                    core = sequence[1000:-1000]
                    row['core_length'] = len(core)
                    if not core:
                        row['status'] = 'empty_core'
                    elif set(core) - set('ACGT'):
                        row['status'] = 'ambiguous_core'
                    else:
                        row['exact_genomic_alleles'] = list(references[gene].get(core, []))
                        row['status'] = 'exact_genomic' if row['exact_genomic_alleles'] else 'no_exact_genomic_match'
                slots.append(row)
            result[donor,gene] = slots
    return result


def truth_slots(records, fields):
    """Unphased diploid truth; suffixes omitted at numeric-field endpoints."""
    if fields not in (2,4):raise ValueError('Unsupported endpoint')
    if len(records) != 2:raise ValueError('Expected two haplotypes')
    slots = []
    for row in records:
        names = row['exact_genomic_alleles']
        if row['status'] != 'exact_genomic' or not names:return None
        projected = set()
        for name in names:
            match = ALLELE.fullmatch(name)
            if match is None:raise ValueError('Invalid truth allele')
            parts = match[2].split(':')
            if len(parts) < fields or (fields == 4 and len(parts) != 4):return None
            projected.add(match[1]+'*'+':'.join(parts[:fields]))
        slots.append(frozenset(projected))
    return tuple(slots)

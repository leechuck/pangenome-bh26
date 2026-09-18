"""Apply development graph pair refinement with native fallback and full provenance."""
import argparse
import collections
import copy
import hashlib
import json
import os
from pathlib import Path
import time
from evidence_io import sha
from build_graph import sequences
from build_panel import GENES
from t1k_reference import resolve_alias
from graph_pair_refinement import coarse,assign_fragment,refine_gene
from run_development_mapping import graph_path


def native_pair(call):
    slots=call.get('resolutions',{}).get('4',[])
    if len(slots)!=2 or any(s['unresolved'] or len(s['alleles'])!=1 for s in slots):return []
    return [s['alleles'][0] for s in slots]


def candidate_metadata(gene,records,seqs,pair,insert):
    families=sorted({coarse(a) for a in pair}) or [gene+'*UNKNOWN']
    candidates={};paths={}
    if {r['path'] for r in records}!=set(seqs):raise ValueError('Reference metadata/path mismatch')
    for record in records:
        path=record['path'];sequence=seqs[path]
        if hashlib.sha256(sequence.encode()).hexdigest()!=record['sequence_sha256']:
            raise ValueError('Path sequence differs from label metadata')
        if len(sequence)-insert+1<=0:raise ValueError('Path shorter than calibrated insert')
        resolved=resolve_alias(path,{path:record},4)
        ids=[]
        for allele in resolved['alleles']:
            if not allele.startswith(gene+'*'):raise ValueError('Cross-gene reference label')
            if coarse(allele) in families:
                candidates.setdefault(allele,dict(id=allele,label=allele,families=[coarse(allele)]))
                ids.append(allele)
        # Every path remains represented. Unknown or incompatible genomic labels
        # become nuisance alternatives; CDS identity never creates a four-field call.
        if resolved['unresolved'] or any(coarse(a) not in families for a in resolved['alleles']) or not ids:
            unknown='UNKNOWN:'+path
            candidates[unknown]=dict(id=unknown,label=None,families=families)
            ids.append(unknown)
        paths[path]=dict(gene=gene,effective_length=len(sequence)-insert+1,candidates=ids)
    return paths,[candidates[k] for k in sorted(candidates)]


def run(baseline,joined,reference,calibration,output):
    if not os.environ.get('SLURM_CPUS_PER_TASK'):raise ValueError('Run under Slurm')
    if output.exists():raise FileExistsError(output)
    bm=json.loads((baseline/'COMPLETE.json').read_text())
    jm=json.loads((joined/'COMPLETE.json').read_text())
    cm=json.loads(calibration.read_text())
    rm=json.loads((reference.parent/'COMPLETE.json').read_text())
    if any(r['status']!='complete' for r in (bm,jm,cm,rm)) or not cm['estimate']['usable']:
        raise ValueError('Incomplete inputs')
    if sha(baseline/'calls.json')!=bm['output_sha256']['calls.json']:
        raise ValueError('Native calls changed')
    calls=json.loads((baseline/'calls.json').read_text())
    reads=sorted(bm['reads_sha256'].values())
    if sorted(cm['reads_sha256'].values())!=reads:raise ValueError('Calibration reads differ')
    mappings={}
    for chain in jm['graph_chains']:
        for folder,digest in chain.items():
            if sha(Path(folder)/'COMPLETE.json')!=digest:raise ValueError('Graph evidence provenance changed')
            if '/mapping-' in folder:
                mm=json.loads((Path(folder)/'COMPLETE.json').read_text())
                if sorted(mm['reads_sha256'].values())!=reads or mm.get('library_calibration_sha256')!=sha(calibration):
                    raise ValueError('Mapping read/calibration mismatch')
                mappings[mm['gene']]=mm
    if set(mappings)!=set(GENES):raise ValueError('Missing competing graph loci')
    paths={};candidates={};pairs={}
    for gene in GENES:
        fasta=reference/('HLA-'+gene+'.fa');metadata=fasta.with_suffix('.json')
        for p in (fasta,metadata):
            if sha(p)!=rm['output_sha256'][str(p.relative_to(reference.parent))]:
                raise ValueError('Observed reference changed')
        gm_path=graph_path(reference.name,gene)/'COMPLETE.json'
        gm=json.loads(gm_path.read_text())
        if (gm['status']!='complete' or sha(gm_path)!=mappings[gene]['graph_manifest_sha256'] or
            gm['source_sha256']!=sha(fasta)):
            raise ValueError('Mapping panel differs from labelled reference')
        pairs[gene]=native_pair(calls.get(gene,{}))
        seqs=sequences(fasta)
        pp,cc=candidate_metadata(gene,json.loads(metadata.read_text()),seqs,pairs[gene],cm['estimate']['fragment_mean'])
        if set(paths)&set(pp):raise ValueError('Repeated path across genes')
        paths.update(pp);candidates[gene]=cc
    source=joined/'evidence.jsonl'
    if sha(source)!=jm['output_sha256']:raise ValueError('Joined fragments changed')
    output.mkdir(parents=True)
    record=dict(status='running',started=time.time(),donor=bm['donor'],reads_sha256=bm['reads_sha256'],
                baseline_sha256=sha(baseline/'COMPLETE.json'),joined_sha256=sha(joined/'COMPLETE.json'),
                reference_manifest_sha256=sha(reference.parent/'COMPLETE.json'),panel=reference.name,
                calibration_sha256=sha(calibration),driver_sha256=sha(Path(__file__)),
                model_sha256=sha(Path(__file__).with_name('graph_pair_refinement.py')),
                parameters=dict(temperature=10,locus_margin=10,minimum_fragments=20,minimum_gap=10,noise=.01),
                scope='Development graph four-field refinement; native two-field calls and unresolved fallbacks retained')
    def save():(output/'manifest.json').write_text(json.dumps(record,indent=2)+'\n')
    save()
    try:
        rows=collections.defaultdict(list);counts=collections.Counter();previous=None
        with source.open() as stream:
            for line in stream:
                item=json.loads(line)
                if previous is not None and item['fragment']<=previous:raise ValueError('Joined fragments duplicated or unsorted')
                previous=item['fragment'];counts['input_fragments']+=1
                assigned=assign_fragment(item,paths)
                if assigned is None:continue
                gene,values=assigned
                rows[gene].append([values.get(c['id'],0) for c in candidates[gene]])
                counts[gene]+=1
        updated=copy.deepcopy(calls);decisions={}
        for gene in GENES:
            result=refine_gene(pairs[gene],candidates[gene],rows[gene]);decisions[gene]=result
            if not result['changed']:continue
            pair=result['pair'];old=pairs[gene]
            if [coarse(a) for a in pair]!=[coarse(a) for a in old]:pair=list(reversed(pair))
            if [coarse(a) for a in pair]!=[coarse(a) for a in old]:raise ValueError('Refinement changed coarse pair')
            updated[gene]['resolutions']['4']=[dict(raw='GRAPH:'+a,alleles=[a],unresolved=False) for a in pair]
            updated[gene]['copy_count']=1 if pair[0]==pair[1] else 2
        for gene in calls:
            if updated[gene]['resolutions']['2']!=calls[gene]['resolutions']['2']:
                raise ValueError('Native two-field calls changed')
        (output/'calls.json').write_text(json.dumps(updated,indent=2)+'\n')
        (output/'decisions.json').write_text(json.dumps(decisions,indent=2)+'\n')
        record.update(status='complete',finished=time.time(),fragment_counts=dict(counts),
                      changed_genes=[g for g,d in decisions.items() if d['changed']],
                      output_sha256={n:sha(output/n) for n in ('calls.json','decisions.json')})
        save();(output/'COMPLETE.json').write_text(json.dumps(record,indent=2)+'\n')
    except BaseException as error:
        record.update(status='failed',finished=time.time(),error=repr(error));save();raise


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('baseline','joined','reference','calibration','output'):p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();run(a.baseline,a.joined,a.reference,a.calibration,a.output)

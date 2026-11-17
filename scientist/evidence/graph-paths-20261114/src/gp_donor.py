"""DEVELOPMENT DONORS ONLY (reads truth for scoring): placement test and haplotype-sampling typer on the
alignment-built graphs, genes DRB1, DQA1, B.
  python gp_donor.py READSET DONOR
Outputs (names, positions, counts, distances, flags only) under work/<readset>/<donor>/:
placement_sites.tsv, sampling.tsv, coverage.json, COMPLETE.json | FAILURE.json.
Sampled walks (node ids only) are kept as GAF; reads subsets stay on IBEX."""
import collections
import json
import re
import shutil
import statistics
import sys
import time
import traceback
from pathlib import Path

import edlib

sys.path.insert(0, str(Path(__file__).resolve().parent))
import gp_lib as G   # noqa: E402
C, GT, PT, RD, L, PR = G.C, G.GT, G.PT, G.RD, G.L, G.PR

# D-sim: giraffe stalls single-threaded estimating the fragment distribution on this graph; give it the
# simulator's fragment model (SIMULATOR_CONTRACT.json: mean 350, sd 50). D-real: estimated from bwa pairs on the
# pass-1 alleles (median, 1.4826*MAD), because giraffe stalls the same way on the small read subsets.
FRAG = {'D-sim': ['--fragment-mean', 350, '--fragment-stdev', 50]}
SAMPLE_TIMEOUT = 300     # seconds per vg haplotypes / vg paths call (normal run < 2 s)
GIRAFFE_TIMEOUT = 900    # seconds per giraffe call (normal run 6-120 s)
MODES = {'dip2': ['--diploid-sampling'], 'hap4': ['--num-haplotypes', 4], 'hap8': ['--num-haplotypes', 8]}


def ed(a, b):
    if a == b:
        return 0
    if not a or not b:
        return max(len(a), len(b))
    return edlib.align(a, b, mode='NW', task='distance')['editDistance']


def body_of(hap_bases, c0, c1):
    return ''.join(b for b, c in hap_bases if c0 <= c <= c1)


def hap_bases(walk, node_seq, node_cols):
    out = []
    for n, o in walk:
        assert o == '+', 'reverse step in sampled path'
        out.extend(zip(node_seq[n], node_cols[n]))
    return out


def count_targets(fqs, targets):
    """targets: {key: set(kmers)} -> {key: Counter over its kmers} (canonical 29-mers)."""
    allk = {}
    for key, ks in targets.items():
        for x in ks:
            allk.setdefault(x, []).append(key)
    cnt = {key: collections.Counter() for key in targets}
    k = 29
    for fq in fqs:
        for _, s in G.read_fastq_seqs(fq):
            s = s.upper()
            r = G.rc(s)
            n = len(s)
            for i in range(n - k + 1):
                x = s[i:i + k]
                y = r[n - i - k:n - i]
                z = x if x < y else y
                keys = allk.get(z)
                if keys:
                    for key in keys:
                        cnt[key][z] += 1
    return cnt


def placement(readset, donor, g, variant, gdir, fq, msa, truth, called_rows, work):
    """GRAPH_TEST.md section 4 measure on this graph: reads carrying the true local version at sites where
    the true copy differs from the nearest called (arm A) allele; placed exactly on a route spelling it?
    Only read pairs in which a mate contains a site's true or called local version (either strand) are given
    to giraffe (the measure counts no other read; mates are kept so pairing still applies)."""
    work.mkdir(parents=True, exist_ok=True)
    gr = GT.Graph(gdir / 'graph.gfa')
    pathset = {s for p, s in gr.pseq.items() if not p.startswith('REF')}
    rows_needed = [t['row'] for t in truth] + called_rows
    pan = PR.Panel([str(r) for r in rows_needed], [msa.aligned[r] for r in rows_needed])
    nT = len(truth)
    named_loci = GT.NAMED if g == 'DRB1' else {}
    sites = []
    for i, t in enumerate(truth):
        Pt = pan.P[i]
        j = min(range(nT, len(rows_needed)), key=lambda j: L.mism(pan.P[j], Pt))
        Pa = pan.P[j]
        for l, rr, a0, b0 in GT.sites(Pt, Pa):
            ts = PT.row_string(Pt, l, rr)
            cu = PT.row_string(Pa, l, rr)
            if ts is None or cu is None:
                continue
            sites.append((i, t, rows_needed[j], a0, b0, ts, cu))
    pats = set()
    for *_, ts, cu in sites:
        if len(ts) <= 120:
            pats |= {ts, G.rc(ts), cu, G.rc(cu)}
    keep = set()
    for f in fq:
        for h, s in G.read_fastq_seqs(f):
            s = s.upper()
            if any(x in s for x in pats):
                keep.add(G.qname(h))
    sub = [work / 'sub_1.fq', work / 'sub_2.fq']
    for m in (0, 1):
        with open(fq[m]) as src, open(sub[m], 'w') as dst:
            while True:
                h = src.readline()
                if not h:
                    break
                rest = [src.readline() for _ in range(3)]
                if G.qname(h) in keep:
                    dst.write(h + ''.join(rest))
    t0 = time.time()
    reads = []
    if keep:
        G.sh([G.VG, 'giraffe', '-Z', gdir / 'graph.gbz', '-d', gdir / 'graph.dist', '-m', gdir / 'graph.min',
              '-z', gdir / 'graph.zipcodes', '-f', sub[0], '-f', sub[1], '-t', 4] + FRAG.get(readset, []), work,
             work / 'mapping.gam', timeout=GIRAFFE_TIMEOUT)
        G.sh([G.VG, 'view', '-aj', 'mapping.gam'], work, work / 'aln.jsonl', timeout=GIRAFFE_TIMEOUT)
        for line in open(work / 'aln.jsonl'):
            a = json.loads(line)
            s = a.get('sequence', '').upper()
            w = []
            for m in a.get('path', {}).get('mapping', []):
                p = m['position']
                ns = gr.seq[(str(p['node_id']), '-' if p.get('is_reverse') else '+')]
                o = int(p.get('offset', 0))
                fl = sum(int(e.get('from_length', 0)) for e in m.get('edit', []))
                w.append(ns[o:o + fl])
            reads.append((s, ''.join(w)))
    secs = round(time.time() - t0, 1)
    for f in ('mapping.gam', 'aln.jsonl', 'sub_1.fq', 'sub_2.fq'):
        (work / f).unlink(missing_ok=True)
    n_mapped = sum(bool(w) for _, w in reads)
    out = []
    for i, t, a_row, a0, b0, ts, cu in sites:
        seq = t['core']
        kind = 'length' if len(ts) != len(cu) else 'substitution'
        named = next((n for n, (c0, c1) in named_loci.items() if a0 <= c1 and b0 >= c0), '')
        in_path = any(ts in s for s in pathset)
        walk = gr.thread(ts) == 0
        ix = seq.find(ts)
        ext = seq[max(0, ix - 60):ix + len(ts) + 60] if ix >= 0 else ''
        win = bool(ext) and any(ext in s for s in pathset)
        win_walk = bool(ext) and gr.thread(ext) == 0
        cnt = collections.Counter()
        if len(ts) <= 120:
            trc, crc = G.rc(ts), G.rc(cu)
            for s, w in reads:
                if ts in s or trc in s:
                    x = ts if ts in s else trc
                    cnt['t'] += 1
                    cnt['t_exact' if x in w else 't_edit' if w else 't_unmapped'] += 1
                elif cu in s or crc in s:
                    cnt['c'] += 1
        out.append(dict(readset=readset, donor=donor, gene=g, variant=variant, copy=t['hap'],
                        cls=t['cls'], ipd_absent=t['absent'], called_row=a_row, col_start=a0, col_end=b0,
                        named=named or '-', kind=kind, len_true=len(ts), len_called=len(cu), testable=len(ts) <= 120,
                        true_in_stored_path=in_path, true_is_walk=walk, reads_true=cnt['t'],
                        reads_true_exact=cnt['t_exact'], reads_true_edit=cnt['t_edit'],
                        reads_true_unmapped=cnt['t_unmapped'], reads_called=cnt['c'],
                        window_in_stored_path=win, window_is_walk=win_walk,
                        reads_given_to_giraffe=len(reads), reads_mapped=n_mapped, giraffe_seconds=secs))
    return out


def main():
    readset, donor = sys.argv[1], sys.argv[2]
    assert readset in ('D-sim', 'D-real')
    dev = {tuple(l.split()) for l in open(G.E / 'mosaic-drb1-20261113/work/dev_tasks.txt') if l.strip()}
    assert (readset, donor) in dev, 'not a development donor'
    out = G.S / 'work' / readset / donor
    out.mkdir(parents=True, exist_ok=True)   # resumable: finished settings live in out/parts
    for f in ('COMPLETE.json', 'FAILURE.json'):
        (out / f).unlink(missing_ok=True)
    parts = out / 'parts'
    parts.mkdir(exist_ok=True)
    failures = []
    t_start = time.time()
    info = dict(readset=readset, donor=donor, times={})
    try:
        records = C.load_records()
        msas = C.load_msas(records, genes=G.GENES3)
        classes = C.load_classes()
        fold = GT.fold_of(donor)
        info['fold'] = fold
        spec = RD.pass1_spec(readset, donor, 'A')
        fq = [spec['run'] / 'calls_candidate_1.fq', spec['run'] / 'calls_candidate_2.fq']
        geno = C.read_t1k_raw(C.verified_table(spec['run']))
        wanted, called = set(), {}
        for g, row in geno.items():
            sl = [s for s in (row['slots'][:row['copies']] if row['copies'] else []) if s['abundance'] > 0 and s['labels']]
            called[g] = [s['labels'][0] for s in sl]
            for s in sl:
                wanted |= set(s['labels'][:3])
        seqs, _ = RD.read_reference(C.GENOMIC_FA, wanted, set())
        target = {n: seqs[n] for n in sorted(wanted) if n in seqs}
        info['called'] = {g: v for g, v in called.items()}
        t = time.time()
        reads, mstats, gene_of_q = L.map_reads2(target, fq[0], fq[1], out / 'map', 4, RD.ENV)
        shutil.rmtree(out / 'map', ignore_errors=True)
        # fragment length, truth-free: bwa pairs with both mates on the same pass-1 allele
        span = collections.defaultdict(list)
        for g_, xs in reads.items():
            for x in xs:
                span[(x['qname'], x['contig'])].append((x['ref_start'], x['ref_end']))
        fr = sorted(max(e for _, e in v) - min(b for b, _ in v) for v in span.values() if len(v) == 2)
        fr = [f for f in fr if 100 <= f <= 2000]
        fmed = statistics.median(fr)
        fsd = 1.4826 * statistics.median([abs(f - fmed) for f in fr])
        info['fragment_estimate'] = dict(pairs=len(fr), median=fmed, robust_sd=round(fsd, 1))
        if readset == 'D-real':
            FRAG['D-real'] = ['--fragment-mean', round(fmed, 1), '--fragment-stdev', round(max(fsd, 10), 1)]
        info['times']['bwa'] = round(time.time() - t, 1)
        info['bwa_stats'] = mstats
        # read input (ii) per gene and k-mer targets for coverage
        inputs = {}
        targets_i, cov = {}, {}
        other_k = {}
        called_k = {n: G.kmers(seqs[n]) for g in called for n in called[g] if n in seqs}
        for g in G.GENES3:
            gd = out / g
            gd.mkdir(exist_ok=True)
            keep = {x['qname'] for x in reads.get(g, []) if not (gene_of_q[x['qname']] - {g})}
            info.setdefault('reads_ii_pairs', {})[g] = len(keep)
            info.setdefault('reads_assigned', {})[g] = len({x['qname'] for x in reads.get(g, [])})
            for m in (1, 2):
                with open(gd / f'ii_{m}.fq', 'w') as f:
                    with open(fq[m - 1]) as src:
                        while True:
                            h = src.readline()
                            if not h:
                                break
                            rest = [src.readline() for _ in range(3)]
                            if G.qname(h) in keep:
                                f.write(h + ''.join(rest))
            inputs[g] = {'i': fq, 'ii': [gd / 'ii_1.fq', gd / 'ii_2.fq']}
            cs = [n for n in called.get(g, []) if n in seqs]
            assert cs, f'no arm-A call for {g}'
            both = set.intersection(*[called_k[n] for n in cs])
            others = set()
            for n, ks in called_k.items():
                if C.gene_of_name(n) != g:
                    others |= ks
            targets_i[g] = both - others
        t = time.time()
        cnt_i = count_targets(fq, targets_i)
        for g in G.GENES3:
            cnt_ii = count_targets(inputs[g]['ii'], {g: targets_i[g]})[g]
            for inp, cn in (('i', cnt_i[g]), ('ii', cnt_ii)):
                v = [cn[x] for x in targets_i[g]]
                med = statistics.median(v) if v else 0
                cov[(g, inp)] = max(1, int(round(med)))
                info.setdefault('coverage', {})[f'{g}:{inp}'] = dict(
                    target_kmers=len(v), median=med, mean=round(sum(v) / len(v), 2) if v else 0,
                    zero_fraction=round(sum(x == 0 for x in v) / len(v), 3) if v else 1, used=cov[(g, inp)])
        info['times']['coverage'] = round(time.time() - t, 1)
        # kmc counts
        kff = {}
        t = time.time()
        for key, files, d in [(('*', 'i'), fq, out / 'kmc_i')] + [((g, 'ii'), inputs[g]['ii'], out / g / 'kmc_ii') for g in G.GENES3]:
            shutil.rmtree(d, ignore_errors=True)
            d.mkdir(parents=True)
            (d / 'tmp').mkdir()
            (d / 'reads.list').write_text('\n'.join(map(str, files)) + '\n')
            G.sh([G.KMC, '-k29', '-m4', '-okff', '-t4', '-hp', '@reads.list', 'counts', 'tmp'], d, log='kmc.log')
            kff[key] = d / 'counts.kff'
        info['times']['kmc'] = round(time.time() - t, 1)
        place_rows, samp_rows = [], []
        for g in G.GENES3:
            msa = msas[g]
            truth = []
            for r in sorted((r for r in records if r['donor'] == donor and r['gene'] == g), key=lambda r: str(r['haplotype'])):
                core = r['core_sequence'].upper()
                ex = r['core_exons']
                row = msa.row_of_source[r['source_id']]
                bc = msa.base_cols(row)
                truth.append(dict(hap=r['haplotype'], core=core, body=core[ex[0][0]:ex[-1][1]], row=row,
                                  c0=bc[ex[0][0]], c1=bc[ex[-1][1] - 1], absent=not r['genomic_labels'],
                                  cls=classes[r['source_id']]['ipd_label_L150']))
            assert len(truth) == 2, f'{g}: {len(truth)} truth copies'
            called_rows = [msa.row_of_allele[n] for n in called[g] if n in msa.row_of_allele]
            ipd_rows = [r for r in range(len(msa.uids)) if msa.ipd_names[r]]
            ipd_bodies = {}
            for tt in truth:
                ipd_bodies[(tt['c0'], tt['c1'])] = {msa.aligned[r][tt['c0']:tt['c1'] + 1].replace('-', '') for r in ipd_rows}
            for variant in G.VARIANTS:
                gdir = G.graph_dir(g, variant, fold)
                assert (gdir / 'COMPLETE.json').exists(), f'graph missing {gdir}'
                rows_in_graph = G.panel_rows(msa, variant, fold)
                stored_bodies = {}
                for tt in truth:
                    stored_bodies[(tt['c0'], tt['c1'])] = {msa.aligned[r][tt['c0']:tt['c1'] + 1].replace('-', '') for r in rows_in_graph}
                node_seq, node_cols = G.load_nodes(gdir)
                wd = out / g / variant
                wd.mkdir(parents=True, exist_ok=True)
                for inp in ('i', 'ii'):
                    kf = kff[('*', 'i')] if inp == 'i' else kff[(g, 'ii')]
                    for Lc in G.SUBCHAINS:
                        for mode, extra in MODES.items():
                            tag = f'{inp}_L{Lc}_{mode}'
                            pf = parts / f'samp_{g}_{variant}_{tag}.json'
                            if pf.exists():
                                prev = C.load_json(pf)
                                samp_rows += prev['rows']
                                if prev.get('failure'):
                                    failures.append(prev['failure'])
                                continue
                            n_before = len(samp_rows)
                            t = time.time()
                            err = ''
                            try:
                                G.sh([G.VG, 'haplotypes', '-t', 4, '-v', 1, '--coverage', cov[(g, inp)], *extra,
                                      '-i', gdir / f'hapl_{Lc}.hapl', '-k', kf, '-g', f'{tag}.gbz', gdir / 'graph.gbz'],
                                     wd, log=f'{tag}.log', timeout=SAMPLE_TIMEOUT)
                                G.sh([G.VG, 'paths', '-x', f'{tag}.gbz', '-A', '-H'], wd, wd / f'{tag}.gaf', log=f'{tag}.log',
                                     timeout=SAMPLE_TIMEOUT)
                            except Exception as e:
                                err = repr(e)[:200]
                            secs = round(time.time() - t, 1)
                            (wd / f'{tag}.gbz').unlink(missing_ok=True)
                            haps = []
                            if not err:
                                for line in open(wd / f'{tag}.gaf'):
                                    if line.startswith('@') or not line.strip():
                                        continue
                                    f = line.rstrip('\n').split('\t')
                                    haps.append((f[0], hap_bases(G.parse_walk(f[5]), node_seq, node_cols)))
                            m = re.search(r'Total: \d+ haplotypes for \d+ chains \((\d+) full, (\d+) subchains, (\d+) fragments\); '
                                          r'connection rate ([0-9.e-]+); used (\d+) kmers', (wd / f'{tag}.log').read_text())
                            base = dict(readset=readset, donor=donor, gene=g, variant=variant, input=inp,
                                        subchain=Lc, mode=mode, coverage=cov[(g, inp)], n_haps=len(haps),
                                        seconds=secs, error=err, n_subchains=m.group(2) if m else '',
                                        connection_rate=m.group(4) if m else '', kmers_used=m.group(5) if m else '')
                            if err or not haps:
                                for tt in truth:
                                    samp_rows.append(dict(base, copy=tt['hap'], ipd_absent=tt['absent'], cls=tt['cls'],
                                                          body_len=len(tt['body']), dist=len(tt['body']), exact=False,
                                                          best_any_dist=len(tt['body']), exact_any=False,
                                                          picked_named=False, picked_stored=False, rank='', picked_len=0))
                                fl = dict(step='sampling', gene=g, variant=variant, setting=tag, error=err or 'no haplotypes', seconds=secs)
                                failures.append(fl)
                                C.save_json(pf, dict(rows=samp_rows[n_before:], failure=fl))
                                continue
                            bodies = [[body_of(h, tt['c0'], tt['c1']) for _, h in haps] for tt in truth]
                            D = [[ed(tt['body'], b) for b in bodies[i]] for i, tt in enumerate(truth)]
                            if mode == 'dip2':
                                h = list(range(len(haps)))
                                if len(h) == 1:
                                    h = h * 2
                                opts = [(D[0][a] + D[1][b], a, b) for a in h for b in h if a != b or len(set(h)) == 1]
                                _, a, b = min(opts)
                                pick = [a, b]
                            else:
                                pick = [min(range(len(haps)), key=lambda j: D[i][j]) for i in range(2)]
                            for i, tt in enumerate(truth):
                                pb = bodies[i][pick[i]]
                                key = (tt['c0'], tt['c1'])
                                samp_rows.append(dict(base, copy=tt['hap'], ipd_absent=tt['absent'], cls=tt['cls'],
                                                      body_len=len(tt['body']), dist=D[i][pick[i]], exact=D[i][pick[i]] == 0,
                                                      best_any_dist=min(D[i]), exact_any=min(D[i]) == 0,
                                                      picked_named=pb in ipd_bodies[key], picked_stored=pb in stored_bodies[key],
                                                      rank=pick[i], picked_len=len(pb)))
                            C.save_json(pf, dict(rows=samp_rows[n_before:]))
                C.save_json(out / 'partial.json', dict(gene=g, variant=variant, stage='sampling done'))
                pf = parts / f'place_{g}_{variant}.json'
                if pf.exists() and not C.load_json(pf).get('failure'):
                    prev = C.load_json(pf)
                else:
                    t = time.time()
                    try:
                        prev = dict(rows=placement(readset, donor, g, variant, gdir, fq, msa, truth, called_rows,
                                                   out / g / variant / 'place'))
                    except Exception as e:
                        prev = dict(rows=[], failure=dict(step='placement', gene=g, variant=variant, setting='giraffe',
                                                          error=repr(e)[:300], seconds=round(time.time() - t, 1)))
                    prev['seconds'] = round(time.time() - t, 1)
                    C.save_json(pf, prev)
                place_rows += prev['rows']
                if prev.get('failure'):
                    failures.append(prev['failure'])
                info['times'][f'place_{g}_{variant}'] = prev.get('seconds')
                C.save_json(out / 'partial.json', dict(gene=g, variant=variant, stage='placement done'))
        for name, rows in (('placement_sites.tsv', place_rows), ('sampling.tsv', samp_rows)):
            cols = list(rows[0]) if rows else []
            with open(out / name, 'w') as f:
                f.write('\t'.join(cols) + '\n')
                for r in rows:
                    f.write('\t'.join(str(r[c]) for c in cols) + '\n')
        info['seconds'] = round(time.time() - t_start, 1)
        info['failures'] = failures
        C.save_json(out / 'coverage.json', info)
        C.save_json(out / 'failures.json', failures)
        C.save_json(out / 'COMPLETE.json', dict(status='complete', readset=readset, donor=donor, fold=fold,
                                                placement_rows=len(place_rows), sampling_rows=len(samp_rows),
                                                failed_settings=len(failures),
                                                sampling_errors=sum(bool(r['error']) for r in samp_rows),
                                                seconds=info['seconds']))
        print('done', readset, donor, len(place_rows), len(samp_rows), info['seconds'])
    except Exception:
        C.save_json(out / 'FAILURE.json', dict(status='failed', error=traceback.format_exc(), **{k: v for k, v in info.items() if k != 'times'}))
        raise


if __name__ == '__main__':
    main()

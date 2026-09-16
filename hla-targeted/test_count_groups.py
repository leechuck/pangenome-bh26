"""Synthetic fragment contracts independent of assay sequences."""
import csv, random, struct, subprocess, tempfile
from pathlib import Path
ROOT=Path(__file__).resolve().parent

def key(s):
    rc=s.translate(str.maketrans('ACGT','TGCA'))[::-1]
    return min(int(''.join(format('ACGT'.index(c),'02b') for c in x),2) for x in (s,rc))
def main():
    rng=random.Random(126)
    probes=[''.join(rng.choices('ACGT',k=31)) for _ in range(4)]
    with tempfile.TemporaryDirectory() as td:
        d=Path(td); exe=d/'counter'
        subprocess.run(['g++','-O2','-std=c++17',str(ROOT/'count_groups.cpp'),'-o',str(exe)],check=True)
        (d/'markers').write_text(''.join(f'{key(s)}\t{m}\n' for s,m in zip(probes,[1,2,64,16])))
        records=[]
        def add(name,s,q=None): records.append(f'@{name}\n{s}\n+\n{q or "I"*len(s)}\n')
        add('dedup/1',probes[0]+'N'+probes[0]);add('dedup/2',probes[0].translate(str.maketrans('ACGT','TGCA'))[::-1])
        add('mixed/1',probes[0]);add('mixed/2',probes[1])
        add('quality',probes[1],'I'*15+'4'+'I'*15) # Q19 invalidates the only marker
        add('threshold',probes[1],'5'*31) # Q20 accepted
        add('other',probes[2]);add('short',probes[3]);add('ambiguous','N'*31)
        result=subprocess.run([str(exe),str(d/'markers'),str(d/'counts'),str(d/'groups')],input=''.join(records),text=True)
        assert result.returncode==0
        assert struct.unpack('<4I',(d/'counts').read_bytes())==(2,2,1,1)
        groups={r['metric']:int(r['fragments']) for r in csv.DictReader(open(d/'groups'),delimiter='\t')}
        assert groups['all_fragments']==7,groups
        assert [groups[k] for k in ['A_only','B_only','OTHER_only','AB_ambiguous','SHORT']]==[1,1,1,1,1],groups
        bad=subprocess.run([str(exe),str(d/'markers'),str(d/'counts'),str(d/'groups')],input='@broken\nACTG\n+\nI\n',text=True)
        assert bad.returncode==3
    print('PASS: fragment deduplication, reverse complements, mixed mates, Q19/Q20, noncanonical probes, missing bases, malformed FASTQ')
if __name__=='__main__':main()

#!/bin/bash
set -euo pipefail
base=/home/leechuck/hla/dogohla-series20260918
real=/home/leechuck/hla/mm/envs/asian50-spechla
cp /home/leechuck/hla/portable-bin/inchworm /home/leechuck/hla/portable-bin/less /home/leechuck/hla/portable-bin/make "$real/bin/"
export LD_LIBRARY_PATH="$real/lib"
python3 -c 'import Bio, edlib, pysam, scipy, numpy'
command -v make inchworm less samtools bowtie2 makeblastdb fermi2 seqtk bwa perl
python3 "$base/series/configure_series.py"
python3 "$base/series/prepare_series.py" "$base" --threads 8

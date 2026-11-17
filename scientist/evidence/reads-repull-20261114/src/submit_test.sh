#!/bin/bash
# Step 2b test submission (login-node safe: sbatch only). Refuses to run without FREEZE_TEST.json.
# Groups: 40 donors with baseline full-WGS candidates (work/treal_cand40.txt), 26 donors needing wgsx
# (work/treal_need26.txt). Per group: selection -> T1K A, T1K B -> builder C (after A), builder D (after B),
# dependencies task by task (aftercorr). Job ids appended to work/TEST_JOBS.txt.
set -euo pipefail
ST=/encrypted/e3008/hohndor-dogohla-20261020/reads-repull-20261114
CF=/encrypted/e3008/hohndor-dogohla-20261020/mosaic-allgenes-20261114/inputs/confounded_v1
cd $ST
[ -s FREEZE_TEST.json ] || { echo "no FREEZE_TEST.json"; exit 1; }
[ ! -e work/TEST_JOBS.txt ] || { echo "already submitted"; exit 1; }
for G in cand40 need26; do
  L=work/treal_$G.txt
  awk '{print "T-real3", $1}' $L > work/test_sel_$G.txt
  awk '{print "T-real3", $1, "A"}' $L > work/test_A_$G.txt
  awk '{print "T-real3", $1, "B"}' $L > work/test_B_$G.txt
  awk -v cf=$CF '{print "T-real3", $1, "C", "test3", "--confounded-dir", cf}' $L > work/test_C_$G.txt
  awk -v cf=$CF '{print "T-real3", $1, "D", "test3", "--confounded-dir", cf}' $L > work/test_D_$G.txt
done
n40=$(($(wc -l < work/treal_cand40.txt)-1)); n26=$(($(wc -l < work/treal_need26.txt)-1))
w=$(sbatch --parsable --job-name=s2b-wgsx --array=0-$n26%12 src/wgsx.sbatch work/treal_need26.txt wgsx)
echo "wgsx need26 $w" >> work/TEST_JOBS.txt
for G in cand40 need26; do
  n=$([ $G = cand40 ] && echo $n40 || echo $n26)
  dep=""; [ $G = need26 ] && dep="--dependency=aftercorr:$w"
  s=$(sbatch --parsable $dep --job-name=s2b-sel --array=0-$n%12 --mem=16G --time=04:00:00 src/job.sbatch select_rule2.py work/test_sel_$G.txt)
  a=$(sbatch --parsable --dependency=aftercorr:$s --job-name=s2b-t1kA --array=0-$n%12 --cpus-per-task=4 --mem=32G --time=06:00:00 src/job.sbatch t1k_arm.py work/test_A_$G.txt)
  b=$(sbatch --parsable --dependency=aftercorr:$s --job-name=s2b-t1kB --array=0-$n%12 --cpus-per-task=4 --mem=32G --time=06:00:00 src/job.sbatch t1k_arm.py work/test_B_$G.txt)
  c=$(sbatch --parsable --dependency=aftercorr:$a --job-name=s2b-agC --array=0-$n%12 --cpus-per-task=4 --mem=32G --time=01:00:00 src/job.sbatch ag3_run.py work/test_C_$G.txt)
  d=$(sbatch --parsable --dependency=aftercorr:$b --job-name=s2b-agD --array=0-$n%12 --cpus-per-task=4 --mem=32G --time=01:00:00 src/job.sbatch ag3_run.py work/test_D_$G.txt)
  echo "sel $G $s" >> work/TEST_JOBS.txt; echo "t1kA $G $a" >> work/TEST_JOBS.txt; echo "t1kB $G $b" >> work/TEST_JOBS.txt
  echo "agC $G $c" >> work/TEST_JOBS.txt; echo "agD $G $d" >> work/TEST_JOBS.txt
done
cat work/TEST_JOBS.txt

#!/bin/bash
# Pre-registered retry rule (DESIGN.md): a task that timed out (task- or gene-level failure_kind=timeout,
# or no COMPLETE/FAILURE record at all, i.e. killed by Slurm) is resubmitted ONCE on a different node
# (--exclude=<first node>). Its first attempt is moved to work/<TAG>_attempt1/. A second timeout is not
# retried; its fallback (pass-1 calls) is scored and counted as a failure. Login-node safe (bash/grep/sbatch).
# usage: resubmit_timeouts.sh TAG TASKFILE ARRAYJOBID JOBNAME
set -u
N=/encrypted/e3008/hohndor-dogohla-20261020/mosaic-allgenes-20261114
TAG=$1; TASKS=$2; JOB=$3; JN=$4
cd $N
idx=0
while read -r rs d arm rest; do
  w=work/$TAG/$rs/$d/$arm
  a=work/${TAG}_attempt1/$rs/$d/$arm
  st=ok
  if [ ! -e $w/COMPLETE.json ]; then
    if [ -e $w/FAILURE.json ] && ! grep -q '"failure_kind": "timeout"' $w/FAILURE.json; then st=failed_other; else st=timeout; fi
  elif grep -q '"failure_kind": "timeout"' $w/*/FAILURE.json 2>/dev/null; then st=timeout; fi
  if [ $st = failed_other ]; then echo "FAILED (not a timeout, not retried): $w"; fi
  if [ $st = timeout ]; then
    if [ -e $a ]; then
      echo "SECOND TIMEOUT (fallback scored, counted as failure): $w"
    else
      node=$(sacct -j ${JOB}_$idx -X -n -o NodeList | tr -d ' ' | tail -1)
      mkdir -p $(dirname $a)
      [ -e $w ] && mv $w $a
      f=work/retry_${TAG}_${rs}_${d}_${arm}.txt
      sed -n "$((idx+1))p" $TASKS > $f
      echo "RETRY $w (first node $node): $(sbatch --parsable --job-name=${JN}-retry --exclude=$node --time=01:00:00 --array=0-0 src/s4.sbatch ag_run.py $f)"
    fi
  fi
  idx=$((idx+1))
done < $TASKS
echo "checked $idx tasks"

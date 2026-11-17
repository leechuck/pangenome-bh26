#!/bin/bash
# Step 2b pre-registered retry rule: a task that timed out (Slurm TIMEOUT, or no COMPLETE/FAILURE record, or for
# the builder a task-/gene-level failure_kind=timeout) is resubmitted ONCE on a different node (--exclude=<first node>);
# its first attempt is moved to <dir>_attempt1. A second timeout is not retried: the fallback is scored and counted as
# a failure. Non-timeout failures are not retried. Login-node safe (bash/grep/sacct/sbatch only).
# usage: resub3.sh KIND TASKFILE ARRAYJOBID JOBNAME     KIND: t1k (lines "RS DONOR ARM") | ag (lines "RS DONOR ARM TAG ...")
set -u
ST=/encrypted/e3008/hohndor-dogohla-20261020/reads-repull-20261114
KIND=$1; TASKS=$2; JOB=$3; JN=$4
cd $ST
idx=0
while read -r rs d arm tag rest; do
  if [ $KIND = t1k ]; then w=work/$rs/$d/$arm; script=t1k_arm.py; tl=06:00:00; extra="--cpus-per-task=4 --mem=32G"
  else w=ag/work/$tag/$rs/$d/$arm; script=ag3_run.py; tl=01:00:00; extra="--cpus-per-task=4 --mem=32G"; fi
  a=${w}_attempt1
  state=$(sacct -j ${JOB}_$idx -X -n -o State | tr -d ' ' | tail -1)
  st=ok
  if [ ! -e $w/COMPLETE.json ]; then
    if [ -e $w/FAILURE.json ] && ! grep -q '"failure_kind": "timeout"' $w/FAILURE.json && [ "$state" != TIMEOUT ]; then st=failed_other; else st=timeout; fi
  elif [ $KIND = ag ] && grep -q '"failure_kind": "timeout"' $w/*/FAILURE.json 2>/dev/null; then st=timeout; fi
  [ $st = failed_other ] && echo "FAILED (not a timeout, not retried): $w ($state)"
  if [ $st = timeout ]; then
    if [ -e $a ]; then echo "SECOND TIMEOUT (fallback, counted as failure): $w"
    else
      node=$(sacct -j ${JOB}_$idx -X -n -o NodeList | tr -d ' ' | tail -1)
      [ -e $w ] && mv $w $a
      f=work/retry_${KIND}_${rs}_${d}_${arm}.txt
      sed -n "$((idx+1))p" $TASKS > $f
      echo "RETRY $w (first node $node, state $state): $(sbatch --parsable --job-name=${JN}-retry --exclude=$node --time=$tl $extra --array=0-0 src/job.sbatch $script $f)"
    fi
  fi
  idx=$((idx+1))
done < $TASKS
echo "checked $idx tasks ($KIND $TASKS)"

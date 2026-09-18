#!/bin/bash
set -euo pipefail
base=/home/leechuck/hla/dogohla-series20260918
names=(full hprc asian_matched)
arm=${names[$SLURM_ARRAY_TASK_ID]}
export DOGOHLA_PANEL="$base/panels/$arm"
export SPECHLA_PG_ROOT="$base/arms/$arm"
if [[ ! -s "$SPECHLA_PG_ROOT/phase/PREPARED.json" ]]; then
python3 "$base/code/build_refs.py" --root "$SPECHLA_PG_ROOT" --folds 0 --threads 8
python3 "$base/code/experiment.py" --root "$SPECHLA_PG_ROOT" --release "$SPECHLA_PG_ROOT/phase" prepare --folds 0
fi
bash "$base/code/build_graph_pggb.sh" 0 DRB1 8

#!/bin/bash
# Transfer measurement outputs only, excluding live collate intermediates.
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p source/read_counts/development source/baseline/development
rsync -az --include='*.bin' --include='*.groups.tsv' --include='*.depth.tsv.gz' --exclude='*' \
 ddbj:/home/leechuck/hla/codex-targeted/counts/development/ source/read_counts/development/
rsync -az --include='*/' --include='C4Investigator_c4_summary.csv' --include='C4Investigator_c4_detailed.csv' --exclude='*' \
 ddbj:/home/leechuck/hla/codex-targeted/baseline/development/ source/baseline/development/

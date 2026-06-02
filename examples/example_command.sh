#!/usr/bin/env bash
set -euo pipefail

: "${BOOM_USERNAME:?Set BOOM_USERNAME before running this script}"
: "${BOOM_PASSWORD:?Set BOOM_PASSWORD before running this script}"

python scripts/boom_pipeline.py \
  --output_dir boom_filter_outputs \
  --filter_file configs/schecks.json \
  --start_time 2026-03-01T12:00:00 \
  --end_time 2026-05-13T12:00:00

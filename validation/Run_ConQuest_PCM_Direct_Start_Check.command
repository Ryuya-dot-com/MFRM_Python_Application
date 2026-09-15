#!/bin/zsh
set -eu
TASK_VALIDATION_DIR="${0:A:h}"
python3 "$TASK_VALIDATION_DIR/mml_pcm_conquest_direct_start.py" --directory "$TASK_VALIDATION_DIR/generated/mml_pcm_conquest_direct_start_20260914"

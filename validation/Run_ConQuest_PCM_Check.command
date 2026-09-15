#!/bin/zsh
set -eu
TASK_VALIDATION_DIR="${0:A:h}"
python3 "$TASK_VALIDATION_DIR/run_pcm_conquest_check.py" --output-dir "$TASK_VALIDATION_DIR/generated/mml_pcm_conquest_terminal_$(date +%Y%m%d_%H%M%S)"

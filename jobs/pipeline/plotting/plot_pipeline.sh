#!/bin/bash
# plot_pipeline.sh — turn an analysis directory's CSVs into figures.
#
# The plotting stage of the pipeline, separate from the analysis stage so that
# figures can be remade — restyled, re-unit'd, re-dated — without recomputing
# anything. The analysis scripts write CSVs and nothing else; every figure comes
# from here.
#
# Usage:
#   ./plot_pipeline.sh                 # plot the current directory
#   ./plot_pipeline.sh <dir>           # plot that analysis directory
#   ./plot_pipeline.sh <dir> --config <file>    # force a specific config
#
# Called by:
#   - by hand, in any directory holding <date>_*.csv from the analysis scripts
#   - distribution_run.sh / distribution_submit.slurm, when RUN_PLOTS=1
#   - submit_pipeline.sh / submit_pipeline_local.sh, via --run-plots 1
#
# CONFIG RESOLUTION. The first of these that exists wins, and the one used is
# printed so it is never a guess:
#   1. --config <file>, if passed
#   2. <target dir>/plot_pipeline.conf   — this run's own settings, so a
#      directory can carry the way its figures should look
#   3. <this dir>/plot_pipeline.conf     — the tracked default in the util
#      checkout, which is what most runs use
# Copy the default into an analysis directory and edit it there when one run
# needs to differ; nothing needs to be edited in tracked code to do that.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON="${PLOT_PYTHON:-python}"

usage() {
  sed -n '2,28p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
  exit "${1:-0}"
}

TARGET_DIR=""
CONFIG_FILE=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help)   usage 0 ;;
    --config)    CONFIG_FILE="$2"; shift 2 ;;
    --style)     PLOT_STYLE_OVERRIDE="$2"; shift 2 ;;
    --date)      PLOT_DATE_OVERRIDE="$2"; shift 2 ;;
    -*)          echo "Error: unknown option $1" >&2; usage 1 ;;
    *)           TARGET_DIR="$1"; shift ;;
  esac
done

TARGET_DIR="${TARGET_DIR:-$PWD}"
if [[ ! -d "$TARGET_DIR" ]]; then
  echo "Error: no such directory: $TARGET_DIR" >&2
  exit 1
fi
TARGET_DIR="$(cd "$TARGET_DIR" && pwd)"

# ---- config resolution, in the order documented above -----------------------
if [[ -n "$CONFIG_FILE" ]]; then
  if [[ ! -f "$CONFIG_FILE" ]]; then
    echo "Error: --config file not found: $CONFIG_FILE" >&2
    exit 1
  fi
  CONFIG_SOURCE="--config"
elif [[ -f "$TARGET_DIR/plot_pipeline.conf" ]]; then
  CONFIG_FILE="$TARGET_DIR/plot_pipeline.conf"
  CONFIG_SOURCE="run directory"
elif [[ -f "$SCRIPT_DIR/plot_pipeline.conf" ]]; then
  CONFIG_FILE="$SCRIPT_DIR/plot_pipeline.conf"
  CONFIG_SOURCE="pipeline default"
else
  echo "Error: no plot_pipeline.conf found." >&2
  echo "  Looked in: $TARGET_DIR/ (this run) and $SCRIPT_DIR/ (the default)." >&2
  echo "  The default ships with the repo; if it is missing, restore it from git." >&2
  exit 1
fi

echo "Config: $CONFIG_FILE  ($CONFIG_SOURCE)"

# Pre-declare every key the .py reads, so a conf predating any of them still
# works under `set -u`. The conf overrides whichever of these it sets, and an
# empty value leaves the .py's own default in effect.
PLOT_STYLE=""
PLOT_DPI=""
PLOT_FORMAT=""
PLOT_DATE=""
PLOT_RDF=""
PLOT_BAD=""
PLOT_DSF=""
PLOT_VDOS=""
PLOT_MSD=""
PLOT_VDOS_DYNMAT=""
PLOT_COMPOSITES=""
PLOT_FREQ_UNIT=""
PLOT_HEATMAP_CMAP=""
PLOT_PYTHON=""

# shellcheck source=/dev/null
source "$CONFIG_FILE"

# Flags win over the config, the same precedence the analysis pipeline uses.
PLOT_STYLE="${PLOT_STYLE_OVERRIDE:-$PLOT_STYLE}"
PLOT_DATE="${PLOT_DATE_OVERRIDE:-$PLOT_DATE}"
PYTHON="${PLOT_PYTHON:-$PYTHON}"

export PLOT_TARGET_DIR="$TARGET_DIR"
for _var in PLOT_STYLE PLOT_DPI PLOT_FORMAT PLOT_DATE \
            PLOT_RDF PLOT_BAD PLOT_DSF PLOT_VDOS PLOT_MSD PLOT_VDOS_DYNMAT \
            PLOT_COMPOSITES PLOT_FREQ_UNIT PLOT_HEATMAP_CMAP; do
  if [[ -n "${!_var}" ]]; then
    export "$_var"
  fi
done
unset _var

exec "$PYTHON" "$SCRIPT_DIR/plot_distributions.py"

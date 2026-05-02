#!/usr/bin/env bash
#
# step9_noninteractive_example.sh — non-interactive wrapper for step9_analyze_results.sh
#
# Edit the CONFIG section below, then run:
#   chmod +x step9_noninteractive_example.sh
#   ./step9_noninteractive_example.sh
#
# Or override from the shell without editing:
#   RESULTS_DIR=/path/to/out PRESET=quick ./step9_noninteractive_example.sh
#   # Same as: STEP9_RESULTS_DIR=/path/to/out PRESET=quick ./step9_noninteractive_example.sh
#

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# =============================================================================
# CONFIG — change these values (or export same names before running)
# =============================================================================

# Path to step9_analyze_results.sh (same folder as this script by default)
STEP9_SCRIPT="${STEP9_SCRIPT:-$HERE/step9_analyze_results.sh}"

# SAIGE-GENE results folder (chromosome *combined_results*.txt, etc.)
# RESULTS_DIR on the command line is honored; else STEP9_RESULTS_DIR; else this script's directory.
RESULTS_DIR="${RESULTS_DIR:-${STEP9_RESULTS_DIR:-$HERE}}"

# Primary operation: ONE name, OR several joined with "+" (no spaces).
# Default preset runs the usual pipeline; override OPERATIONS to pick anything below.
USE_PRESET="${USE_PRESET:-1}"

# If USE_PRESET=1, choose one of: standard | quick | full | merge_only | plots_only | custom_chain
PRESET="${PRESET:-standard}"

# When PRESET=custom_chain, set OPERATIONS explicitly, e.g. "mergeall+listgroups+findsig"
OPERATIONS="${OPERATIONS:-}"

# Optional positional arguments for step9 (arg 2–8). Use empty string '' to skip.
GENE_COORD="${GENE_COORD:-}"       # e.g. /path/to/genes.tsv (chr, start, end, gene)
POSITION_MODE="${POSITION_MODE:-midpoint}"   # start | end | midpoint
FILTER_GROUP="${FILTER_GROUP:-}"   # e.g. lof  or  lof;missense
FILTER_MAF="${FILTER_MAF:-}"       # e.g. 0.01
COORD_SOURCE="${COORD_SOURCE:-sequence}"     # sequence | ensembl
ENSEMBL_BUILD="${ENSEMBL_BUILD:-38}"        # 37 | 38 (only if COORD_SOURCE=ensembl)
ENSEMBL_RELEASE="${ENSEMBL_RELEASE:-}"      # optional e.g. 115, or full REST URL

# Common STEP9_* overrides (uncomment or export before running)
# export STEP9_PLOT_UNFILTERED=1
# export STEP9_MANHATTAN_LABEL_TOP_N=10
# export STEP9_MANHATTAN_FDR_ALPHA=0.1
# export STEP9_BONFERRONI_MAF_TESTS=3
# export STEP9_ENSEMBL_BATCH_SIZE=500
# export STEP9_ENSEMBL_PARALLEL=4

# =============================================================================
# Reference — every operation (copy into OPERATIONS or custom_chain)
# =============================================================================
#
# Data combination:
#   mergechrom       Merge chunked results by chromosome
#   mergeall         Combine all chromosomes into all_results.txt
#   listgroups       List annotation groups
#
# Significance:
#   findsig          Multiple significance files (GWS, suggestive, nominal, …)
#   findgws          p < 5e-8
#   findsug          p < 1e-5
#   findnom          p < 0.05
#
# Gene ranking:
#   top10 | top50 | top100
#
# Summary:
#   chromsum         Per-chromosome summary
#   groupsum         Per–annotation-group summary
#   fullsum          Full text summary report
#
# Plot data / PNG:
#   qqdata           qq_plot_data.txt
#   mandata          manhattan_plot_data.txt
#   plotdata         qq + mandata + makeplots (see step9)
#   makeplots        PNGs + group_statistics.txt (needs Rscript)
#
# Bundled:
#   standard         mergeall, listgroups, findsig, top50, groupsum, makeplots, fullsum
#   quick              mergeall, listgroups, top50, makeplots, fullsum
#   full               mergechrom, mergeall, … long pipeline + makeplots + fullsum
#
# Chaining examples:
#   mergeall+listgroups+findsig
#   mergeall+makeplots+fullsum
#   mergeall+qqdata+mandata
#

# =============================================================================
# Resolve OPERATIONS from preset (explicit OPERATIONS wins over PRESET)
# =============================================================================

if [ -z "${OPERATIONS}" ]; then
  if [ "${USE_PRESET}" != "1" ]; then
    echo "ERROR: Set OPERATIONS=... or USE_PRESET=1 with PRESET=standard|quick|..." >&2
    exit 2
  fi
  case "${PRESET}" in
    standard)
      OPERATIONS="standard"
      ;;
    quick)
      OPERATIONS="quick"
      ;;
    full)
      OPERATIONS="full"
      ;;
    merge_only)
      OPERATIONS="mergeall+listgroups"
      ;;
    plots_only)
      OPERATIONS="makeplots"
      ;;
    custom_chain)
      echo "ERROR: PRESET=custom_chain requires OPERATIONS to be set in CONFIG." >&2
      exit 2
      ;;
    *)
      echo "ERROR: Unknown PRESET=${PRESET}" >&2
      exit 2
      ;;
  esac
fi

if [ ! -f "${STEP9_SCRIPT}" ]; then
  echo "ERROR: step9 script not found: ${STEP9_SCRIPT}" >&2
  exit 2
fi

echo "=========================================="
echo "step9 non-interactive example"
echo "=========================================="
echo "STEP9_SCRIPT:  ${STEP9_SCRIPT}"
echo "RESULTS_DIR:   ${RESULTS_DIR}"
echo "OPERATIONS:    ${OPERATIONS}"
echo "POSITION_MODE: ${POSITION_MODE}"
echo "COORD_SOURCE:  ${COORD_SOURCE}"
echo ""

if [ "${VERBOSE:-0}" = "1" ]; then
  set -x
fi

exec bash "${STEP9_SCRIPT}" \
  -d "${RESULTS_DIR}" \
  "${OPERATIONS}" \
  "${GENE_COORD}" \
  "${POSITION_MODE}" \
  "${FILTER_GROUP}" \
  "${FILTER_MAF}" \
  "${COORD_SOURCE}" \
  "${ENSEMBL_BUILD}" \
  "${ENSEMBL_RELEASE}"

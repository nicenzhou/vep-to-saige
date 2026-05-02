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
# When all_results.txt includes Pvalue_Burden / Pvalue_SKAT (SAIGE-GENE), significance summaries
# (findsig, fullsum, chromsum, groupsum, findgws/findsug/findnom) report all present P columns.
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
# export STEP9_MANHATTAN_LABEL_EXTRA=PCSK9,LDLR,APOB
# export STEP9_MANHATTAN_LABEL_EXTRA_FILE=/path/to/gene_list.tsv
# QQ/Manhattan / qqdata/mandata P-value column (when merged results include Pvalue_Burden / Pvalue_SKAT):
#   STEP9_MANHATTAN_P_MODE=pvalue | burden | skat | all
#     pvalue — combined Pvalue only (default if unset)
#     burden — Pvalue_Burden column only
#     skat   — Pvalue_SKAT column only
#     all    — emit qq_plot_data*.txt / manhattan_plot_data*.txt and PNGs for each present column
# Interactive runs prompt when Burden/SKAT exist; batch runs use this env var.
# STEP9_MANHATTAN_P_MODE="${STEP9_MANHATTAN_P_MODE:-all}"

# Cauchy rows (all_results.txt is never modified; .all_results_no_cauchy.txt is built when needed):
#   STEP9_CAUCHY_MODE=off|plots|pipeline|full
#     off       — keep Cauchy everywhere
#     plots     — omit from QQ/Manhattan/Bonferroni gene count only
#     pipeline  — also omit from findsig/top/chromsum/groupsum/listgroups
#     full      — also omit Cauchy from full summary report totals
# export STEP9_CAUCHY_MODE=pipeline
# export STEP9_EXCLUDE_CAUCHY=1   # legacy: same as plots if STEP9_CAUCHY_MODE is unset
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
#   findsig          Tier files + console counts for Pvalue, Pvalue_Burden, Pvalue_SKAT (when columns exist)
#   findgws          p < 5e-8 (combined + optional Burden/SKAT counts)
#   findsug          p < 1e-5 (combined + optional Burden/SKAT counts)
#   findnom          p < 0.05 (combined + optional Burden/SKAT counts)
#
# Gene ranking:
#   top10 | top50 | top100
#
# Summary:
#   chromsum         Per-chromosome summary (sections per P column when Burden/SKAT exist)
#   groupsum         Per–annotation-group summary (sections per P column when Burden/SKAT exist)
#   fullsum          Full text summary report (threshold tables + group tables per P when present)
#
# Plot data / PNG:
#   qqdata           qq_plot_data.txt (optional qq_plot_data_burden.txt / *_skat.txt if STEP9_MANHATTAN_P_MODE=all)
#   mandata          manhattan_plot_data.txt (optional *_burden / *_skat per STEP9_MANHATTAN_P_MODE)
#   plotdata         qq + mandata + makeplots (see step9)
#   makeplots        PNGs + group_statistics.txt (needs Rscript; STEP9_MANHATTAN_P_MODE for multi-P)
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
echo "STEP9_CAUCHY_MODE: ${STEP9_CAUCHY_MODE:-off}  (export STEP9_CAUCHY_MODE or legacy STEP9_EXCLUDE_CAUCHY=1)"
echo "STEP9_MANHATTAN_P_MODE: ${STEP9_MANHATTAN_P_MODE:-<unset → pvalue in step9>}"
echo ""

if [ "${VERBOSE:-0}" = "1" ]; then
  set -x
fi

# Pass through to step9 only when set (e.g. STEP9_MANHATTAN_P_MODE=all)
[ -n "${STEP9_MANHATTAN_P_MODE:-}" ] && export STEP9_MANHATTAN_P_MODE

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

#!/usr/bin/env bash
# Run step9_analyze_results.sh for each cohort subfolder that contains
# chr*_combined_results.txt, then copy plots (*.png) and main tables into SummaryResults/
# with filenames prefixed by the folder name so cohorts do not overwrite each other.
#
# Quick usage (see --help):
#   ./step10_run_step9_all_folders.sh                         # all cohort folders → step9 standard
#   ./step10_run_step9_all_folders.sh AFR EUR -- quick        # selected folders → step9 quick
#   ./step10_run_step9_all_folders.sh -i                      # interactive prompts
#
# After "--", remaining arguments are passed to step9 for every cohort (do not pass --results-dir).
# Same env vars as step9; coordinate cache: STEP9_COORD_CACHE_DIR, reuse via STEP9_REUSE_* (auto).

set -euo pipefail

step10_show_help() {
    cat <<'EOF'
step10_run_step9_all_folders.sh — Step 10: batch Step 9 across multiple SAIGE result directories

WHAT IT DOES
  1) Finds cohort folders under the cohort root (each must contain chr*_combined_results.txt).
  2) Runs step9_analyze_results.sh once per folder (sets STEP9_RESULTS_DIR automatically).
  3) Copies plots + key tables into SummaryResults/ with names like COHORT_qq_plot.png.

USAGE
  step10_run_step9_all_folders.sh [OPTIONS] [COHORT_DIR_NAME ...] [-- STEP9_ARGUMENTS...]

OPTIONS
  -i, --interactive   Prompt for operations, coords, filters, Cauchy, plots (see step9).
                      With 2+ cohorts, asks whether to use one shared parameter set.
  -h, --help          Show this help.

COHORT_DIR_NAME
  Optional; each name must be a direct child folder of the cohort root (e.g. AFR EUR).
  If omitted, every qualifying subfolder of the cohort root is used (except SummaryResults).

STEP9_ARGUMENTS (after --)
  Same order as step9_analyze_results.sh. Omit -- entirely to run preset "standard".
  Do not pass --results-dir here (this script sets results dir per cohort).

COHORT ROOT (where AFR/, EUR/, ... live)
  • Script next to step9: cohort root = script dir, unless the script is under codes/ → parent dir.
  • Override: export STEP9_COHORT_ROOT=/path/to/parent

OUTPUT COPIES
  Default: STEP9_SUMMARY_DIR=<cohort_root>/SummaryResults

COORDINATE CACHE (gene file or Ensembl Manhattan)
  First cohort builds .gene_coords_lookup.txt; later cohorts reuse it when step9 args match.
  Cache: <cohort_root>/.step9_batch_coord_cache/<hash>/  (override STEP9_COORD_CACHE_DIR)

EXAMPLES
  ./step10_run_step9_all_folders.sh
  ./step10_run_step9_all_folders.sh Pooled -- quick
  STEP9_COHORT_ROOT=$PWD ./step10_run_step9_all_folders.sh AFR AMR -- standard
  nohup ./step10_run_step9_all_folders.sh > step10_batch.log 2>&1 &
EOF
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STEP9="${SCRIPT_DIR}/step9_analyze_results.sh"

# Where cohort folders (AFR/, EUR/, …) live: override with STEP9_COHORT_ROOT.
# Default: same directory as this script; if the script sits in codes/, use the parent directory (repo root).
if [[ -n "${STEP9_COHORT_ROOT:-}" ]]; then
    COHORT_ROOT="$(cd "${STEP9_COHORT_ROOT}" && pwd)"
elif [[ "$(basename "$SCRIPT_DIR")" == "codes" ]] && [[ -d "${SCRIPT_DIR}/.." ]]; then
    COHORT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
else
    COHORT_ROOT="$SCRIPT_DIR"
fi

SUMMARY_DIR="${STEP9_SUMMARY_DIR:-${COHORT_ROOT}/SummaryResults}"

INTERACTIVE_BATCH=false
USER_COHORT_TAGS=()
# Arguments forwarded to step9 for each cohort (non-interactive). Empty => "standard" only.
STEP9_PASS_ARGS=()

while [[ "$#" -gt 0 ]]; do
    case "$1" in
        -h|--help)
            step10_show_help
            exit 0
            ;;
        -i|--interactive)
            INTERACTIVE_BATCH=true
            shift
            ;;
        --)
            shift
            STEP9_PASS_ARGS+=("$@")
            break
            ;;
        *)
            USER_COHORT_TAGS+=("$1")
            shift
            ;;
    esac
done

if [[ "$INTERACTIVE_BATCH" == true ]] && [[ ${#STEP9_PASS_ARGS[@]} -gt 0 ]]; then
    echo "ERROR: use either --interactive (-i) or '--' step9 arguments, not both." >&2
    exit 1
fi

if [[ ! -f "$STEP9" ]]; then
    echo "ERROR: step9_analyze_results.sh not found at $STEP9" >&2
    exit 1
fi

mkdir -p "$SUMMARY_DIR"
START_SEC=$(date +%s)

if [[ "$INTERACTIVE_BATCH" == true ]] && [[ ! -t 0 ]]; then
    echo "ERROR: --interactive requires a terminal (stdin). Run from Terminal/iTerm, not from a non-interactive context." >&2
    exit 1
fi

# Y/y/empty with default y -> true
batch_prompt_is_yes() {
    local raw="${1-}"
    local default="${2:-y}"
    local s
    s=$(printf '%s' "$raw" | tr '[:upper:]' '[:lower:]' | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')
    if [[ -z "$s" ]]; then
        [[ "$default" == "y" ]]
        return
    fi
    case "$s" in
        y|yes) return 0 ;;
        n|no) return 1 ;;
        *) [[ "$default" == "y" ]] ;;
    esac
}

# True if step9 will build .gene_coords_lookup.txt from a local file or Ensembl (not sequence-only).
batch_argv_uses_coord_resources() {
    local -a a=( "$@" )
    [[ ${#a[@]} -lt 2 ]] && return 1
    if [[ -n "${a[1]}" ]] && [[ -f "${a[1]}" ]]; then
        return 0
    fi
    if [[ ${#a[@]} -ge 6 ]] && [[ "${a[5]}" == "ensembl" ]]; then
        return 0
    fi
    return 1
}

batch_argv_canonical_for_hash() {
    local -a a=( "$@" )
    if [[ ${#a[@]} -ge 2 ]] && [[ -n "${a[1]}" ]] && [[ -f "${a[1]}" ]]; then
        local d b
        d="$(cd "$(dirname "${a[1]}")" && pwd)"
        b="$(basename "${a[1]}")"
        a[1]="${d}/${b}"
    fi
    printf '%s\n' "${a[@]}"
}

batch_coord_cache_key() {
    local h
    if command -v shasum >/dev/null 2>&1; then
        h=$(batch_argv_canonical_for_hash "$@" | shasum -a 256 2>/dev/null | awk '{print $1}')
    else
        h=$(batch_argv_canonical_for_hash "$@" | sha256sum 2>/dev/null | awk '{print $1}')
    fi
    [[ -n "$h" ]] || h=$(batch_argv_canonical_for_hash "$@" | cksum | awk '{print $1}')
    printf '%s' "$h"
}

batch_coord_cache_clear_env() {
    unset STEP9_REUSE_COORD_LOOKUP STEP9_REUSE_ENSEMBL_REPORT 2>/dev/null || true
}

batch_coord_cache_apply_env() {
    local slot="$1"
    export STEP9_REUSE_COORD_LOOKUP="${slot}/gene_coords_lookup.txt"
    if [[ -f "${slot}/ensembl_lookup_report.txt" ]]; then
        export STEP9_REUSE_ENSEMBL_REPORT="${slot}/ensembl_lookup_report.txt"
    else
        unset STEP9_REUSE_ENSEMBL_REPORT 2>/dev/null || true
    fi
}

batch_coord_cache_save() {
    local src_dir="$1"
    local slot="$2"
    [[ -f "${src_dir}/.gene_coords_lookup.txt" ]] || return 0
    mkdir -p "$slot"
    cp -f "${src_dir}/.gene_coords_lookup.txt" "${slot}/gene_coords_lookup.txt"
    if [[ -f "${src_dir}/.ensembl_lookup_report.txt" ]]; then
        cp -f "${src_dir}/.ensembl_lookup_report.txt" "${slot}/ensembl_lookup_report.txt"
    fi
    echo ""
    echo " -> Saved coordinate lookup cache for later cohorts: ${slot}"
}

# Prompt for step9 options; sets BATCH_ARGV and BATCH_CAUCHY, BATCH_PLOT, BATCH_PMODE.
# step9 order: OPERATIONS [gene_coord] [position] [filter_group] [filter_maf] [coord_source] [ensembl_build] [ensembl_release]
step9_interactive_collect_params() {
    local title="$1"
    echo ""
    echo "────────────────────────────────────────"
    echo " Step 9 parameters: ${title}"
    echo "────────────────────────────────────────"
    local ops gc pm fg fm cs eb er pplot cauchy pmode
    read -r -p "  Operations [standard] (e.g. standard, quick, full, mergeall+findsig+top50): " ops
    ops="${ops:-standard}"
    read -r -p "  Gene coordinate file (tab: chr start end gene) [Enter = none]: " gc
    read -r -p "  Position mode (start|end|midpoint) [midpoint]: " pm
    pm="${pm:-midpoint}"
    read -r -p "  Filter annotation group(s) [Enter = none]: " fg
    read -r -p "  Filter max_MAF ≤ threshold [Enter = none]: " fm
    if [[ -z "$gc" ]]; then
        read -r -p "  Coordinate source when no file (sequence|ensembl) [sequence]: " cs
        cs="${cs:-sequence}"
        if [[ "${cs}" == "ensembl" ]]; then
            read -r -p "  Ensembl genome build (37|38) [38]: " eb
            eb="${eb:-38}"
            read -r -p "  Ensembl release or URL [Enter = skip]: " er
        else
            eb="38"
            er=""
        fi
    else
        cs="sequence"
        eb="38"
        er=""
    fi
    read -r -p "  Cauchy mode (off|plots|pipeline|full) [off]: " cauchy
    cauchy="${cauchy:-off}"
    read -r -p "  Generate QQ/Manhattan for all groups×MAF (1|0) [1]: " pplot
    pplot="${pplot:-1}"
    read -r -p "  Manhattan / QQ P column (pvalue|burden|skat|all) [pvalue]: " pmode
    pmode="${pmode:-pvalue}"

    BATCH_CAUCHY="$cauchy"
    BATCH_PLOT="$pplot"
    BATCH_PMODE="$pmode"

    local -a argv=( "$ops" )
    if [[ -n "$gc" ]]; then
        argv+=( "$gc" "$pm" )
        [[ -n "$fg" || -n "$fm" ]] && argv+=( "$fg" "$fm" )
    else
        if [[ "$cs" == "ensembl" ]]; then
            argv+=( "" "$pm" "${fg}" "${fm}" "$cs" "$eb" "${er:-}" )
        elif [[ -n "$fg" || -n "$fm" ]]; then
            argv+=( "" "$pm" "$fg" "$fm" )
        fi
    fi
    BATCH_ARGV=( "${argv[@]}" )
}

step9_run_batch_argv() {
    export STEP9_CAUCHY_MODE="$BATCH_CAUCHY"
    export STEP9_PLOT_UNFILTERED="$BATCH_PLOT"
    export STEP9_MANHATTAN_P_MODE="$BATCH_PMODE"
    printf '  Running: STEP9_RESULTS_DIR=%s bash step9' "${STEP9_RESULTS_DIR}"
    printf ' %q' "${BATCH_ARGV[@]}"
    printf '\n'
    bash "$STEP9" "${BATCH_ARGV[@]}"
}

step10_interactive_for_cohort() {
    local tag="$1"
    step9_interactive_collect_params "cohort ${tag}"
    step9_run_batch_argv
}

shopt -s nullglob

step9_copy_summary_artifacts() {
    local cohort_dir="$1"
    local tag="$2"
    local dest="$3"

    # PNG plots (QQ, Manhattan, optional burden/SKAT)
    local f
    for f in "${cohort_dir}"/*.png; do
        [[ -f "$f" ]] || continue
        cp -f "$f" "${dest}/${tag}_$(basename "$f")"
    done

    # Core tables from standard = mergeall + findsig + top50 + groupsum + makeplots + fullsum
    local base
    for base in \
        analysis_summary.txt \
        top50_genes.txt \
        annotation_group_summary.txt \
        group_statistics.txt \
        genome_wide_sig.txt \
        suggestive_sig.txt \
        nominal_sig.txt \
        sig_p001.txt \
        chromosome_summary.txt \
        annotation_groups.txt
    do
        if [[ -f "${cohort_dir}/${base}" ]]; then
            cp -f "${cohort_dir}/${base}" "${dest}/${tag}_${base}"
        fi
    done

    # Optional Pvalue_Burden / Pvalue_SKAT significance splits
    for f in \
        "${cohort_dir}/genome_wide_sig_burden.txt" \
        "${cohort_dir}/suggestive_sig_burden.txt" \
        "${cohort_dir}/nominal_sig_burden.txt" \
        "${cohort_dir}/sig_p001_burden.txt" \
        "${cohort_dir}/genome_wide_sig_skat.txt" \
        "${cohort_dir}/suggestive_sig_skat.txt" \
        "${cohort_dir}/nominal_sig_skat.txt" \
        "${cohort_dir}/sig_p001_skat.txt"
    do
        [[ -f "$f" ]] || continue
        cp -f "$f" "${dest}/${tag}_$(basename "$f")"
    done

    # QQ / Manhattan intermediate tables (optional multi-mode)
    for f in "${cohort_dir}"/qq_plot_data*.txt "${cohort_dir}"/manhattan_plot_data*.txt; do
        [[ -f "$f" ]] || continue
        cp -f "$f" "${dest}/${tag}_$(basename "$f")"
    done
}

ran_any=false

COHORT_DIRS=()
if [[ ${#USER_COHORT_TAGS[@]} -gt 0 ]]; then
    for tag in "${USER_COHORT_TAGS[@]}"; do
        cohort_dir="${COHORT_ROOT%/}/${tag}/"
        if [[ ! -d "$cohort_dir" ]]; then
            echo "ERROR: cohort folder not found: $cohort_dir" >&2
            exit 1
        fi
        if ! compgen -G "${cohort_dir}"chr*_combined_results.txt >/dev/null; then
            echo "ERROR: no chr*_combined_results.txt in ${cohort_dir}" >&2
            exit 1
        fi
        COHORT_DIRS+=("$cohort_dir")
    done
else
    for cohort_dir in "${COHORT_ROOT}"/*/; do
        [[ -d "$cohort_dir" ]] || continue
        tag="$(basename "${cohort_dir%/}")"
        if [[ "$tag" == "SummaryResults" ]]; then
            continue
        fi
        if ! compgen -G "${cohort_dir}"chr*_combined_results.txt >/dev/null; then
            continue
        fi
        COHORT_DIRS+=("$cohort_dir")
    done
fi

if [[ ${#COHORT_DIRS[@]} -eq 0 ]]; then
    echo "ERROR: No cohort folders with chr*_combined_results.txt under ${COHORT_ROOT}" >&2
    echo "  Tip: set STEP9_COHORT_ROOT to the directory that contains AFR/, EUR/, …" >&2
    exit 1
fi

_cohort_list_note="auto-discovered"
[[ ${#USER_COHORT_TAGS[@]} -gt 0 ]] && _cohort_list_note="$(printf '%s ' "${USER_COHORT_TAGS[@]}")"

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo " Step 10 — batch Step 9 (step10_run_step9_all_folders.sh)"
echo "   Cohort root:     ${COHORT_ROOT}"
echo "   Summary copies:  ${SUMMARY_DIR}"
echo "   Cohort folders:  ${#COHORT_DIRS[@]}  (${_cohort_list_note})"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

SAME_PARAMS_FOR_ALL=false
if [[ "$INTERACTIVE_BATCH" == true ]] && [[ ${#COHORT_DIRS[@]} -gt 0 ]]; then
    if [[ ${#COHORT_DIRS[@]} -eq 1 ]]; then
        SAME_PARAMS_FOR_ALL=true
    else
        echo ""
        read -r -p "Use the same step9 parameters for all ${#COHORT_DIRS[@]} cohort folders? [Y/n]: " _same_ans
        if batch_prompt_is_yes "${_same_ans}" "y"; then
            SAME_PARAMS_FOR_ALL=true
        fi
    fi
    if [[ "$SAME_PARAMS_FOR_ALL" == true ]]; then
        step9_interactive_collect_params "all ${#COHORT_DIRS[@]} cohort folders"
    fi
fi

# Same argv across cohorts → reuse one coordinate lookup (local file or Ensembl) for cohorts 2..N.
declare -a ARGV_REF=()
if [[ "$INTERACTIVE_BATCH" == true ]] && [[ "$SAME_PARAMS_FOR_ALL" == true ]]; then
    ARGV_REF=( "${BATCH_ARGV[@]}" )
elif [[ "$INTERACTIVE_BATCH" != true ]]; then
    if [[ ${#STEP9_PASS_ARGS[@]} -eq 0 ]]; then
        ARGV_REF=( standard )
    else
        ARGV_REF=( "${STEP9_PASS_ARGS[@]}" )
    fi
fi

USE_COORD_CACHE=false
CACHE_SLOT=""
if [[ ${#COHORT_DIRS[@]} -gt 1 ]] && { [[ "$INTERACTIVE_BATCH" != true ]] || [[ "$SAME_PARAMS_FOR_ALL" == true ]]; }; then
    if [[ ${#ARGV_REF[@]} -gt 0 ]] && batch_argv_uses_coord_resources "${ARGV_REF[@]}"; then
        COORD_CACHE_DIR="${STEP9_COORD_CACHE_DIR:-${COHORT_ROOT}/.step9_batch_coord_cache}"
        mkdir -p "$COORD_CACHE_DIR"
        CACHE_SLOT="${COORD_CACHE_DIR}/$(batch_coord_cache_key "${ARGV_REF[@]}")"
        USE_COORD_CACHE=true
        echo ""
        echo "Coordinate lookup cache (shared across cohorts): ${CACHE_SLOT}"
    fi
fi

cohort_idx=0
for cohort_dir in "${COHORT_DIRS[@]}"; do
    [[ -d "$cohort_dir" ]] || continue
    tag="$(basename "${cohort_dir%/}")"
    if [[ "$tag" == "SummaryResults" ]]; then
        continue
    fi
    if ! compgen -G "${cohort_dir}"chr*_combined_results.txt >/dev/null; then
        echo "SKIP (no chr*_combined_results.txt): ${cohort_dir}" >&2
        continue
    fi

    ran_any=true
    echo ""
    echo "=========================================="
    echo " Cohort folder: ${tag}"
    echo " Results dir:   ${cohort_dir}"
    echo "=========================================="

    cd "$SCRIPT_DIR"
    export STEP9_RESULTS_DIR="$cohort_dir"

    batch_coord_cache_clear_env
    if [[ "$USE_COORD_CACHE" == true ]] && [[ "$cohort_idx" -gt 0 ]] && [[ -n "$CACHE_SLOT" ]] && [[ -f "${CACHE_SLOT}/gene_coords_lookup.txt" ]]; then
        batch_coord_cache_apply_env "$CACHE_SLOT"
        echo " -> Reusing cached coordinate lookup (skips repeated Ensembl REST / local rebuild)."
    fi

    if [[ "$INTERACTIVE_BATCH" == true ]]; then
        if [[ "$SAME_PARAMS_FOR_ALL" == true ]]; then
            step9_run_batch_argv
        else
            step10_interactive_for_cohort "$tag"
        fi
    else
        export STEP9_PLOT_UNFILTERED="${STEP9_PLOT_UNFILTERED:-1}"
        if [[ ${#STEP9_PASS_ARGS[@]} -eq 0 ]]; then
            bash "$STEP9" standard
        else
            printf '  Running: STEP9_RESULTS_DIR=%s bash step9' "${STEP9_RESULTS_DIR}"
            printf ' %q' "${STEP9_PASS_ARGS[@]}"
            printf '\n'
            bash "$STEP9" "${STEP9_PASS_ARGS[@]}"
        fi
    fi

    echo " -> Copying plots and tables to ${SUMMARY_DIR} as ${tag}_*"
    step9_copy_summary_artifacts "$cohort_dir" "$tag" "$SUMMARY_DIR"

    if [[ "$USE_COORD_CACHE" == true ]] && [[ "$cohort_idx" -eq 0 ]] && [[ -n "$CACHE_SLOT" ]]; then
        batch_coord_cache_save "$cohort_dir" "$CACHE_SLOT"
    fi
    cohort_idx=$((cohort_idx + 1))
done

shopt -u nullglob

if [[ "$ran_any" != true ]]; then
    echo "ERROR: No cohort runs completed under ${COHORT_ROOT}" >&2
    exit 1
fi

_end=$(date +%s)
_elapsed=$((_end - START_SEC))
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
printf ' Done: %s cohort(s) in %ss\n' "${cohort_idx}" "${_elapsed}"
echo " Prefixed plots & tables: ${SUMMARY_DIR}"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

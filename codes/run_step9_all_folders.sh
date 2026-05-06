#!/usr/bin/env bash
# Run step9_analyze_results.sh for each cohort subfolder that contains
# chr*_combined_results.txt, then copy plots (*.png) and main tables into SummaryResults/
# with filenames prefixed by the folder name so cohorts do not overwrite each other.
#
# Usage:
#   ./run_step9_all_folders.sh                  # all cohort folders; fixed pipeline: standard
#   ./run_step9_all_folders.sh AFR EUR          # only these folders (same fixed pipeline)
#   ./run_step9_all_folders.sh -i               # prompt for operations/filters/env per cohort
#   ./run_step9_all_folders.sh --interactive AFR EUR
#
# Non-interactive exports (optional): STEP9_CAUCHY_MODE, STEP9_PLOT_UNFILTERED,
# STEP9_MANHATTAN_P_MODE (same for every cohort unless you use -i).

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STEP9="${SCRIPT_DIR}/step9_analyze_results.sh"
SUMMARY_DIR="${SCRIPT_DIR}/SummaryResults"

INTERACTIVE_BATCH=false
USER_COHORT_TAGS=()
while [[ "$#" -gt 0 ]]; do
    case "$1" in
        -i|--interactive)
            INTERACTIVE_BATCH=true
            shift
            ;;
        --)
            shift
            USER_COHORT_TAGS+=("$@")
            break
            ;;
        *)
            USER_COHORT_TAGS+=("$1")
            shift
            ;;
    esac
done

if [[ ! -f "$STEP9" ]]; then
    echo "ERROR: step9_analyze_results.sh not found at $STEP9" >&2
    exit 1
fi

mkdir -p "$SUMMARY_DIR"

if [[ "$INTERACTIVE_BATCH" == true ]] && [[ ! -t 0 ]]; then
    echo "ERROR: --interactive requires a terminal (stdin). Run from Terminal/iTerm, not from a non-interactive context." >&2
    exit 1
fi

# Prompt once per cohort; builds positional args for step9 + exports Cauchy/plot/P-mode env vars.
# step9 order: OPERATIONS [gene_coord] [position] [filter_group] [filter_maf] [coord_source] [ensembl_build] [ensembl_release]
run_step9_interactive_for_cohort() {
    local tag="$1"
    echo ""
    echo "────────────────────────────────────────"
    echo " Parameters for cohort: ${tag}"
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

    export STEP9_CAUCHY_MODE="$cauchy"
    export STEP9_PLOT_UNFILTERED="$pplot"
    export STEP9_MANHATTAN_P_MODE="$pmode"

    # Positional args must match step9's optional chain (see step9_analyze_results.sh).
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

    printf '  Running: STEP9_RESULTS_DIR=%s bash step9' "${STEP9_RESULTS_DIR}"
    printf ' %q' "${argv[@]}"
    printf '\n'
    bash "$STEP9" "${argv[@]}"
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
        cohort_dir="${SCRIPT_DIR%/}/${tag}/"
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
    for cohort_dir in "${SCRIPT_DIR}"/*/; do
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
    if [[ "$INTERACTIVE_BATCH" == true ]]; then
        run_step9_interactive_for_cohort "$tag"
    else
        # Honor caller env for every cohort (optional).
        export STEP9_PLOT_UNFILTERED="${STEP9_PLOT_UNFILTERED:-1}"
        bash "$STEP9" standard
    fi

    echo " -> Copying plots and tables to ${SUMMARY_DIR} as ${tag}_*"
    step9_copy_summary_artifacts "$cohort_dir" "$tag" "$SUMMARY_DIR"
done

shopt -u nullglob

if [[ "$ran_any" != true ]]; then
    echo "ERROR: No subfolders under ${SCRIPT_DIR} contain chr*_combined_results.txt" >&2
    exit 1
fi

echo ""
echo "Finished. Summary copies are under: ${SUMMARY_DIR}"

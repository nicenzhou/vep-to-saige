#!/usr/bin/env bash
# SAIGE group -> regenie (bash only; no awk). See Milind et al. medRxiv 2024.11.11.24317065v2
# for REGENIE burden settings (LoF MAF < 1%, etc.).
#
# Edit paths in the block below, then run:  bash saige_group_to_regenie_plain.sh
# Environment variables still override these defaults when set.
# If you see "set: pipefail" / "invalid option", the file has Windows (CRLF) line endings — run: dos2unix this_file
#
# ========== EDIT DEFAULT PATHS (optional) ==========
SAIGE_GROUP_DIR_DEFAULT=""
OUT_DIR_DEFAULT=""
# ====================================================

set -euo pipefail

SAIGE_GROUP_DIR="${SAIGE_GROUP_DIR:-${SAIGE_GROUP_DIR_DEFAULT:-./saige_groups}}"
OUT_DIR="${OUT_DIR:-${OUT_DIR_DEFAULT:-./regenie_masks_from_saige}}"
STRIP_CHR_PREFIX="${STRIP_CHR_PREFIX:-0}"

# --- bash helpers (SAIGE variant IDs like chr1:pos:ref:alt) ---
saige_norm_vid() {
  local s="$1"
  if [[ "${STRIP_CHR_PREFIX}" == "1" ]]; then
    case "${s}" in
      [Cc][Hh][Rr]*) s="${s:3}" ;;
    esac
  fi
  printf '%s' "${s}"
}

# Sets _saige_chr and _saige_pos (integer string) from variant ID, or returns 1 if unparseable
# (same rule as legacy awk: need at least four colon-separated fields chr:pos:ref:alt…)
saige_chr_pos_from_vid() {
  local v="$1" s
  local -a ap
  _saige_chr=""
  _saige_pos=""
  s="$(saige_norm_vid "${v}")"
  IFS=: read -r -a ap <<< "${s}" || true
  ((${#ap[@]} >= 4)) || return 1
  case "${ap[0]}" in
    [Cc][Hh][Rr]*) ap[0]="${ap[0]:3}" ;;
  esac
  _saige_chr="${ap[0]}"
  _saige_pos="${ap[1]}"
  return 0
}

# Read SAIGE group file; append regenie anno + setlist (tab-separated columns).
saige_to_regenie_files() {
  local in_path="$1" annof="$2" setf="$3"
  local want=0 skipw=0 lastgene="" gene="" nv=0
  local -a vars=()
  local -a annos=()

  : >"${annof}"
  : >"${setf}"

  while IFS= read -r line || [[ -n "${line}" ]]; do
    [[ -z "${line//[[:space:]]/}" ]] && continue
    read -r -a tok <<< "${line}"
    ((${#tok[@]} >= 2)) || continue

    if [[ "${tok[1]}" == "var" ]]; then
      skipw=0
      gene="${tok[0]}"
      vars=("${tok[@]:2}")
      nv=${#vars[@]}
      want=1
      continue
    fi

    if (( want == 1 )); then
      if [[ "${tok[1]}" == "anno" && "${tok[0]}" == "${gene}" ]]; then
        annos=("${tok[@]:2}")
        if ((${#annos[@]} != nv)); then
          echo "${in_path}: variant/annotation count mismatch for ${gene}" >&2
          want=0
          continue
        fi
        local pmin=-1 cmin="" vlist="" sep="" i vv
        for ((i = 0; i < nv; i++)); do
          vv="$(saige_norm_vid "${vars[i]}")"
          printf '%s\t%s\t%s\n' "${vv}" "${gene}" "${annos[i]}" >>"${annof}"
          if saige_chr_pos_from_vid "${vars[i]}"; then
            if [[ "${_saige_pos}" =~ ^[0-9]+$ ]]; then
              if (( pmin < 0 || _saige_pos < pmin )); then
                pmin="${_saige_pos}"
                cmin="${_saige_chr}"
              fi
            fi
          fi
          vlist+="${sep}${vv}"
          sep=","
        done
        if (( pmin >= 0 )); then
          printf '%s\t%s\t%s\t%s\n' "${gene}" "${cmin}" "${pmin}" "${vlist}" >>"${setf}"
        fi
        want=0
        skipw=1
        lastgene="${gene}"
        continue
      fi
      echo "${in_path}: expected anno line after var for ${gene}" >&2
      want=0
      continue
    fi

    if (( skipw == 1 )) && [[ "${tok[1]}" == "weight" && "${tok[0]}" == "${lastgene}" ]]; then
      continue
    fi
    if (( skipw == 1 )); then
      skipw=0
    fi
  done <"${in_path}"
}

mkdir -p "${OUT_DIR}"

for n in $(seq 1 22); do
  tag="$(printf '%02d' "${n}")"
  in_path=""
  if [[ -f "${SAIGE_GROUP_DIR}/chr${tag}_group.txt" ]]; then
    in_path="${SAIGE_GROUP_DIR}/chr${tag}_group.txt"
  elif [[ -f "${SAIGE_GROUP_DIR}/chr${n}_group.txt" ]]; then
    in_path="${SAIGE_GROUP_DIR}/chr${n}_group.txt"
  fi
  if [[ -z "${in_path}" ]]; then
    echo "warning: missing group file for chr ${tag}" >&2
    continue
  fi
  sub="${OUT_DIR}/chr${tag}"
  mkdir -p "${sub}"
  annof="${sub}/chr${tag}.regenie.annotations.txt"
  setf="${sub}/chr${tag}.regenie.setlist.txt"
  saige_to_regenie_files "${in_path}" "${annof}" "${setf}"
  echo "wrote ${annof} and ${setf}"
done

{
  printf '%s\n' "mask_lof	lof"
  printf '%s\n' "mask_missense	missense"
  printf '%s\n' "mask_synonymous	synonymous"
  printf '%s\n' "mask_lof_missense	lof,missense"
  printf '%s\n' "# Edit mask rows so the second column lists annotation tokens exactly as they appear in the SAIGE anno lines."
} >"${OUT_DIR}/mask_def.example.txt"
echo "wrote example mask definition -> ${OUT_DIR}/mask_def.example.txt"

echo "Done. Point regenie --anno-file / --set-list to per-chr *.txt files under ${OUT_DIR}/chrNN/."
echo "Example mask file: ${OUT_DIR}/mask_def.example.txt (merge with your real annotation spellings)."

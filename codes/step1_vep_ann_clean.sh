#!/bin/bash

# Script: step1_vep_ann_clean.sh
# Description: VEP annotation extraction with configurable LoF definition and deduplication
# Usage: vep_ann_clean.sh input.vcf.gz output.txt [threads] [lof_mode]

set -euo pipefail
trap '' PIPE

#==========================================
# ARGUMENTS
#==========================================

INPUT="${1:-}"
OUTPUT="${2:-}"
THREADS="${3:-1}"
LOF_MODE="${4:-vep_lof}"

#==========================================
# VALIDATION
#==========================================

if [ -z "$INPUT" ] || [ -z "$OUTPUT" ]; then
    cat << 'EOF'
Usage: step1_vep_ann_clean.sh <input.vcf.gz> <output.txt> [threads] [lof_mode]

Description:
  Extract VEP annotations with configurable LoF definition.
  Filters out LC LoF variants and NA group variants.
  Deduplicates identical annotations for the same variant.

Arguments:
  input.vcf.gz    Input VCF file with VEP annotations
  output.txt      Output annotation file
  threads         Number of parallel threads (default: 1)
  lof_mode        LoF grouping mode (default: vep_lof)

LoF Modes:
  loftee_only     Only LOFTEE HC variants → group='lof'
  vep_lof         VEP LoF consequences → group='lof' (Default)
  any             LOFTEE HC OR VEP LoF consequences → group='lof'

VEP LoF Consequences (vep_lof mode):
  - frameshift_variant
  - stop_gained
  - stop_lost
  - splice_acceptor_variant
  - splice_donor_variant
  - inframe_deletion
  - inframe_insertion

Deduplication:
  If same variant has identical rows (same gene, consequence, LoF, flags, group),
  only one row is kept to minimize file size.

Examples:
  ./step1_vep_ann_clean.sh input.vcf.gz output.txt 1
  ./step1_vep_ann_clean.sh input.vcf.gz output.txt 1 loftee_only
  ./step1_vep_ann_clean.sh input.vcf.gz output.txt 4 vep_lof
EOF
    exit 1
fi

[ ! -f "$INPUT" ] && { echo "Error: Input file '$INPUT' not found!" >&2; exit 1; }

if ! command -v bcftools &> /dev/null; then
    echo "Error: bcftools not found. Please install bcftools." >&2
    exit 1
fi

case "$LOF_MODE" in
    loftee_only|vep_lof|any)
        ;;
    *)
        echo "Error: Invalid lof_mode '$LOF_MODE'" >&2
        echo "Valid modes: loftee_only, vep_lof, any" >&2
        exit 1
        ;;
esac

#==========================================
# MAIN PROCESSING
#==========================================

echo "[$(date)] Starting VEP annotation extraction..." >&2
echo "  Input: $INPUT" >&2
echo "  Output: $OUTPUT" >&2
echo "  Threads: $THREADS" >&2
echo "  LoF mode: $LOF_MODE" >&2
echo "  Deduplication: ON" >&2

{
    echo -e "Variant\tGene\tConsequence\tLoF\tLoF_flags\tGroup\tTotal_Anno\tNON_CAN_SPLICE_Count"
    
    bcftools view "$INPUT" --exclude 'CSQ~"LoF=LC"' --threads "$THREADS" 2>/dev/null | \
      bcftools +split-vep -d \
        -f '%CHROM:%POS:%REF:%ALT\t%SYMBOL\t%Consequence\t%LoF\t%LoF_flags\n' 2>/dev/null | \
      awk -F'\t' -v lof_mode="$LOF_MODE" '
      BEGIN {
        prev_variant = ""
        n = 0
        nc = 0
        duplicates_removed = 0
      }
      
      {
        variant = $1
        symbol = $2
        consequence = $3
        lof = $4
        lof_flags = $5
        
        # When variant changes, output previous variant data
        if (variant != prev_variant && prev_variant != "") {
          for (i = 1; i <= n; i++) {
            if (groups[i] != "NA") {
              print lines[i] "\t" n "\t" nc
            }
          }
          delete lines
          delete groups
          delete row_keys
          n = 0
          nc = 0
        }
        
        # Determine functional group based on lof_mode
        group = "NA"
        
        is_vep_lof = (consequence ~ /frameshift_variant|stop_gained|stop_lost|splice_acceptor_variant|splice_donor_variant|inframe_deletion|inframe_insertion/)
        is_loftee_hc = (lof == "HC")
        
        if (lof_mode == "loftee_only") {
          if (is_loftee_hc) {
            group = "lof"
          }
        }
        else if (lof_mode == "vep_lof") {
          if (is_vep_lof) {
            group = "lof"
          }
        }
        else if (lof_mode == "any") {
          if (is_loftee_hc || is_vep_lof) {
            group = "lof"
          }
        }
        
        if (group == "NA" && consequence ~ /missense_variant/) {
          group = "missense"
        }
        
        if (group == "NA" && consequence ~ /synonymous_variant/) {
          group = "synonymous"
        }
        
        # Create unique key for deduplication
        row_key = variant SUBSEP symbol SUBSEP consequence SUBSEP lof SUBSEP lof_flags SUBSEP group
        
        # Check if this exact row has been seen for this variant
        if (row_key in row_keys) {
          duplicates_removed++
        }
        else {
          row_keys[row_key] = 1
          
          lines[++n] = variant "\t" symbol "\t" consequence "\t" lof "\t" lof_flags "\t" group
          groups[n] = group
          
          if (lof_flags ~ /NON_CAN_SPLICE/) {
            nc++
          }
        }
        
        prev_variant = variant
      }
      
      END {
        # Output last variant
        if (n > 0) {
          for (i = 1; i <= n; i++) {
            if (groups[i] != "NA") {
              print lines[i] "\t" n "\t" nc
            }
          }
        }
        
        if (duplicates_removed > 0) {
          print "  Duplicates removed: " duplicates_removed > "/dev/stderr"
        }
      }'
} > "$OUTPUT"

#==========================================
# SUMMARY
#==========================================

if [ -f "$OUTPUT" ]; then
    TOTAL_LINES=$(wc -l < "$OUTPUT")
    TOTAL_ANNOS=$((TOTAL_LINES - 1))
    FILE_SIZE=$(du -h "$OUTPUT" | cut -f1)
    
    echo "[$(date)] Processing complete!" >&2
    echo "  Total annotations: $TOTAL_ANNOS" >&2
    echo "  Output size: $FILE_SIZE" >&2
    
    echo "  Group breakdown:" >&2
    awk -F'\t' 'NR>1 {print $6}' "$OUTPUT" | sort | uniq -c | sort -rn | while read count group; do
        echo "    $group: $count" >&2
    done
    
    if [ "$LOF_MODE" != "loftee_only" ]; then
        echo "  LoF group detail:" >&2
        awk -F'\t' 'NR>1 && $6=="lof" {print $4}' "$OUTPUT" | sort | uniq -c | sort -rn | while read count lof_val; do
            lof_label=$lof_val
            [ "$lof_val" = "." ] && lof_label="VEP_only"
            [ "$lof_val" = "HC" ] && lof_label="LOFTEE_HC"
            echo "    $lof_label: $count" >&2
        done
    fi
    
    echo "  Output file: $OUTPUT" >&2
else
    echo "[$(date)] Error: Output file was not created!" >&2
    exit 1
fi

exit 0

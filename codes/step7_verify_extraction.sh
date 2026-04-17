#!/bin/bash

# Script: step7_verify_extraction.sh
# Description: Verify PLINK2 extraction quality and completeness
# Usage: step7_verify_extraction.sh <gene_bfiles_dir> <regions_dir> [output_report]

set -euo pipefail

BFILES_DIR="${1:-gene_bfiles}"
REGIONS_DIR="${2:-plink_regions}"
REPORT="${3:-${BFILES_DIR}/verification_report.txt}"

[ ! -d "$BFILES_DIR" ] && { echo "Error: BFILEs directory not found: $BFILES_DIR" >&2; exit 1; }
[ ! -d "$REGIONS_DIR" ] && { echo "Error: Regions directory not found: $REGIONS_DIR" >&2; exit 1; }

echo "[$(date)] Verifying extraction..." >&2

{
  echo "=========================================="
  echo "PLINK2 Extraction Verification Report"
  echo "Generated: $(date)"
  echo "=========================================="
  echo ""
  
  for chr in {1..22} X Y; do
    BFILE="$BFILES_DIR/chr${chr}_genes"
    REGIONS="$REGIONS_DIR/chr${chr}_regions.txt"
    
    if [ ! -f "${BFILE}.bed" ]; then
      continue
    fi
    
    echo "Chromosome $chr"
    echo "----------------------------------------"
    
    # Expected genes
    if [ -f "$REGIONS" ]; then
      EXPECTED_GENES=$(wc -l < "$REGIONS")
      echo "  Expected genes: $EXPECTED_GENES"
    fi
    
    # Actual variants
    VARIANTS=$(wc -l < "${BFILE}.bim")
    echo "  Variants extracted: $VARIANTS"
    
    # Samples
    SAMPLES=$(wc -l < "${BFILE}.fam")
    echo "  Samples: $SAMPLES"
    
    # MAF distribution
    if command -v plink2 &> /dev/null; then
      plink2 --bfile "$BFILE" --freq --out "${BFILE}_freq" --silent 2>/dev/null
      
      if [ -f "${BFILE}_freq.afreq" ]; then
        echo "  MAF distribution:"
        awk 'NR>1 {
          maf = ($5 < 0.5) ? $5 : 1-$5
          if (maf < 0.001) rare++
          else if (maf < 0.01) low++
          else if (maf < 0.05) common++
          else vcommon++
        }
        END {
          print "    Ultra-rare (MAF<0.001): " rare
          print "    Rare (0.001≤MAF<0.01): " low
          print "    Low-freq (0.01≤MAF<0.05): " common
          print "    Common (MAF≥0.05): " vcommon
        }' "${BFILE}_freq.afreq"
        
        rm "${BFILE}_freq.afreq" "${BFILE}_freq.log"
      fi
    fi
    
    # File sizes
    echo "  File sizes:"
    du -h "${BFILE}.bed" | awk '{print "    .bed: " $1}'
    du -h "${BFILE}.bim" | awk '{print "    .bim: " $1}'
    du -h "${BFILE}.fam" | awk '{print "    .fam: " $1}'
    
    echo ""
  done
  
  echo "=========================================="
  echo "Summary Statistics"
  echo "=========================================="
  
  TOTAL_VARS=0
  TOTAL_CHRS=0
  
  for bim in "$BFILES_DIR"/chr*_genes.bim; do
    if [ -f "$bim" ]; then
      VARS=$(wc -l < "$bim")
      TOTAL_VARS=$((TOTAL_VARS + VARS))
      TOTAL_CHRS=$((TOTAL_CHRS + 1))
    fi
  done
  
  echo "  Total chromosomes: $TOTAL_CHRS"
  echo "  Total variants: $TOTAL_VARS"
  
  if [ $TOTAL_CHRS -gt 0 ]; then
    echo "  Average variants/chr: $((TOTAL_VARS / TOTAL_CHRS))"
  fi
  
  echo ""
  echo "=========================================="
  
} > "$REPORT"

cat "$REPORT"

echo "" >&2
echo "Report saved to: $REPORT" >&2

exit 0

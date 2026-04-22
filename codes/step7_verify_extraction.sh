#!/bin/bash

# Script: step7_verify_extraction.sh
# Description: Verify PLINK2a extraction quality and completeness (supports chunked files)
# Usage: step7_verify_extraction.sh <gene_bfiles_dir> <regions_dir> [output_report] [output_format]

set -euo pipefail

BFILES_DIR="${1:-gene_bfiles}"
REGIONS_DIR="${2:-plink_regions}"
REPORT="${3:-${BFILES_DIR}/verification_report.txt}"
OUTPUT_FORMAT="${4:-bfile}"  # Options: bfile, pgen, vcf, bgen

[ ! -d "$BFILES_DIR" ] && { echo "Error: Output directory not found: $BFILES_DIR" >&2; exit 1; }
[ ! -d "$REGIONS_DIR" ] && { echo "Error: Regions directory not found: $REGIONS_DIR" >&2; exit 1; }

# Validate output format
OUTPUT_FORMAT=$(echo "$OUTPUT_FORMAT" | tr '[:upper:]' '[:lower:]')
if [[ ! "$OUTPUT_FORMAT" =~ ^(bfile|pgen|vcf|bgen)$ ]]; then
    echo "ERROR: Invalid output_format. Must be: bfile, pgen, vcf, or bgen" >&2
    exit 1
fi

# Set file extensions based on format
case "$OUTPUT_FORMAT" in
    bfile)
        MAIN_EXT="bed"
        VAR_EXT="bim"
        SAMPLE_EXT="fam"
        ;;
    pgen)
        MAIN_EXT="pgen"
        VAR_EXT="pvar"
        SAMPLE_EXT="psam"
        ;;
    vcf)
        MAIN_EXT="vcf.gz"
        VAR_EXT="vcf.gz"
        SAMPLE_EXT="vcf.gz"
        ;;
    bgen)
        MAIN_EXT="bgen"
        VAR_EXT="bgen"
        SAMPLE_EXT="sample"
        ;;
esac

echo "[$(date)] Verifying extraction..." >&2
echo "  Output format: $OUTPUT_FORMAT" >&2
echo "  Checking for chunked and non-chunked files..." >&2

{
  echo "=========================================="
  echo "PLINK2a Extraction Verification Report"
  echo "Generated: $(date)"
  echo "Output format: $OUTPUT_FORMAT"
  echo "=========================================="
  echo ""
  
  TOTAL_VARS=0
  TOTAL_CHRS=0
  TOTAL_CHUNKS=0
  TOTAL_SAMPLES=0
  
  for chr in {1..22} X Y; do
    
    # Check for non-chunked files
    BFILE="$BFILES_DIR/chr${chr}_genes"
    REGIONS="$REGIONS_DIR/chr${chr}_regions.txt"
    
    # Check for chunked files
    CHUNK_FILES=($BFILES_DIR/chr${chr}_genes_chunk*.${MAIN_EXT} 2>/dev/null)
    HAS_CHUNKS=false
    
    if [ -e "${CHUNK_FILES[0]}" ]; then
      HAS_CHUNKS=true
    fi
    
    # Skip if neither chunked nor non-chunked files exist
    if [ ! -f "${BFILE}.${MAIN_EXT}" ] && [ "$HAS_CHUNKS" = false ]; then
      continue
    fi
    
    echo "Chromosome $chr"
    echo "----------------------------------------"
    
    # Expected genes from regions file
    if [ -f "$REGIONS" ]; then
      EXPECTED_GENES=$(wc -l < "$REGIONS")
      echo "  Expected genes: $EXPECTED_GENES"
    else
      EXPECTED_GENES="N/A"
      echo "  Expected genes: N/A (no region file)"
    fi
    
    TOTAL_CHRS=$((TOTAL_CHRS + 1))
    
    #==========================================
    # Process Non-Chunked Files
    #==========================================
    
    if [ -f "${BFILE}.${MAIN_EXT}" ]; then
      echo ""
      echo "  Non-chunked file: chr${chr}_genes"
      echo "  ........................................"
      
      # Count variants
      case "$OUTPUT_FORMAT" in
        bfile)
          VARIANTS=$(wc -l < "${BFILE}.bim")
          SAMPLES=$(wc -l < "${BFILE}.fam")
          ;;
        pgen)
          VARIANTS=$(grep -v "^#" "${BFILE}.pvar" | wc -l)
          SAMPLES=$(tail -n +2 "${BFILE}.psam" | wc -l)
          ;;
        vcf)
          VARIANTS=$(zgrep -v "^#" "${BFILE}.vcf.gz" | wc -l)
          SAMPLES=$(zgrep "^#CHROM" "${BFILE}.vcf.gz" | awk '{print NF-9}')
          ;;
        bgen)
          # Get from log or count via bgenix
          if [ -f "${BFILE}.log" ]; then
            VARIANTS=$(grep -i "variants written" "${BFILE}.log" | grep -oP '\d+' | head -1 || echo "0")
          else
            VARIANTS="N/A"
          fi
          if [ -f "${BFILE}.sample" ]; then
            SAMPLES=$(tail -n +3 "${BFILE}.sample" | wc -l)
          else
            SAMPLES="N/A"
          fi
          ;;
      esac
      
      echo "    Variants extracted: $VARIANTS"
      echo "    Samples: $SAMPLES"
      TOTAL_VARS=$((TOTAL_VARS + VARIANTS))
      TOTAL_SAMPLES=$SAMPLES
      
      # MAF distribution (only for bfile/pgen)
      if [[ "$OUTPUT_FORMAT" =~ ^(bfile|pgen)$ ]] && command -v plink2a &> /dev/null; then
        if [ "$OUTPUT_FORMAT" = "bfile" ]; then
          plink2a --bfile "$BFILE" --freq --out "${BFILE}_freq" --silent 2>/dev/null
        else
          plink2a --pfile "$BFILE" --freq --out "${BFILE}_freq" --silent 2>/dev/null
        fi
        
        if [ -f "${BFILE}_freq.afreq" ]; then
          echo "    MAF distribution:"
          awk 'NR>1 {
            maf = ($5 < 0.5) ? $5 : 1-$5
            if (maf < 0.001) rare++
            else if (maf < 0.01) low++
            else if (maf < 0.05) common++
            else vcommon++
          }
          END {
            print "      Ultra-rare (MAF<0.001): " rare+0
            print "      Rare (0.001≤MAF<0.01): " low+0
            print "      Low-freq (0.01≤MAF<0.05): " common+0
            print "      Common (MAF≥0.05): " vcommon+0
          }' "${BFILE}_freq.afreq"
          
          rm -f "${BFILE}_freq.afreq" "${BFILE}_freq.log"
        fi
      fi
      
      # File sizes
      echo "    File sizes:"
      case "$OUTPUT_FORMAT" in
        bfile)
          du -h "${BFILE}.bed" 2>/dev/null | awk '{print "      .bed: " $1}'
          du -h "${BFILE}.bim" 2>/dev/null | awk '{print "      .bim: " $1}'
          du -h "${BFILE}.fam" 2>/dev/null | awk '{print "      .fam: " $1}'
          ;;
        pgen)
          du -h "${BFILE}.pgen" 2>/dev/null | awk '{print "      .pgen: " $1}'
          du -h "${BFILE}.pvar" 2>/dev/null | awk '{print "      .pvar: " $1}'
          du -h "${BFILE}.psam" 2>/dev/null | awk '{print "      .psam: " $1}'
          ;;
        vcf)
          du -h "${BFILE}.vcf.gz" 2>/dev/null | awk '{print "      .vcf.gz: " $1}'
          [ -f "${BFILE}.vcf.gz.tbi" ] && du -h "${BFILE}.vcf.gz.tbi" 2>/dev/null | awk '{print "      .tbi: " $1}'
          ;;
        bgen)
          du -h "${BFILE}.bgen" 2>/dev/null | awk '{print "      .bgen: " $1}'
          du -h "${BFILE}.sample" 2>/dev/null | awk '{print "      .sample: " $1}'
          [ -f "${BFILE}.bgen.bgi" ] && du -h "${BFILE}.bgen.bgi" 2>/dev/null | awk '{print "      .bgi: " $1}'
          ;;
      esac
    fi
    
    #==========================================
    # Process Chunked Files
    #==========================================
    
    if [ "$HAS_CHUNKS" = true ]; then
      echo ""
      echo "  Chunked files detected"
      echo "  ........................................"
      
      CHUNK_COUNT=0
      CHR_CHUNK_VARS=0
      
      for chunk_file in "$BFILES_DIR"/chr${chr}_genes_chunk*.${MAIN_EXT}; do
        [ ! -f "$chunk_file" ] && continue
        
        CHUNK_COUNT=$((CHUNK_COUNT + 1))
        CHUNK_BASE="${chunk_file%.${MAIN_EXT}}"
        CHUNK_NAME=$(basename "$CHUNK_BASE")
        
        # Count variants per chunk
        case "$OUTPUT_FORMAT" in
          bfile)
            CHUNK_VARS=$(wc -l < "${CHUNK_BASE}.bim")
            CHUNK_SAMPLES=$(wc -l < "${CHUNK_BASE}.fam")
            ;;
          pgen)
            CHUNK_VARS=$(grep -v "^#" "${CHUNK_BASE}.pvar" | wc -l)
            CHUNK_SAMPLES=$(tail -n +2 "${CHUNK_BASE}.psam" | wc -l)
            ;;
          vcf)
            CHUNK_VARS=$(zgrep -v "^#" "${CHUNK_BASE}.vcf.gz" | wc -l)
            CHUNK_SAMPLES=$(zgrep "^#CHROM" "${CHUNK_BASE}.vcf.gz" | awk '{print NF-9}')
            ;;
          bgen)
            if [ -f "${CHUNK_BASE}.log" ]; then
              CHUNK_VARS=$(grep -i "variants written" "${CHUNK_BASE}.log" | grep -oP '\d+' | head -1 || echo "0")
            else
              CHUNK_VARS="N/A"
            fi
            if [ -f "${CHUNK_BASE}.sample" ]; then
              CHUNK_SAMPLES=$(tail -n +3 "${CHUNK_BASE}.sample" | wc -l)
            else
              CHUNK_SAMPLES="N/A"
            fi
            ;;
        esac
        
        echo "    Chunk $CHUNK_COUNT: $CHUNK_NAME"
        echo "      Variants: $CHUNK_VARS"
        echo "      Samples: $CHUNK_SAMPLES"
        
        if [[ "$CHUNK_VARS" =~ ^[0-9]+$ ]]; then
          CHR_CHUNK_VARS=$((CHR_CHUNK_VARS + CHUNK_VARS))
        fi
        TOTAL_SAMPLES=$CHUNK_SAMPLES
      done
      
      echo "    ........................................"
      echo "    Total chunks: $CHUNK_COUNT"
      echo "    Total variants (all chunks): $CHR_CHUNK_VARS"
      
      TOTAL_VARS=$((TOTAL_VARS + CHR_CHUNK_VARS))
      TOTAL_CHUNKS=$((TOTAL_CHUNKS + CHUNK_COUNT))
      
      # Calculate genes per chunk
      if [ "$EXPECTED_GENES" != "N/A" ] && [ $CHUNK_COUNT -gt 0 ]; then
        AVG_GENES_PER_CHUNK=$((EXPECTED_GENES / CHUNK_COUNT))
        echo "    Average genes per chunk: ~$AVG_GENES_PER_CHUNK"
      fi
    fi
    
    echo ""
  done
  
  echo "=========================================="
  echo "Summary Statistics"
  echo "=========================================="
  echo "  Output format: $OUTPUT_FORMAT"
  echo "  Total chromosomes: $TOTAL_CHRS"
  
  if [ $TOTAL_CHUNKS -gt 0 ]; then
    echo "  Total chunks: $TOTAL_CHUNKS"
    echo "  Average chunks per chromosome: $((TOTAL_CHUNKS / TOTAL_CHRS))"
  else
    echo "  Chunking: No (single file per chromosome)"
  fi
  
  echo "  Total variants: $TOTAL_VARS"
  echo "  Total samples: $TOTAL_SAMPLES"
  
  if [ $TOTAL_CHRS -gt 0 ]; then
    echo "  Average variants/chr: $((TOTAL_VARS / TOTAL_CHRS))"
  fi
  
  echo ""
  echo "=========================================="
  echo "File Inventory"
  echo "=========================================="
  
  if [ $TOTAL_CHUNKS -eq 0 ]; then
    echo "  Non-chunked files:"
    ls -lh "$BFILES_DIR"/chr*_genes.${MAIN_EXT} 2>/dev/null | awk '{print "    " $9 " (" $5 ")"}'
  else
    echo "  Chunked files:"
    ls -lh "$BFILES_DIR"/chr*_genes_chunk*.${MAIN_EXT} 2>/dev/null | awk '{print "    " $9 " (" $5 ")"}'
  fi
  
  echo ""
  echo "=========================================="
  echo "Verification Complete"
  echo "=========================================="
  
} > "$REPORT"

cat "$REPORT"

echo "" >&2
echo "Report saved to: $REPORT" >&2

exit 0

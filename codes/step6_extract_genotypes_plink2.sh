#!/bin/bash

# Script: step6_extract_genotypes_plink2.sh
# Description: Extract genotypes for gene regions using PLINK2 (one BFILE per chromosome)
# Usage: step6_extract_genotypes_plink2.sh <bfile_dir> <regions_dir> <output_dir> [threads]

set -euo pipefail

BFILE_DIR="${1:-}"
REGIONS_DIR="${2:-}"
OUTPUT_DIR="${3:-gene_bfiles}"
THREADS="${4:-4}"

if [ -z "$BFILE_DIR" ] || [ -z "$REGIONS_DIR" ]; then
    cat << 'EOF'
Usage: step6_extract_genotypes_plink2.sh <bfile_dir> <regions_dir> <output_dir> [threads]

Description:
  Extract genotypes for genes using PLINK2.
  Creates ONE BFILE PER CHROMOSOME containing all genes on that chromosome.
  This is more efficient for SAIGE-GENE than per-gene BFILEs.

Arguments:
  bfile_dir      Directory with chromosome BFILE sets (chr1.bed/bim/fam, chr2.bed/bim/fam, ...)
  regions_dir    Directory with chr*_regions.txt files (from match_genes_to_groups.sh)
  output_dir     Output directory for extracted BFILEs
  threads        Number of threads for PLINK2 (default: 4)

Input BFILE naming:
  Expected: chr1.bed/bim/fam, chr2.bed/bim/fam, ..., chr22.bed/bim/fam
  Or: chr1, chr2, ..., chr22 (without extension)

Output:
  - chr1_genes.bed/bim/fam   Genotypes for all genes on chr1
  - chr2_genes.bed/bim/fam   Genotypes for all genes on chr2
  - ...
  - extraction_summary.txt   Summary of extraction

Examples:
  # Standard usage
  ./step6_extract_genotypes_plink2.sh /data/bfiles plink_regions gene_bfiles

  # With more threads
  ./step6_extract_genotypes_plink2.sh /data/bfiles plink_regions gene_bfiles 16

  # Use in SAIGE
  Rscript step2_SPAtests.R \
    --bgenFile=gene_bfiles/chr1_genes.bgen \
    --groupFile=chr1_groups.txt \
    ...
EOF
    exit 1
fi

[ ! -d "$BFILE_DIR" ] && { echo "Error: BFILE directory not found: $BFILE_DIR" >&2; exit 1; }
[ ! -d "$REGIONS_DIR" ] && { echo "Error: Regions directory not found: $REGIONS_DIR" >&2; exit 1; }

# Check for PLINK2
if ! command -v plink2 &> /dev/null; then
    echo "Error: plink2 not found. Please install PLINK2." >&2
    echo "Download: https://www.cog-genomics.org/plink/2.0/" >&2
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

echo "[$(date)] Extracting genotypes with PLINK2..." >&2
echo "  BFILE directory: $BFILE_DIR" >&2
echo "  Regions directory: $REGIONS_DIR" >&2
echo "  Output directory: $OUTPUT_DIR" >&2
echo "  Threads: $THREADS" >&2
echo "" >&2

# Summary file
SUMMARY="$OUTPUT_DIR/extraction_summary.txt"
echo -e "Chr\tGenes\tVariants_Input\tVariants_Output\tSamples\tStatus" > "$SUMMARY"

TOTAL_EXTRACTED=0
TOTAL_FAILED=0

#==========================================
# Process each chromosome
#==========================================

for chr in {1..22} X Y; do
  REGION_FILE="$REGIONS_DIR/chr${chr}_regions.txt"
  
  # Skip if no region file
  if [ ! -f "$REGION_FILE" ]; then
    echo "Skipping chr${chr}: no region file" >&2
    continue
  fi
  
  # Check for BFILE (try multiple naming conventions)
  BFILE=""
  if [ -f "$BFILE_DIR/chr${chr}.bed" ]; then
    BFILE="$BFILE_DIR/chr${chr}"
  elif [ -f "$BFILE_DIR/chr${chr}.bgen" ]; then
    BFILE="$BFILE_DIR/chr${chr}"
  elif [ -f "$BFILE_DIR/chr${chr}_merged.bed" ]; then
    BFILE="$BFILE_DIR/chr${chr}_merged"
  else
    echo "WARNING: BFILE not found for chr${chr}" >&2
    echo -e "${chr}\tNA\tNA\tNA\tNA\tNOT_FOUND" >> "$SUMMARY"
    continue
  fi
  
  GENE_COUNT=$(wc -l < "$REGION_FILE")
  
  echo "[$(date)] Processing chr${chr}..." >&2
  echo "  BFILE: $BFILE" >&2
  echo "  Genes: $GENE_COUNT" >&2
  
  # Extract using PLINK2
  OUTPUT_PREFIX="$OUTPUT_DIR/chr${chr}_genes"
  
  # Method 1: Use --extract bed1 (more efficient)
  plink2 --bfile "$BFILE" \
         --extract bed1 "$REGION_FILE" \
         --make-bed \
         --threads "$THREADS" \
         --out "$OUTPUT_PREFIX" \
         --silent 2>&1 | grep -v "^PLINK" || true
  
  # Check if extraction succeeded
  if [ -f "${OUTPUT_PREFIX}.bed" ]; then
    # Count variants and samples
    VARS_OUTPUT=$(wc -l < "${OUTPUT_PREFIX}.bim")
    SAMPLES=$(wc -l < "${OUTPUT_PREFIX}.fam")
    
    # Get input variant count
    VARS_INPUT=$(wc -l < "${BFILE}.bim" 2>/dev/null || echo "NA")
    
    echo "  ✓ Extracted: $VARS_OUTPUT variants, $SAMPLES samples" >&2
    echo -e "${chr}\t${GENE_COUNT}\t${VARS_INPUT}\t${VARS_OUTPUT}\t${SAMPLES}\tSUCCESS" >> "$SUMMARY"
    
    TOTAL_EXTRACTED=$((TOTAL_EXTRACTED + 1))
  else
    echo "  ✗ Extraction failed" >&2
    echo -e "${chr}\t${GENE_COUNT}\tNA\tNA\tNA\tFAILED" >> "$SUMMARY"
    TOTAL_FAILED=$((TOTAL_FAILED + 1))
  fi
  
  echo "" >&2
done

#==========================================
# Summary
#==========================================

echo "==========================================" >&2
echo "Genotype Extraction Complete" >&2
echo "==========================================" >&2
echo "Chromosomes processed: $((TOTAL_EXTRACTED + TOTAL_FAILED))" >&2
echo "Successfully extracted: $TOTAL_EXTRACTED" >&2
echo "Failed: $TOTAL_FAILED" >&2
echo "" >&2

echo "Output BFILEs:" >&2
ls -lh "$OUTPUT_DIR"/chr*_genes.{bed,bim,fam} 2>/dev/null | \
  awk '{print "  " $9 " (" $5 ")"}' >&2

echo "" >&2
echo "Summary table:" >&2
column -t -s $'\t' "$SUMMARY" | sed 's/^/  /' >&2

echo "" >&2
echo "==========================================" >&2
echo "Next steps:" >&2
echo "" >&2
echo "1. Verify extraction:" >&2
echo "   for chr in {1..22}; do" >&2
echo "     plink2 --bfile $OUTPUT_DIR/chr\${chr}_genes --freq --out $OUTPUT_DIR/chr\${chr}_freq" >&2
echo "   done" >&2
echo "" >&2
echo "2. Use in SAIGE-GENE:" >&2
echo "   Rscript step2_SPAtests.R \\" >&2
echo "     --bgenFile=$OUTPUT_DIR/chr1_genes.bgen \\" >&2
echo "     --groupFile=chr1_groups.txt \\" >&2
echo "     --sampleFile=sample.txt \\" >&2
echo "     ..." >&2
echo "==========================================" >&2

exit 0

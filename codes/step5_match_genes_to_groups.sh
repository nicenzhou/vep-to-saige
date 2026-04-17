#!/bin/bash

# Script: step5_match_genes_to_groups.sh
# Description: Match gene coordinates with SAIGE group files and create PLINK2 region files
# Usage: step5_match_genes_to_groups.sh <gene_coords_dir> <group_file> <output_dir> [buffer_kb]

set -euo pipefail

GENE_COORDS_DIR="${1:-}"
GROUP_FILE="${2:-}"
OUTPUT_DIR="${3:-plink_regions}"
BUFFER_KB="${4:-10}"

if [ -z "$GENE_COORDS_DIR" ] || [ -z "$GROUP_FILE" ]; then
    cat << 'EOF'
Usage: step5_match_genes_to_groups.sh <gene_coords_dir> <group_file> <output_dir> [buffer_kb]

Description:
  Match genes in SAIGE group file with Ensembl coordinates.
  Create PLINK2-compatible region files for genotype extraction.

Arguments:
  gene_coords_dir    Directory with chr*_genes.txt files
  group_file         SAIGE group file (all chromosomes or single chromosome)
  output_dir         Output directory for PLINK2 region files
  buffer_kb          Buffer around genes in kb (default: 10)

Input Files (gene_coords_dir):
  - chr1_genes.txt, chr2_genes.txt, ..., chr22_genes.txt
  - Format: Gene_Symbol Gene_ID Chr Start End Strand Biotype

Input Files (group_file):
  - SAIGE group file format (space-delimited)
  - GENE var variant1 variant2 ...
  - GENE anno lof missense ...

Output Files:
  - chr*_regions.txt        PLINK2 region files (chr:start-end format)
  - chr*_gene_list.txt      Gene symbols for each chromosome
  - matched_genes.txt       All matched genes with coordinates
  - missing_genes.txt       Genes in groups but not in coordinates

Examples:
  # Process all chromosomes
  ./step5_match_genes_to_groups.sh gene_coords all_genes_groups.txt plink_regions

  # Single chromosome
  ./step5_match_genes_to_groups.sh gene_coords chr1_groups.txt plink_regions

  # Custom buffer
  ./step5_match_genes_to_groups.sh gene_coords all_genes_groups.txt plink_regions 50
EOF
    exit 1
fi

[ ! -d "$GENE_COORDS_DIR" ] && { echo "Error: Gene coords directory not found: $GENE_COORDS_DIR" >&2; exit 1; }
[ ! -f "$GROUP_FILE" ] && { echo "Error: Group file not found: $GROUP_FILE" >&2; exit 1; }

mkdir -p "$OUTPUT_DIR"

echo "[$(date)] Matching genes to coordinates..." >&2
echo "  Gene coordinates: $GENE_COORDS_DIR" >&2
echo "  Group file: $GROUP_FILE" >&2
echo "  Output: $OUTPUT_DIR" >&2
echo "  Buffer: ${BUFFER_KB}kb" >&2
echo "" >&2

#==========================================
# Extract gene list from group file
#==========================================

echo "Extracting genes from group file..." >&2

awk '{if(NR%2==1) print $1}' "$GROUP_FILE" | \
  sort -u > "$OUTPUT_DIR/genes_in_groups.txt"

TOTAL_IN_GROUPS=$(wc -l < "$OUTPUT_DIR/genes_in_groups.txt")
echo "  Genes in group file: $TOTAL_IN_GROUPS" >&2

#==========================================
# Match with coordinates (per chromosome)
#==========================================

echo "" >&2
echo "Matching genes with coordinates by chromosome..." >&2

TOTAL_MATCHED=0
TOTAL_MISSING=0
> "$OUTPUT_DIR/matched_genes.txt"
> "$OUTPUT_DIR/missing_genes.txt"

# Process each chromosome coordinate file
for COORD_FILE in "$GENE_COORDS_DIR"/chr*_genes.txt; do
  [ ! -f "$COORD_FILE" ] && continue
  
  CHR=$(basename "$COORD_FILE" | sed 's/chr//; s/_genes.txt//')
  
  echo "  Processing chr${CHR}..." >&2
  
  # Match genes
  awk -F'\t' -v genes="$OUTPUT_DIR/genes_in_groups.txt" -v chr="$CHR" '
  BEGIN {
    # Load genes to match
    while ((getline < genes) > 0) {
      wanted[$1] = 1
    }
    close(genes)
  }
  
  NR==1 {next}  # Skip header
  
  {
    gene = $1
    if (gene in wanted) {
      print $0
      matched[gene] = 1
    }
  }
  
  END {
    # Track missing genes
    for (g in wanted) {
      if (!(g in matched)) {
        missing[g] = 1
      }
    }
    
    # Print missing to stderr for this chr
    n_missing = 0
    for (g in missing) {
      n_missing++
    }
    if (n_missing > 0) {
      print "    Missing: " n_missing " genes" > "/dev/stderr"
    }
  }
  ' "$COORD_FILE" > "$OUTPUT_DIR/chr${CHR}_matched.txt"
  
  MATCHED=$(wc -l < "$OUTPUT_DIR/chr${CHR}_matched.txt")
  
  if [ "$MATCHED" -gt 0 ]; then
    TOTAL_MATCHED=$((TOTAL_MATCHED + MATCHED))
    
    # Append to all matched genes
    cat "$OUTPUT_DIR/chr${CHR}_matched.txt" >> "$OUTPUT_DIR/matched_genes.txt"
    
    # Create PLINK2 region file
    awk -F'\t' -v buffer="$BUFFER_KB" -v chr="$CHR" '
    BEGIN {
      buffer_bp = buffer * 1000
    }
    {
      gene = $1
      start = $4 - buffer_bp
      end = $5 + buffer_bp
      
      if (start < 1) start = 1
      
      # PLINK2 format: chr:start-end
      print chr ":" start "-" end "\t" gene
    }
    ' "$OUTPUT_DIR/chr${CHR}_matched.txt" > "$OUTPUT_DIR/chr${CHR}_regions.txt"
    
    # Create gene list
    cut -f1 "$OUTPUT_DIR/chr${CHR}_matched.txt" > "$OUTPUT_DIR/chr${CHR}_gene_list.txt"
    
    echo "    Matched: $MATCHED genes → chr${CHR}_regions.txt" >&2
  else
    echo "    Matched: 0 genes (skipped)" >&2
  fi
  
  # Cleanup temp file
  rm -f "$OUTPUT_DIR/chr${CHR}_matched.txt"
done

#==========================================
# Find missing genes across all chromosomes
#==========================================

echo "" >&2
echo "Identifying missing genes..." >&2

# All genes that were matched
cat "$OUTPUT_DIR/matched_genes.txt" | cut -f1 | sort -u > "$OUTPUT_DIR/all_matched.txt"

# Find genes in groups but not matched
comm -23 \
  <(sort "$OUTPUT_DIR/genes_in_groups.txt") \
  <(sort "$OUTPUT_DIR/all_matched.txt") \
  > "$OUTPUT_DIR/missing_genes.txt"

TOTAL_MISSING=$(wc -l < "$OUTPUT_DIR/missing_genes.txt")

#==========================================
# Summary
#==========================================

echo "" >&2
echo "==========================================" >&2
echo "Gene Matching Complete" >&2
echo "==========================================" >&2
echo "Genes in group file: $TOTAL_IN_GROUPS" >&2
echo "Genes matched: $TOTAL_MATCHED" >&2
echo "Genes missing: $TOTAL_MISSING" >&2
echo "" >&2
echo "Output files:" >&2
ls -1 "$OUTPUT_DIR"/chr*_regions.txt 2>/dev/null | while read f; do
  count=$(wc -l < "$f")
  echo "  $(basename $f): $count genes" >&2
done
echo "" >&2

if [ "$TOTAL_MISSING" -gt 0 ]; then
  echo "WARNING: Missing genes saved to:" >&2
  echo "  $OUTPUT_DIR/missing_genes.txt" >&2
  echo "" >&2
  echo "First 10 missing genes:" >&2
  head -10 "$OUTPUT_DIR/missing_genes.txt" | sed 's/^/  /' >&2
  echo "" >&2
fi

echo "Sample region file (chr1, first 5):" >&2
if [ -f "$OUTPUT_DIR/chr1_regions.txt" ]; then
  head -5 "$OUTPUT_DIR/chr1_regions.txt" | sed 's/^/  /' >&2
fi

echo "==========================================" >&2
echo "" >&2
echo "Next step: Extract genotypes with PLINK2" >&2
echo "  ./extract_genotypes_plink2.sh <bfile_prefix> $OUTPUT_DIR <output_dir>" >&2

# Cleanup
rm -f "$OUTPUT_DIR/genes_in_groups.txt"
rm -f "$OUTPUT_DIR/all_matched.txt"

exit 0

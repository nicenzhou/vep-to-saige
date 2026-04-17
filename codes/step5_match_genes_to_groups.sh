#!/bin/bash

# Script: step5_match_genes_to_groups.sh
# Description: Match gene coordinates with SAIGE group files and create PLINK2 region files
# Usage: step5_match_genes_to_groups.sh <gene_coords_dir> <group_file> <output_dir> [buffer_kb] [force_regen]

set -euo pipefail

#==========================================
# Parse Arguments
#==========================================

GENE_COORDS_DIR="${1:-}"
GROUP_FILE="${2:-}"
OUTPUT_DIR="${3:-plink_regions}"
BUFFER_KB="${4:-10}"
FORCE_REGEN="${5:-no}"

#==========================================
# Help Message
#==========================================

if [ -z "$$GENE_COORDS_DIR" ] || [ -z "$$GROUP_FILE" ]; then
    cat << 'EOF'
Usage: step5_match_genes_to_groups.sh <gene_coords_dir> <group_file> <output_dir> [buffer_kb] [force_regen]

Description:
  Match genes in SAIGE group file with Ensembl coordinates.
  Create PLINK2-compatible region files for genotype extraction.
  Optionally regenerate missing genes from variant positions.

Arguments:
  gene_coords_dir    Directory with chr*_genes.txt files
  group_file         SAIGE group file (from Step 3)
  output_dir         Output directory for region files (default: plink_regions)
  buffer_kb          Buffer around genes in kb (default: 10)
  force_regen        Force regenerate missing genes: yes/no (default: no)

Input Files (gene_coords_dir):
  - chr1_genes.txt, chr2_genes.txt, ..., chr22_genes.txt
  - Format: Gene_Symbol Gene_ID Chr Start End Strand Biotype

Input Files (group_file):
  - SAIGE group file format (space-delimited)
  - Line 1: GENE var variant1 variant2 ...
  - Line 2: GENE anno lof missense ...
  - Variant format: chr:pos:ref:alt

Output Files:
  - chr*_regions.txt        PLINK2 region files (chr:start-end format)
  - chr*_gene_list.txt      Gene symbols for each chromosome
  - matched_genes.txt       All matched genes with coordinates
  - missing_genes.txt       Genes in groups but not in coordinates
  - regenerated_genes.txt   Regenerated coordinates for missing genes (if force_regen=yes)
  - matching_summary.txt    Detailed summary statistics

Force Regenerate Mode:
  When force_regen=yes, missing genes are regenerated using variant positions:
  - Start position: First variant position - buffer_kb - 10kb (soft buffer)
  - End position: Last variant position + buffer_kb + 10kb (soft buffer)
  - Total buffer: User-defined buffer + 10kb safety margin on each side
  - Missing genes are added to regenerated_genes.txt

Examples:
  # Basic usage with default settings
  ./step5_match_genes_to_groups.sh gene_coords all_genes_groups.txt plink_regions

  # Custom buffer of 50kb
  ./step5_match_genes_to_groups.sh gene_coords all_genes_groups.txt plink_regions 50

  # Force regenerate missing genes with default 10kb buffer
  ./step5_match_genes_to_groups.sh gene_coords all_genes_groups.txt plink_regions 10 yes

  # Force regenerate with 50kb buffer
  ./step5_match_genes_to_groups.sh gene_coords all_genes_groups.txt plink_regions 50 yes

Next Steps:
  After running this script, extract genotypes using:
  ./extract_genotypes_plink2.sh <bfile_prefix> plink_regions <output_dir>

EOF
    exit 1
fi

#==========================================
# Validate Inputs
#==========================================

if [ ! -d "$GENE_COORDS_DIR" ]; then
    echo "ERROR: Gene coordinates directory not found: $GENE_COORDS_DIR" >&2
    exit 1
fi

if [ ! -f "$GROUP_FILE" ]; then
    echo "ERROR: Group file not found: $GROUP_FILE" >&2
    exit 1
fi

# Validate buffer is numeric
if ! [[ "$$BUFFER_KB" =~ ^[0-9]+$$ ]]; then
    echo "ERROR: Buffer must be a positive integer (got: $BUFFER_KB)" >&2
    exit 1
fi

# Validate force_regen option
FORCE_REGEN=$$(echo "$$FORCE_REGEN" | tr '[:upper:]' '[:lower:]')
if [[ ! "$$FORCE_REGEN" =~ ^(yes|no)$$ ]]; then
    echo "ERROR: force_regen must be 'yes' or 'no' (got: $FORCE_REGEN)" >&2
    exit 1
fi

# Check for coordinate files
COORD_COUNT=$$(find "$$GENE_COORDS_DIR" -name "chr*_genes.txt" -type f | wc -l)
if [ "$COORD_COUNT" -eq 0 ]; then
    echo "ERROR: No chr*_genes.txt files found in $GENE_COORDS_DIR" >&2
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

#==========================================
# Log Header
#==========================================

LOG_FILE="$OUTPUT_DIR/matching.log"
SUMMARY_FILE="$OUTPUT_DIR/matching_summary.txt"

exec > >(tee -a "$LOG_FILE") 2>&1

echo "=========================================="
echo "Gene-to-Group Matching Pipeline"
echo "=========================================="
echo "Started: $(date)"
echo ""
echo "Configuration:"
echo "  Gene coordinates: $GENE_COORDS_DIR"
echo "  Group file: $GROUP_FILE"
echo "  Output directory: $OUTPUT_DIR"
echo "  Buffer size: $${BUFFER_KB}kb (±$$(echo "$BUFFER_KB * 1000" | bc) bp)"
echo "  Force regenerate: $FORCE_REGEN"
if [ "$FORCE_REGEN" = "yes" ]; then
    echo "  Regenerate buffer: $${BUFFER_KB}kb + 10kb safety = $$((BUFFER_KB + 10))kb total"
fi
echo "  Coordinate files found: $COORD_COUNT"
echo ""

#==========================================
# Extract Gene List from Group File
#==========================================

echo "Step 1: Extracting genes from group file..."

# Extract gene names (every odd line in SAIGE format)
awk '{if(NR%2==1) print \$1}' "$GROUP_FILE" | \
  grep -v '^$' | \
  sort -u > "$OUTPUT_DIR/genes_in_groups.txt"

TOTAL_IN_GROUPS=$$(wc -l < "$$OUTPUT_DIR/genes_in_groups.txt")
echo "  Total unique genes in group file: $TOTAL_IN_GROUPS"
echo ""

if [ "$TOTAL_IN_GROUPS" -eq 0 ]; then
    echo "ERROR: No genes extracted from group file" >&2
    exit 1
fi

#==========================================
# Match Genes with Coordinates by Chromosome
#==========================================

echo "Step 2: Matching genes with coordinates by chromosome..."
echo ""

TOTAL_MATCHED=0
TOTAL_MISSING=0

# Initialize output files
> "$OUTPUT_DIR/matched_genes.txt"
> "$OUTPUT_DIR/missing_genes.txt"

# Create summary header
{
    echo "Chromosome Matching Summary"
    echo "Generated: $(date)"
    echo ""
    printf "%-5s %-12s %-12s %-12s\n" "Chr" "Genes_Found" "Regions" "Status"
    echo "--------------------------------------------------------"
} > "$SUMMARY_FILE"

# Process each chromosome coordinate file
for COORD_FILE in "$GENE_COORDS_DIR"/chr*_genes.txt; do
    [ ! -f "$COORD_FILE" ] && continue
    
    # Extract chromosome number/name
    CHR=$$(basename "$$COORD_FILE" | sed 's/chr//; s/_genes.txt//')
    
    echo "  Processing chr${CHR}..."
    
    # Match genes for this chromosome
    awk -F'\t' -v genes_file="$OUTPUT_DIR/genes_in_groups.txt" '
    BEGIN {
        # Load genes to match
        while ((getline < genes_file) > 0) {
            wanted[\$1] = 1
        }
        close(genes_file)
    }
    
    NR==1 {next}  # Skip header line
    
    {
        gene = \$1
        if (gene in wanted) {
            print \$0
            matched[gene] = 1
        }
    }
    ' "$$COORD_FILE" > "$$OUTPUT_DIR/chr${CHR}_matched_temp.txt"
    
    MATCHED=$$(wc -l < "$$OUTPUT_DIR/chr${CHR}_matched_temp.txt")
    
    if [ "$MATCHED" -gt 0 ]; then
        TOTAL_MATCHED=$((TOTAL_MATCHED + MATCHED))
        
        # Append to all matched genes
        cat "$$OUTPUT_DIR/chr$${CHR}_matched_temp.txt" >> "$OUTPUT_DIR/matched_genes.txt"
        
        # Create PLINK2 region file with buffer
        awk -F'\t' -v buffer="$$BUFFER_KB" -v chr="$$CHR" '
        BEGIN {
            buffer_bp = buffer * 1000
        }
        {
            gene = \$1
            gene_id = \$2
            start = \$4 - buffer_bp
            end = \$5 + buffer_bp
            
            # Ensure start is not negative
            if (start < 1) start = 1
            
            # PLINK2 format: chr:start-end
            region = chr ":" start "-" end
            
            # Output: region gene_symbol gene_id
            print region "\t" gene "\t" gene_id
        }
        ' "$$OUTPUT_DIR/chr$${CHR}_matched_temp.txt" > "$$OUTPUT_DIR/chr$${CHR}_regions.txt"
        
        # Create simple gene list
        cut -f1 "$$OUTPUT_DIR/chr$${CHR}_matched_temp.txt" > "$$OUTPUT_DIR/chr$${CHR}_gene_list.txt"
        
        REGION_COUNT=$$(wc -l < "$$OUTPUT_DIR/chr${CHR}_regions.txt")
        
        echo "    ✓ Matched: $$MATCHED genes → $$REGION_COUNT regions"
        printf "%-5s %-12s %-12s %-12s\n" "$$CHR" "$$MATCHED" "$$REGION_COUNT" "Success" >> "$$SUMMARY_FILE"
    else
        echo "    ✗ Matched: 0 genes (skipped)"
        printf "%-5s %-12s %-12s %-12s\n" "$$CHR" "0" "0" "No matches" >> "$$SUMMARY_FILE"
    fi
    
    # Cleanup temp file
    rm -f "$$OUTPUT_DIR/chr$${CHR}_matched_temp.txt"
done

echo ""

#==========================================
# Identify Missing Genes
#==========================================

echo "Step 3: Identifying missing genes..."

# Extract all matched gene symbols
if [ -s "$OUTPUT_DIR/matched_genes.txt" ]; then
    cut -f1 "$$OUTPUT_DIR/matched_genes.txt" | sort -u > "$$OUTPUT_DIR/all_matched_genes.txt"
    
    # Find genes in group file but not matched
    comm -23 \
        <(sort "$OUTPUT_DIR/genes_in_groups.txt") \
        <(sort "$OUTPUT_DIR/all_matched_genes.txt") \
        > "$OUTPUT_DIR/missing_genes.txt"
else
    # No matches at all
    cp "$$OUTPUT_DIR/genes_in_groups.txt" "$$OUTPUT_DIR/missing_genes.txt"
fi

TOTAL_MISSING=$$(wc -l < "$$OUTPUT_DIR/missing_genes.txt")

echo "  Missing genes: $TOTAL_MISSING"
echo ""

#==========================================
# Force Regenerate Missing Genes
#==========================================

TOTAL_REGENERATED=0

if [ "$$FORCE_REGEN" = "yes" ] && [ "$$TOTAL_MISSING" -gt 0 ]; then
    echo "Step 4: Force regenerating missing genes from variant positions..."
    echo ""
    
    # Calculate total buffer (user buffer + 10kb safety)
    SOFT_BUFFER_KB=$((BUFFER_KB + 10))
    SOFT_BUFFER_BP=$((SOFT_BUFFER_KB * 1000))
    
    echo "  Total buffer for regeneration: $${SOFT_BUFFER_KB}kb ($${BUFFER_KB}kb user + 10kb safety)"
    echo ""
    
    > "$OUTPUT_DIR/regenerated_genes.txt"
    > "$OUTPUT_DIR/regenerated_regions_temp.txt"
    
    # Process group file to extract variant positions for missing genes
    awk -v missing_file="$OUTPUT_DIR/missing_genes.txt" \
        -v buffer="$SOFT_BUFFER_BP" \
        -v output="$OUTPUT_DIR/regenerated_genes.txt" \
        -v regions="$OUTPUT_DIR/regenerated_regions_temp.txt" '
    BEGIN {
        # Load missing genes
        while ((getline < missing_file) > 0) {
            missing[\$1] = 1
        }
        close(missing_file)
    }
    
    {
        # Process variant lines (odd line numbers)
        if (NR % 2 == 1) {
            gene = \$1
            
            # Only process missing genes
            if (gene in missing) {
                delete positions
                n_vars = 0
                
                # Extract variant positions (skip first two fields: gene name and "var")
                for (i = 3; i <= NF; i++) {
                    variant = $i
                    
                    # Parse chr:pos:ref:alt
                    split(variant, parts, ":")
                    if (length(parts) >= 2) {
                        chr = parts[1]
                        pos = parts[2]
                        
                        # Store position
                        positions[n_vars++] = pos
                    }
                }
                
                # If we have variants, calculate region
                if (n_vars > 0) {
                    # Find min and max positions
                    min_pos = positions[0]
                    max_pos = positions[0]
                    
                    for (j = 1; j < n_vars; j++) {
                        if (positions[j] < min_pos) min_pos = positions[j]
                        if (positions[j] > max_pos) max_pos = positions[j]
                    }
                    
                    # Apply buffer
                    start = min_pos - buffer
                    end = max_pos + buffer
                    
                    if (start < 1) start = 1
                    
                    # Output regenerated gene info
                    print gene "\tREGENERATED\t" chr "\t" start "\t" end "\tNA\tNA\t" n_vars >> output
                    
                    # Output region for PLINK2
                    print chr ":" start "-" end "\t" gene "\tREGENERATED" >> regions
                    
                    regen_count++
                }
            }
        }
    }
    
    END {
        print "  Regenerated " regen_count " genes from variant positions" > "/dev/stderr"
    }
    ' "$GROUP_FILE"
    
    TOTAL_REGENERATED=$$(wc -l < "$$OUTPUT_DIR/regenerated_genes.txt")
    
    echo ""
    
    if [ "$TOTAL_REGENERATED" -gt 0 ]; then
        # Merge regenerated regions into chromosome-specific files
        while IFS=$'\t' read -r region gene source; do
            # Extract chromosome from region (format: chr:start-end)
            CHR=$$(echo "$$region" | cut -d':' -f1)
            
            # Append to chromosome region file
            echo -e "$${region}\t$${gene}\t${source

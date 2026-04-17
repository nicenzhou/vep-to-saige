#!/bin/bash

# Script: step5_match_genes_to_groups.sh
# Description: Match gene coordinates with SAIGE group files and create PLINK2 region files
# Usage: step5_match_genes_to_groups.sh <gene_coords_dir> <group_file> <output_dir> [buffer_kb] [force_regen] [merge_regen]

set -euo pipefail

#==========================================
# Parse Arguments
#==========================================

GENE_COORDS_DIR="${1:-}"
GROUP_FILE="${2:-}"
OUTPUT_DIR="${3:-plink_regions}"
BUFFER_KB="${4:-10}"
FORCE_REGEN="${5:-no}"
MERGE_REGEN="${6:-yes}"  # New parameter: merge regenerated genes or keep separate

#==========================================
# Help Message
#==========================================

if [ -z "$GENE_COORDS_DIR" ] || [ -z "$GROUP_FILE" ]; then
    cat << 'EOF'
Usage: step5_match_genes_to_groups.sh <gene_coords_dir> <group_file> <output_dir> [buffer_kb] [force_regen] [merge_regen]

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
  merge_regen        Merge regenerated genes with matched: yes/no (default: yes)

Input Files (gene_coords_dir):
  - chr1_genes.txt, chr2_genes.txt, ..., chr22_genes.txt
  - Format: Gene_Symbol Gene_ID Chr Start End Strand Biotype

Input Files (group_file):
  - SAIGE group file format (space-delimited)
  - Line 1: GENE var variant1 variant2 ...
  - Line 2: GENE anno lof missense ...
  - Variant format: chr:pos:ref:alt

Output Files (merge_regen=yes):
  - chr*_regions.txt        PLINK2 region files (chr:start-end format)
  - chr*_gene_list.txt      Gene symbols for each chromosome
  - matched_genes.txt       All matched genes with coordinates
  - missing_genes.txt       Genes in groups but not in coordinates
  - regenerated_genes.txt   Regenerated coordinates for missing genes (if force_regen=yes)
  - matching_summary.txt    Detailed summary statistics

Output Files (merge_regen=no):
  - chr*_regions.txt              PLINK2 region files (matched genes only)
  - chr*_regions_recovered.txt    PLINK2 region files (regenerated genes only)
  - chr*_gene_list.txt            Gene symbols (matched genes only)
  - chr*_gene_list_recovered.txt  Gene symbols (regenerated genes only)
  - matched_genes.txt             Matched genes with coordinates
  - regenerated_genes.txt         Regenerated genes with coordinates
  - missing_genes.txt             Genes still missing
  - matching_summary.txt          Detailed summary statistics

Force Regenerate Mode:
  When force_regen=yes, missing genes are regenerated using variant positions:
  - Start position: First variant position - buffer_kb - 10kb (soft buffer)
  - End position: Last variant position + buffer_kb + 10kb (soft buffer)
  - Total buffer: User-defined buffer + 10kb safety margin on each side
  
  When merge_regen=yes: Regenerated genes are merged into main output files
  When merge_regen=no:  Regenerated genes are kept in separate *_recovered.txt files

Examples:
  # Basic usage with default settings (merge regenerated genes)
  ./step5_match_genes_to_groups.sh gene_coords all_genes_groups.txt plink_regions

  # Force regenerate and keep separate
  ./step5_match_genes_to_groups.sh gene_coords all_genes_groups.txt plink_regions 10 yes no

  # Force regenerate and merge
  ./step5_match_genes_to_groups.sh gene_coords all_genes_groups.txt plink_regions 10 yes yes

  # Custom buffer, regenerate, keep separate
  ./step5_match_genes_to_groups.sh gene_coords all_genes_groups.txt plink_regions 50 yes no

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
if ! [[ "$BUFFER_KB" =~ ^[0-9]+$ ]]; then
    echo "ERROR: Buffer must be a positive integer (got: $BUFFER_KB)" >&2
    exit 1
fi

# Validate force_regen option
FORCE_REGEN=$(echo "$FORCE_REGEN" | tr '[:upper:]' '[:lower:]')
if [[ ! "$FORCE_REGEN" =~ ^(yes|no)$ ]]; then
    echo "ERROR: force_regen must be 'yes' or 'no' (got: $FORCE_REGEN)" >&2
    exit 1
fi

# Validate merge_regen option
MERGE_REGEN=$(echo "$MERGE_REGEN" | tr '[:upper:]' '[:lower:]')
if [[ ! "$MERGE_REGEN" =~ ^(yes|no)$ ]]; then
    echo "ERROR: merge_regen must be 'yes' or 'no' (got: $MERGE_REGEN)" >&2
    exit 1
fi

# Check for coordinate files
COORD_COUNT=$(find "$GENE_COORDS_DIR" -name "chr*_genes.txt" -type f | wc -l)
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
echo "  Buffer size: ${BUFFER_KB}kb (±$((BUFFER_KB * 1000)) bp)"
echo "  Force regenerate: $FORCE_REGEN"
echo "  Merge regenerated: $MERGE_REGEN"
if [ "$FORCE_REGEN" = "yes" ]; then
    echo "  Regenerate buffer: ${BUFFER_KB}kb + 10kb safety = $((BUFFER_KB + 10))kb total"
    if [ "$MERGE_REGEN" = "yes" ]; then
        echo "  Regenerated genes will be merged into main output files"
    else
        echo "  Regenerated genes will be saved to separate *_recovered.txt files"
    fi
fi
echo "  Coordinate files found: $COORD_COUNT"
echo ""

#==========================================
# Extract Gene List from Group File
#==========================================

echo "Step 1: Extracting genes from group file..."

# Extract gene names (every odd line in SAIGE format)
awk '{if(NR%2==1) print $1}' "$GROUP_FILE" | \
  grep -v '^$' | \
  sort -u > "$OUTPUT_DIR/genes_in_groups.txt"

TOTAL_IN_GROUPS=$(wc -l < "$OUTPUT_DIR/genes_in_groups.txt")
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
    CHR=$(basename "$COORD_FILE" | sed 's/chr//; s/_genes.txt//')
    
    echo "  Processing chr${CHR}..."
    
    # Match genes for this chromosome
    awk -F'\t' -v genes_file="$OUTPUT_DIR/genes_in_groups.txt" '
    BEGIN {
        # Load genes to match
        while ((getline < genes_file) > 0) {
            wanted[$1] = 1
        }
        close(genes_file)
    }
    
    NR==1 {next}  # Skip header line
    
    {
        gene = $1
        if (gene in wanted) {
            print $0
            matched[gene] = 1
        }
    }
    ' "$COORD_FILE" > "$OUTPUT_DIR/chr${CHR}_matched_temp.txt"
    
    MATCHED=$(wc -l < "$OUTPUT_DIR/chr${CHR}_matched_temp.txt")
    
    if [ "$MATCHED" -gt 0 ]; then
        TOTAL_MATCHED=$((TOTAL_MATCHED + MATCHED))
        
        # Append to all matched genes
        cat "$OUTPUT_DIR/chr${CHR}_matched_temp.txt" >> "$OUTPUT_DIR/matched_genes.txt"
        
        # Create PLINK2 region file with buffer
        awk -F'\t' -v buffer="$BUFFER_KB" -v chr="$CHR" '
        BEGIN {
            buffer_bp = buffer * 1000
        }
        {
            gene = $1
            gene_id = $2
            start = $4 - buffer_bp
            end = $5 + buffer_bp
            
            # Ensure start is not negative
            if (start < 1) start = 1
            
            # PLINK2 format: chr:start-end
            region = chr ":" start "-" end
            
            # Output: region gene_symbol gene_id
            print region "\t" gene "\t" gene_id
        }
        ' "$OUTPUT_DIR/chr${CHR}_matched_temp.txt" > "$OUTPUT_DIR/chr${CHR}_regions.txt"
        
        # Create simple gene list
        cut -f1 "$OUTPUT_DIR/chr${CHR}_matched_temp.txt" > "$OUTPUT_DIR/chr${CHR}_gene_list.txt"
        
        REGION_COUNT=$(wc -l < "$OUTPUT_DIR/chr${CHR}_regions.txt")
        
        echo "   Matched: $MATCHED genes → $REGION_COUNT regions"
        printf "%-5s %-12s %-12s %-12s\n" "$CHR" "$MATCHED" "$REGION_COUNT" "Success" >> "$SUMMARY_FILE"
    else
        echo "   Matched: 0 genes (skipped)"
        printf "%-5s %-12s %-12s %-12s\n" "$CHR" "0" "0" "No matches" >> "$SUMMARY_FILE"
    fi
    
    # Cleanup temp file
    rm -f "$OUTPUT_DIR/chr${CHR}_matched_temp.txt"
done

echo ""

#==========================================
# Identify Missing Genes
#==========================================

echo "Step 3: Identifying missing genes..."

# Extract all matched gene symbols
if [ -s "$OUTPUT_DIR/matched_genes.txt" ]; then
    cut -f1 "$OUTPUT_DIR/matched_genes.txt" | sort -u > "$OUTPUT_DIR/all_matched_genes.txt"
    
    # Find genes in group file but not matched
    comm -23 \
        <(sort "$OUTPUT_DIR/genes_in_groups.txt") \
        <(sort "$OUTPUT_DIR/all_matched_genes.txt") \
        > "$OUTPUT_DIR/missing_genes.txt"
else
    # No matches at all
    cp "$OUTPUT_DIR/genes_in_groups.txt" "$OUTPUT_DIR/missing_genes.txt"
fi

TOTAL_MISSING=$(wc -l < "$OUTPUT_DIR/missing_genes.txt")

echo "  Missing genes: $TOTAL_MISSING"
echo ""

#==========================================
# Force Regenerate Missing Genes
#==========================================

TOTAL_REGENERATED=0

if [ "$FORCE_REGEN" = "yes" ] && [ "$TOTAL_MISSING" -gt 0 ]; then
    echo "Step 4: Force regenerating missing genes from variant positions..."
    echo ""
    
    # Calculate total buffer (user buffer + 10kb safety)
    SOFT_BUFFER_KB=$((BUFFER_KB + 10))
    SOFT_BUFFER_BP=$((SOFT_BUFFER_KB * 1000))
    
    echo "  Total buffer for regeneration: ${SOFT_BUFFER_KB}kb (${BUFFER_KB}kb user + 10kb safety)"
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
            missing[$1] = 1
        }
        close(missing_file)
        regen_count = 0
    }
    
    {
        # Process variant lines (odd line numbers)
        if (NR % 2 == 1) {
            gene = $1
            
            # Only process missing genes
            if (gene in missing) {
                delete positions
                delete chrs
                n_vars = 0
                
                # Extract variant positions (skip first two fields: gene name and "var")
                for (i = 3; i <= NF; i++) {
                    variant = $i
                    
                    # Parse chr:pos:ref:alt
                    split(variant, parts, ":")
                    if (length(parts) >= 2) {
                        chr = parts[1]
                        pos = parts[2]
                        
                        # Remove "chr" prefix if present
                        gsub(/^chr/, "", chr)
                        
                        # Store position and chromosome
                        positions[n_vars] = pos
                        chrs[n_vars] = chr
                        n_vars++
                    }
                }
                
                # If we have variants, calculate region
                if (n_vars > 0) {
                    # Find min and max positions
                    min_pos = positions[0]
                    max_pos = positions[0]
                    chr = chrs[0]
                    
                    for (j = 1; j < n_vars; j++) {
                        if (positions[j] < min_pos) min_pos = positions[j]
                        if (positions[j] > max_pos) max_pos = positions[j]
                    }
                    
                    # Apply buffer
                    start = min_pos - buffer
                    end = max_pos + buffer
                    
                    if (start < 1) start = 1
                    
                    # Output regenerated gene info (chr WITHOUT prefix)
                    print gene "\tREGENERATED\t" chr "\t" start "\t" end "\tNA\tNA\t" n_vars >> output
                    
                    # Output region for PLINK2 (chr WITHOUT prefix)
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
	    
    TOTAL_REGENERATED=$(wc -l < "$OUTPUT_DIR/regenerated_genes.txt")
    
    echo ""
    
    if [ "$TOTAL_REGENERATED" -gt 0 ]; then
        
        if [ "$MERGE_REGEN" = "yes" ]; then
            #==========================================
            # MERGE MODE: Add regenerated genes to main files
            #==========================================
            echo "  Merging regenerated genes into main output files..."
    
            # Merge regenerated regions into chromosome-specific files
            while IFS=$'\t' read -r region gene source; do
                # Extract chromosome from region (format: chr:start-end)
                # Remove 'chr' prefix and get just the number/letter
                CHR=$(echo "$region" | sed 's/^chr//; s/:.*//')
        
                # Append to chromosome region file (create if doesn't exist)
                echo -e "${region}\t${gene}\t${source}" >> "$OUTPUT_DIR/chr${CHR}_regions_temp.txt"
        
            done < "$OUTPUT_DIR/regenerated_regions_temp.txt"
    
            # Merge with existing files and sort by coordinates
            for chr in {1..22} X Y; do
                MAIN_REGION="$OUTPUT_DIR/chr${chr}_regions.txt"
                TEMP_REGION="$OUTPUT_DIR/chr${chr}_regions_temp.txt"
                MAIN_GENES="$OUTPUT_DIR/chr${chr}_gene_list.txt"
        
                # Merge and sort region files by start coordinate
                if [ -f "$MAIN_REGION" ] || [ -f "$TEMP_REGION" ]; then
                    {
                        [ -f "$MAIN_REGION" ] && cat "$MAIN_REGION"
                        [ -f "$TEMP_REGION" ] && cat "$TEMP_REGION"
                    } | awk -F'\t' '
                    {
                        # Extract start position from region (format: chr:start-end)
                        region = $1
                        split(region, parts, ":")
                        split(parts[2], coords, "-")
                        start = coords[1]
                
                        # Store line with start position for sorting
                        print start "\t" $0
                    }
                    ' | sort -k1,1n | cut -f2- | awk '!seen[$0]++' > "$MAIN_REGION.sorted"
            
                    mv "$MAIN_REGION.sorted" "$MAIN_REGION"
                    rm -f "$TEMP_REGION"
            
                    # Regenerate gene list from sorted regions (preserves genomic order)
                    cut -f2 "$MAIN_REGION" > "$MAIN_GENES"
                fi
            done
    
            # Add regenerated genes to matched genes
            cat "$OUTPUT_DIR/regenerated_genes.txt" >> "$OUTPUT_DIR/matched_genes.txt"
    
            # Update totals
            TOTAL_MATCHED=$((TOTAL_MATCHED + TOTAL_REGENERATED))
            TOTAL_MISSING=$((TOTAL_MISSING - TOTAL_REGENERATED))
    
            echo "  ✓ Successfully merged $TOTAL_REGENERATED genes into main files"
            echo ""
            
        else
            #==========================================
            # SEPARATE MODE: Keep regenerated genes in separate files
            #==========================================
            echo "  Saving regenerated genes to separate *_recovered.txt files..."
    
            # Create separate chromosome-specific files for regenerated genes
            while IFS=$'\t' read -r region gene source; do
                # Extract chromosome from region (format: chr:start-end)
                # Remove 'chr' prefix and get just the number/letter
                CHR=$(echo "$region" | sed 's/^chr//; s/:.*//')
        
                # Append to chromosome RECOVERED region file
                echo -e "${region}\t${gene}\t${source}" >> "$OUTPUT_DIR/chr${CHR}_regions_recovered.txt"
        
            done < "$OUTPUT_DIR/regenerated_regions_temp.txt"
    
            # Sort recovered chromosome files by coordinates
            for chr in {1..22} X Y; do
                RECOVERED_REGION="$OUTPUT_DIR/chr${chr}_regions_recovered.txt"
                RECOVERED_GENES="$OUTPUT_DIR/chr${chr}_gene_list_recovered.txt"
        
                if [ -f "$RECOVERED_REGION" ]; then
                    # Sort by start coordinate
                    awk -F'\t' '
                    {
                        # Extract start position from region (format: chr:start-end)
                        region = $1
                        split(region, parts, ":")
                        split(parts[2], coords, "-")
                        start = coords[1]
                
                        # Store line with start position for sorting
                        print start "\t" $0
                    }
                    ' "$RECOVERED_REGION" | sort -k1,1n | cut -f2- | awk '!seen[$0]++' > "$RECOVERED_REGION.sorted"
            
                    mv "$RECOVERED_REGION.sorted" "$RECOVERED_REGION"
            
                    # Regenerate gene list from sorted regions (preserves genomic order)
                    cut -f2 "$RECOVERED_REGION" > "$RECOVERED_GENES"
                fi
            done
    
            echo "  ✓ Successfully saved $TOTAL_REGENERATED genes to separate files"
            echo ""
        fi
        
        # Update summary file with regenerated info
        {
            echo ""
            echo "Regenerated Genes:"
            echo "  Total regenerated: $TOTAL_REGENERATED"
            echo "  Buffer applied: ${SOFT_BUFFER_KB}kb (${BUFFER_KB}kb + 10kb safety)"
            echo "  Output mode: $([ "$MERGE_REGEN" = "yes" ] && echo "Merged with matched genes" || echo "Separate *_recovered.txt files")"
            echo ""
        } >> "$SUMMARY_FILE"
        
        # Cleanup temp file
        rm -f "$OUTPUT_DIR/regenerated_regions_temp.txt"
    else
        echo "  No genes could be regenerated from variant positions"
        echo ""
    fi
fi

#==========================================
# Update Missing Genes List After Regeneration
#==========================================

if [ "$FORCE_REGEN" = "yes" ] && [ "$TOTAL_REGENERATED" -gt 0 ] && [ "$MERGE_REGEN" = "yes" ]; then
    echo "Step 5: Updating missing genes list after regeneration..."
    
    # Re-extract all matched gene symbols (including regenerated)
    cut -f1 "$OUTPUT_DIR/matched_genes.txt" | sort -u > "$OUTPUT_DIR/all_matched_genes.txt"
    
    # Find genes still missing
    comm -23 \
        <(sort "$OUTPUT_DIR/genes_in_groups.txt") \
        <(sort "$OUTPUT_DIR/all_matched_genes.txt") \
        > "$OUTPUT_DIR/missing_genes.txt"
    
    TOTAL_MISSING=$(wc -l < "$OUTPUT_DIR/missing_genes.txt")
    
    echo "  Updated missing genes: $TOTAL_MISSING"
    echo ""
fi

#==========================================
# Generate Summary Report
#==========================================

{
    echo ""
    echo "=========================================="
    echo "Overall Summary"
    echo "=========================================="
    echo "Genes in group file:  $TOTAL_IN_GROUPS"
    
    if [ "$MERGE_REGEN" = "yes" ]; then
        echo "Genes matched:        $TOTAL_MATCHED"
        if [ "$FORCE_REGEN" = "yes" ] && [ "$TOTAL_REGENERATED" -gt 0 ]; then
            echo "  - From coordinates: $((TOTAL_MATCHED - TOTAL_REGENERATED))"
            echo "  - Regenerated:      $TOTAL_REGENERATED"
        fi
        echo "Genes missing:        $TOTAL_MISSING"
        echo "Match rate:           $(awk -v m="$TOTAL_MATCHED" -v t="$TOTAL_IN_GROUPS" 'BEGIN {printf "%.2f%%", (m/t)*100}')"
    else
        # In separate mode, show different breakdown
        MATCHED_FROM_COORDS=$TOTAL_MATCHED
        echo "Genes from coordinates: $MATCHED_FROM_COORDS"
        if [ "$FORCE_REGEN" = "yes" ] && [ "$TOTAL_REGENERATED" -gt 0 ]; then
            echo "Genes regenerated:      $TOTAL_REGENERATED (saved separately)"
            echo "Total recovered:        $((MATCHED_FROM_COORDS + TOTAL_REGENERATED))"
        fi
        echo "Genes still missing:    $TOTAL_MISSING"
        TOTAL_RECOVERED=$((MATCHED_FROM_COORDS + TOTAL_REGENERATED))
        echo "Recovery rate:          $(awk -v m="$TOTAL_RECOVERED" -v t="$TOTAL_IN_GROUPS" 'BEGIN {printf "%.2f%%", (m/t)*100}')"
    fi
    echo ""
} >> "$SUMMARY_FILE"

#==========================================
# Final Report
#==========================================

echo "=========================================="
echo "Gene Matching Complete"
echo "=========================================="
echo "Genes in group file:  $TOTAL_IN_GROUPS"

if [ "$MERGE_REGEN" = "yes" ]; then
    echo "Genes matched:        $TOTAL_MATCHED ($(awk -v m="$TOTAL_MATCHED" -v t="$TOTAL_IN_GROUPS" 'BEGIN {printf "%.1f%%", (m/t)*100}'))"
    if [ "$FORCE_REGEN" = "yes" ] && [ "$TOTAL_REGENERATED" -gt 0 ]; then
        echo "  - From coordinates: $((TOTAL_MATCHED - TOTAL_REGENERATED))"
        echo "  - Regenerated:      $TOTAL_REGENERATED"
    fi
    echo "Genes missing:        $TOTAL_MISSING"
else
    MATCHED_FROM_COORDS=$TOTAL_MATCHED
    echo "Genes from coordinates: $MATCHED_FROM_COORDS"
    if [ "$FORCE_REGEN" = "yes" ] && [ "$TOTAL_REGENERATED" -gt 0 ]; then
        echo "Genes regenerated:      $TOTAL_REGENERATED (in *_recovered.txt files)"
        TOTAL_RECOVERED=$((MATCHED_FROM_COORDS + TOTAL_REGENERATED))
        echo "Total recovered:        $TOTAL_RECOVERED ($(awk -v m="$TOTAL_RECOVERED" -v t="$TOTAL_IN_GROUPS" 'BEGIN {printf "%.1f%%", (m/t)*100}'))"
    fi
    echo "Genes still missing:    $TOTAL_MISSING"
fi
echo ""

echo "Output Files Generated:"
echo "----------------------------------------"

# List region files
if [ "$MERGE_REGEN" = "yes" ]; then
    REGION_FILES=$(ls -1 "$OUTPUT_DIR"/chr*_regions.txt 2>/dev/null | wc -l)
    if [ "$REGION_FILES" -gt 0 ]; then
        echo "  Region files (PLINK2 format):"
        ls -1 "$OUTPUT_DIR"/chr*_regions.txt 2>/dev/null | while read f; do
            count=$(wc -l < "$f")
            printf "    %-35s %6s regions\n" "$(basename $f)" "$count"
        done
        echo ""
    fi
else
    # Show both matched and recovered files
    REGION_FILES=$(ls -1 "$OUTPUT_DIR"/chr*_regions.txt 2>/dev/null | wc -l)
    if [ "$REGION_FILES" -gt 0 ]; then
        echo "  Region files (matched genes):"
        ls -1 "$OUTPUT_DIR"/chr*_regions.txt 2>/dev/null | while read f; do
            count=$(wc -l < "$f")
            printf "    %-35s %6s regions\n" "$(basename $f)" "$count"
        done
        echo ""
    fi
    
    RECOVERED_FILES=$(ls -1 "$OUTPUT_DIR"/chr*_regions_recovered.txt 2>/dev/null | wc -l)
    if [ "$RECOVERED_FILES" -gt 0 ]; then
        echo "  Region files (regenerated genes):"
        ls -1 "$OUTPUT_DIR"/chr*_regions_recovered.txt 2>/dev/null | while read f; do
            count=$(wc -l < "$f")
            printf "    %-35s %6s regions\n" "$(basename $f)" "$count"
        done
        echo ""
    fi
fi

# List other outputs
echo "  Summary files:"

if [ "$MERGE_REGEN" = "yes" ]; then
    printf "    %-35s %6s genes\n" "matched_genes.txt" "$TOTAL_MATCHED"
else
    printf "    %-35s %6s genes\n" "matched_genes.txt" "$MATCHED_FROM_COORDS"
fi

printf "    %-35s %6s genes\n" "missing_genes.txt" "$TOTAL_MISSING"

if [ "$FORCE_REGEN" = "yes" ] && [ -f "$OUTPUT_DIR/regenerated_genes.txt" ]; then
    printf "    %-35s %6s genes\n" "regenerated_genes.txt" "$TOTAL_REGENERATED"
fi

printf "    %-35s\n" "matching_summary.txt"
printf "    %-35s\n" "matching.log"
echo ""

#==========================================
# Missing Genes Report
#==========================================

if [ "$TOTAL_MISSING" -gt 0 ]; then
    echo "⚠ WARNING: $TOTAL_MISSING genes could not be matched"
    echo ""
    echo "  Possible reasons:"
    echo "    - Gene symbols don't match Ensembl nomenclature"
    echo "    - Genes are on non-autosomal chromosomes (X, Y, MT)"
    echo "    - Genes are not in the annotation release used"
    echo "    - Typos or outdated gene symbols"
    if [ "$FORCE_REGEN" = "yes" ]; then
        echo "    - No variants found for these genes in the group file"
    fi
    echo ""
    echo "  Missing genes saved to:"
    echo "    $OUTPUT_DIR/missing_genes.txt"
    echo ""
    
    if [ "$TOTAL_MISSING" -le 20 ]; then
        echo "  All missing genes:"
        cat "$OUTPUT_DIR/missing_genes.txt" | sed 's/^/    /'
    else
        echo "  First 20 missing genes:"
        head -20 "$OUTPUT_DIR/missing_genes.txt" | sed 's/^/    /'
        echo "    ... and $((TOTAL_MISSING - 20)) more"
    fi
    echo ""
    
    if [ "$FORCE_REGEN" = "no" ]; then
        echo "  TIP: Try running with force_regen=yes to regenerate missing genes"
        echo "     from variant positions:"
        echo ""
        echo "     # Merge regenerated genes with matched:"
        echo "     ./step5_match_genes_to_groups.sh \\"
        echo "       $GENE_COORDS_DIR \\"
        echo "       $GROUP_FILE \\"
        echo "       $OUTPUT_DIR \\"
        echo "       $BUFFER_KB \\"
        echo "       yes \\"
        echo "       yes"
        echo ""
        echo "     # Keep regenerated genes separate:"
        echo "     ./step5_match_genes_to_groups.sh \\"
        echo "       $GENE_COORDS_DIR \\"
        echo "       $GROUP_FILE \\"
        echo "       $OUTPUT_DIR \\"
        echo "       $BUFFER_KB \\"
        echo "       yes \\"
        echo "       no"
        echo ""
    fi
fi

#==========================================
# Regenerated Genes Report
#==========================================

if [ "$FORCE_REGEN" = "yes" ] && [ "$TOTAL_REGENERATED" -gt 0 ]; then
    echo "=========================================="
    echo "Regenerated Genes Summary"
    echo "=========================================="
    echo "Total regenerated: $TOTAL_REGENERATED genes"
    echo "Buffer applied:    ${SOFT_BUFFER_KB}kb (${BUFFER_KB}kb user + 10kb safety)"
    echo "Output mode:       $([ "$MERGE_REGEN" = "yes" ] && echo "Merged with matched genes" || echo "Separate *_recovered.txt files")"
    echo ""
    
    if [ "$TOTAL_REGENERATED" -le 20 ]; then
        echo "All regenerated genes:"
        cut -f1 "$OUTPUT_DIR/regenerated_genes.txt" | sed 's/^/  /'
    else
        echo "First 20 regenerated genes:"
        cut -f1 "$OUTPUT_DIR/regenerated_genes.txt" | head -20 | sed 's/^/  /'
        echo "  ... and $((TOTAL_REGENERATED - 20)) more"
    fi
    echo ""
    echo "Details saved to:"
    echo "  $OUTPUT_DIR/regenerated_genes.txt"
    echo ""
    
    if [ "$MERGE_REGEN" = "no" ]; then
        echo "Regenerated gene files by chromosome:"
        ls -1 "$OUTPUT_DIR"/chr*_regions_recovered.txt 2>/dev/null | while read f; do
            count=$(wc -l < "$f")
            printf "  %-40s %6s genes\n" "$(basename $f)" "$count"
        done
        echo ""
    fi
fi

#==========================================
# Sample Output Preview
#==========================================

echo "=========================================="
echo "Sample Region File (chr1, first 5 lines):"
echo "----------------------------------------"
if [ -f "$OUTPUT_DIR/chr1_regions.txt" ]; then
    head -5 "$OUTPUT_DIR/chr1_regions.txt" | awk '{printf "  %-35s %-15s %s\n", $1, $2, $3}'
    echo ""
else
    echo "  (No chr1 regions generated)"
    echo ""
fi

if [ "$MERGE_REGEN" = "no" ] && [ -f "$OUTPUT_DIR/chr1_regions_recovered.txt" ]; then
    echo "Sample Recovered Region File (chr1, first 5 lines):"
    echo "----------------------------------------"
    head -5 "$OUTPUT_DIR/chr1_regions_recovered.txt" | awk '{printf "  %-35s %-15s %s\n", $1, $2, $3}'
    echo ""
fi

echo "=========================================="
echo "Completed: $(date)"
echo "=========================================="

#==========================================
# Cleanup Temporary Files
#==========================================

rm -f "$OUTPUT_DIR/genes_in_groups.txt"
rm -f "$OUTPUT_DIR/all_matched_genes.txt"

exit 0

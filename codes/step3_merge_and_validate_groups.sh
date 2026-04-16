#!/bin/bash

# Script: merge_and_validate_groups.sh
# Description: Merge all chromosome group files in order
# Usage: merge_and_validate_groups.sh output_file [input_dir] [pattern]

set -euo pipefail

#==========================================
# ARGUMENTS
#==========================================

if [ $# -lt 1 ]; then
    cat << 'EOF'
Usage: merge_and_validate_groups.sh <output_file> [input_dir] [pattern]

Description:
  Merge all chromosome group files in chromosomal order (1-22 or 01-22).
  No validation - simple concatenation of files.

Arguments:
  output_file    Final merged output file
  input_dir      Directory containing group files (default: current directory)
  pattern        File pattern (default: auto-detect chr*_groups.txt)

Pattern Options:
  - Auto-detect: script will find chr1-chr22 or chr01-chr22 automatically
  - Custom: specify your own pattern (e.g., "gene_groups_chr*.txt")

Examples:
  # Auto-detect in current directory
  ./merge_and_validate_groups.sh all_genes_groups.txt

  # Specify input directory
  ./merge_and_validate_groups.sh all_genes.txt ./results

  # Specify directory and custom pattern
  ./merge_and_validate_groups.sh output.txt ./data "custom_chr*_groups.txt"

  # Use full paths
  ./merge_and_validate_groups.sh /path/to/output.txt /path/to/input

Output:
  - All chromosome files merged in order
  - Summary report printed to stderr
EOF
    exit 1
fi

OUTPUT=$1
INPUT_DIR=${2:-.}
PATTERN=${3:-""}

#==========================================
# STEP 1: FIND INPUT FILES
#==========================================

echo "[$(date)] Merging chromosome group files..." >&2
echo "  Input directory: $INPUT_DIR" >&2

# Change to input directory
if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: Input directory not found: $INPUT_DIR" >&2
    exit 1
fi

cd "$INPUT_DIR"

# Auto-detect pattern if not specified
if [ -z "$PATTERN" ]; then
    # Try to detect chr01-chr22 pattern first
    if ls chr01_groups.txt >/dev/null 2>&1; then
        echo "  Detected zero-padded pattern (chr01-chr22)" >&2
        PATTERN="zero-padded"
    # Then try chr1-chr22 pattern
    elif ls chr1_groups.txt >/dev/null 2>&1; then
        echo "  Detected non-padded pattern (chr1-chr22)" >&2
        PATTERN="non-padded"
    # Fallback to wildcard
    else
        echo "  Using wildcard pattern (chr*_groups.txt)" >&2
        PATTERN="chr*_groups.txt"
    fi
else
    echo "  Using custom pattern: $PATTERN" >&2
fi

# Build file list based on pattern
declare -a INPUT_FILES

if [ "$PATTERN" = "zero-padded" ]; then
    # chr01, chr02, ..., chr22
    for i in {01..22}; do
        file="chr${i}_groups.txt"
        if [ -f "$file" ]; then
            INPUT_FILES+=("$file")
        else
            echo "  Warning: $file not found" >&2
        fi
    done
elif [ "$PATTERN" = "non-padded" ]; then
    # chr1, chr2, ..., chr22
    for i in {1..22}; do
        file="chr${i}_groups.txt"
        if [ -f "$file" ]; then
            INPUT_FILES+=("$file")
        else
            echo "  Warning: $file not found" >&2
        fi
    done
else
    # Use provided pattern
    mapfile -t INPUT_FILES < <(ls -1 $PATTERN 2>/dev/null | sort -V)
fi

# Check if any files found
if [ ${#INPUT_FILES[@]} -eq 0 ]; then
    echo "Error: No files found matching pattern in: $INPUT_DIR" >&2
    echo "  Pattern: $PATTERN" >&2
    echo "  Available files:" >&2
    ls -1 *_groups.txt 2>/dev/null | head -10 | sed 's/^/    /' >&2
    exit 1
fi

echo "  Found ${#INPUT_FILES[@]} files in order:" >&2
for f in "${INPUT_FILES[@]}"; do
    lines=$(wc -l < "$f" 2>/dev/null || echo "0")
    genes=$((lines / 2))
    echo "    - $f ($genes genes, $lines lines)" >&2
done
echo "" >&2

# Return to original directory for output
cd - >/dev/null

#==========================================
# STEP 2: MERGE FILES
#==========================================

echo "[$(date)] Merging files..." >&2

# Build full paths for merging
FULL_PATH_FILES=()
for f in "${INPUT_FILES[@]}"; do
    FULL_PATH_FILES+=("${INPUT_DIR}/$f")
done

# Simple concatenation - merge all files
cat "${FULL_PATH_FILES[@]}" > "$OUTPUT"

#==========================================
# STEP 3: SUMMARY
#==========================================

if [ -f "$OUTPUT" ]; then
    TOTAL_LINES=$(wc -l < "$OUTPUT")
    TOTAL_GENES=$((TOTAL_LINES / 2))
    FILE_SIZE=$(du -h "$OUTPUT" | cut -f1)
    
    echo "" >&2
    echo "==========================================" >&2
    echo "Merge Complete" >&2
    echo "==========================================" >&2
    echo "Output file: $OUTPUT" >&2
    echo "Total lines: $TOTAL_LINES" >&2
    echo "Total genes: $TOTAL_GENES" >&2
    echo "File size: $FILE_SIZE" >&2
    echo "" >&2
    
    # Count unique genes
    UNIQUE_GENES=$(awk '{print $1}' "$OUTPUT" | sort -u | wc -l)
    echo "Unique gene symbols: $UNIQUE_GENES" >&2
    
    # Annotation breakdown (for space-delimited format)
    if [ $TOTAL_GENES -gt 0 ]; then
        echo "" >&2
        echo "Annotation breakdown:" >&2
        awk '$2 == "anno" {for(i=3;i<=NF;i++) print $i}' "$OUTPUT" | \
            sort | uniq -c | sort -rn | head -10 | sed 's/^/  /' >&2
    fi
    
    echo "" >&2
    echo "Sample output (first 4 lines):" >&2
    head -4 "$OUTPUT" | sed 's/^/  /' >&2
    
    echo "==========================================" >&2
else
    echo "Error: Output file was not created!" >&2
    exit 1
fi

echo "" >&2
echo "[$(date)] Done!" >&2

exit 0

#!/bin/bash

# Script: create_gene_groups.sh
# Description: Convert VEP annotations to SAIGE gene group format with quality filters and deduplication
# Usage: create_gene_groups.sh input_annotations.txt output_groups.txt [annotation_filter] [priority]

set -euo pipefail

#==========================================
# USAGE
#==========================================
a
if [ $# -lt 2 ]; then
    cat << 'EOF'
Usage: create_gene_groups.sh <input_annotations.txt> <output_groups.txt> [annotation_filter] [priority]

Description:
  Convert VEP annotation file to SAIGE gene group format.
  Output uses SPACE-delimited format for var and anno rows.
  
  Handling rules:
  1. Variants with missing gene symbols → REMOVE
  2. Exact duplicate variants (same annotation) → keep one
  3. Duplicate variants (different annotations) → keep highest priority or flag conflict
  4. LoF variants where ALL annotations have NON_CAN_SPLICE → REMOVE
  5. Variants with NON_CAN_SPLICE_Count > Total_Anno → REMOVE

Arguments:
  input_annotations.txt  Input VEP annotation file (tab-delimited)
  output_groups.txt      Output gene group file (space-delimited)
  annotation_filter      Which annotations to include (default: all)
                         Options: all, lof, missense, synonymous, lof+missense
  priority               Priority order for conflict resolution (default: lof,missense,synonymous)
                         Format: comma-separated, highest to lowest
                         Special value: "keepall" - keep all annotations, no priority resolution
                         Examples: 
                           lof,missense,synonymous (default)
                           missense,lof,synonymous
                           synonymous,missense,lof
                           lof,synonymous,missense
                           keepall (keep all, create separate entries)

Priority Resolution:
  When same variant has different annotations in same gene:
  - Higher priority annotation is kept
  - Lower priority is discarded
  - Same priority → conflict (both removed, logged to .conflict file)
  - "keepall" → keep all annotations as separate entries

Input format:
  Tab-delimited file with columns:
  Variant Gene Consequence LoF LoF_flags Group Total_Anno NON_CAN_SPLICE_Count

Output format (space-delimited):
  GENE1 var variant1 variant2 variant3 ...
  GENE1 anno group1 group2 group3 ...

Conflict file:
  output_groups.txt.conflict - Variants with conflicting annotations (same priority)

Examples:
  # All annotations, default priority (lof > missense > synonymous)
  ./create_gene_groups.sh chr01_clean.txt chr01_groups.txt

  # Only LoF variants
  ./create_gene_groups.sh chr01_clean.txt chr01_lof.txt lof

  # All annotations, custom priority (missense highest)
  ./create_gene_groups.sh chr01_clean.txt chr01_groups.txt all missense,lof,synonymous

  # LoF + missense, prioritize LoF
  ./create_gene_groups.sh chr01_clean.txt chr01_coding.txt lof+missense lof,missense

  # Keep all annotations, no priority resolution
  ./create_gene_groups.sh chr01_clean.txt chr01_groups.txt all keepall

  # Keep all LoF annotations
  ./create_gene_groups.sh chr01_clean.txt chr01_lof.txt lof keepall

  # Batch process chromosomes
  for chr in {01..22}; do
    ./create_gene_groups.sh /data/data1/chr$${chr}.txt output_chr$${chr}.txt
  done
EOF
    exit 1
fi

#==========================================
# ARGUMENTS
#==========================================

INPUT="$1"
OUTPUT="$2"
ANNO_FILTER="${3:-all}"
PRIORITY_STR="${4:-lof,missense,synonymous}"

CONFLICT_FILE="${OUTPUT}.conflict"

#==========================================
# VALIDATION
#==========================================

[ ! -f "$INPUT" ] && { echo "Error: Input file '$INPUT' not found!" >&2; exit 1; }
[ ! -r "$INPUT" ] && { echo "Error: Input file '$INPUT' not readable!" >&2; exit 1; }

# Validate annotation filter
case "$ANNO_FILTER" in
    all|lof|missense|synonymous|lof+missense)
        ;;
    *)
        echo "Error: Invalid annotation filter '$ANNO_FILTER'" >&2
        echo "Valid options: all, lof, missense, synonymous, lof+missense" >&2
        exit 1
        ;;
esac

# Validate priority string
if [ "$PRIORITY_STR" != "keepall" ]; then
    if ! echo "$PRIORITY_STR" | grep -qE '^(lof|missense|synonymous)(,(lof|missense|synonymous))*$'; then
        echo "Error: Invalid priority string '$PRIORITY_STR'" >&2
        echo "Format: comma-separated list of lof, missense, synonymous OR 'keepall'" >&2
        echo "Example: lof,missense,synonymous" >&2
        echo "Example: keepall" >&2
        exit 1
    fi
fi

#==========================================
# PROCESSING
#==========================================

echo "[$(date)] Creating gene groups..." >&2
echo "  Input: $INPUT" >&2
echo "  Output: $OUTPUT" >&2
echo "  Format: SPACE-delimited" >&2
echo "  Annotation filter: $ANNO_FILTER" >&2

if [ "$PRIORITY_STR" = "keepall" ]; then
    echo "  Priority mode: KEEP ALL (no priority resolution)" >&2
else
    echo "  Priority order: $PRIORITY_STR (highest to lowest)" >&2
    echo "  Conflict file: $CONFLICT_FILE" >&2
fi

echo "  Filters:" >&2
echo "    - Variants with missing gene symbols" >&2
echo "    - LoF with all annotations NON_CAN_SPLICE" >&2
echo "    - Variants with NON_CAN_SPLICE_Count > Total_Anno" >&2

if [ "$PRIORITY_STR" = "keepall" ]; then
    echo "    - Exact duplicates (same variant, same annotation)" >&2
    echo "    - Keeping all different annotations for same variant" >&2
else
    echo "    - Duplicate variants with priority resolution" >&2
fi

awk -F'\t' -v anno_filter="$ANNO_FILTER" -v conflict_file="$CONFLICT_FILE" -v priority_str="$PRIORITY_STR" '
BEGIN {
    skipped_all_non_can = 0
    skipped_count_error = 0
    skipped_missing_gene = 0
    skipped_filtered = 0
    exact_duplicates = 0
    conflicts_resolved = 0
    conflicts_unresolved = 0
    kept_all_count = 0
    
    # Check if keepall mode
    keep_all_mode = (priority_str == "keepall")
    
    if (!keep_all_mode) {
        # Parse priority string and assign priority values
        # Higher number = higher priority
        n_priorities = split(priority_str, priority_array, ",")
        
        for (i = 1; i <= n_priorities; i++) {
            anno_type = priority_array[i]
            priority[anno_type] = n_priorities - i + 1  # Reverse: first = highest
        }
        
        # Debug: print priority assignments
        print "  Priority assignments:" > "/dev/stderr"
        for (anno in priority) {
            print "    " anno " = " priority[anno] > "/dev/stderr"
        }
        print "" > "/dev/stderr"
    } else {
        print "  Mode: KEEP ALL annotations (no priority resolution)" > "/dev/stderr"
        print "" > "/dev/stderr"
    }
}

# Skip header
NR == 1 {
    next
}

{
    variant = $1
    gene = $2
    consequence = $3
    lof = $4
    lof_flags = $5
    group = $6
    total_anno = $7
    non_can_splice_count = $8
    
    # Skip if gene symbol is missing, empty, or placeholder
    if (gene == "" || gene == "." || gene == "NA" || gene == "na") {
        skipped_missing_gene++
        if (skipped_missing_gene <= 10) {
            print "  WARNING: Skipping " variant " - missing gene symbol (\"" gene "\")" > "/dev/stderr"
        }
        next
    }
    
    # Skip if NON_CAN_SPLICE_Count > Total_Anno (data quality issue)
    if (non_can_splice_count > total_anno) {
        skipped_count_error++
        if (skipped_count_error <= 10) {
            print "  WARNING: Skipping " variant " in " gene " - NON_CAN_SPLICE_Count (" non_can_splice_count ") > Total_Anno (" total_anno ")" > "/dev/stderr"
        }
        next
    }
    
    # Skip LoF variants where ALL annotations have NON_CAN_SPLICE
    if (group == "lof" && total_anno > 0 && non_can_splice_count == total_anno) {
        skipped_all_non_can++
        next
    }
    
    # Apply annotation filter
    if (anno_filter == "lof" && group != "lof") {
        skipped_filtered++
        next
    }
    if (anno_filter == "missense" && group != "missense") {
        skipped_filtered++
        next
    }
    if (anno_filter == "synonymous" && group != "synonymous") {
        skipped_filtered++
        next
    }
    if (anno_filter == "lof+missense" && group != "lof" && group != "missense") {
        skipped_filtered++
        next
    }
    
    # Track genes in order
    if (!(gene in gene_seen)) {
        genes[++gene_count] = gene
        gene_seen[gene] = 1
    }
    
    # Create unique key for variant-gene combination
    var_gene_key = gene SUBSEP variant
    
    if (keep_all_mode) {
        # KEEP ALL MODE: Create unique key including annotation type
        var_gene_anno_key = gene SUBSEP variant SUBSEP group
        
        # Check for exact duplicates only
        if (var_gene_anno_key in variant_groups) {
            # Exact duplicate - skip
            exact_duplicates++
            next
        } else {
            # Keep this variant-annotation combination
            variant_groups[var_gene_anno_key] = group
            kept_all_count++
        }
    } else {
        # PRIORITY MODE: Original behavior
        # Check for duplicates
        if (var_gene_key in variant_groups) {
            # Duplicate variant for this gene
            existing_group = variant_groups[var_gene_key]
            
            if (existing_group == group) {
                # Exact duplicate - skip
                exact_duplicates++
                next
            } else {
                # Different annotations - apply priority
                existing_priority = priority[existing_group]
                current_priority = priority[group]
                
                if (current_priority > existing_priority) {
                    # Current has higher priority - replace
                    variant_groups[var_gene_key] = group
                    conflicts_resolved++
                    
                    # Log the resolution
                    if (conflicts_resolved <= 10) {
                        print "  RESOLVED: " gene " " variant " - keeping " group " (priority " current_priority ") over " existing_group " (priority " existing_priority ")" > "/dev/stderr"
                    }
                } else if (current_priority < existing_priority) {
                    # Existing has higher priority - keep existing, skip current
                    conflicts_resolved++
                    
                    if (conflicts_resolved <= 10) {
                        print "  RESOLVED: " gene " " variant " - keeping " existing_group " (priority " existing_priority ") over " group " (priority " current_priority ")" > "/dev/stderr"
                    }
                    next
                } else {
                    # Same priority - conflict!
                    conflicts_unresolved++
                    print gene "\t" variant "\t" existing_group "\t" group "\tPRIORITY=" current_priority >> conflict_file
                    
                    # Remove both from output
                    delete variant_groups[var_gene_key]
                    conflict_variants[var_gene_key] = 1
                    next
                }
            }
        } else {
            # New variant
            variant_groups[var_gene_key] = group
        }
    }
}

END {
    # Build gene arrays
    if (keep_all_mode) {
        # KEEP ALL MODE: Process all variant-annotation combinations
        for (vga_key in variant_groups) {
            split(vga_key, parts, SUBSEP)
            gene = parts[1]
            variant = parts[2]
            group = parts[3]
            
            # Append to gene arrays (using space as separator)
            gene_vars[gene] = gene_vars[gene] (gene_vars[gene] ? " " : "") variant
            gene_annos[gene] = gene_annos[gene] (gene_annos[gene] ? " " : "") group
            
            gene_var_count[gene]++
        }
    } else {
        # PRIORITY MODE: Original behavior
        for (vg_key in variant_groups) {
            if (vg_key in conflict_variants) {
                continue  # Skip conflicts
            }
            
            split(vg_key, parts, SUBSEP)
            gene = parts[1]
            variant = parts[2]
            group = variant_groups[vg_key]
            
            # Append to gene arrays (using space as separator)
            gene_vars[gene] = gene_vars[gene] (gene_vars[gene] ? " " : "") variant
            gene_annos[gene] = gene_annos[gene] (gene_annos[gene] ? " " : "") group
            
            gene_var_count[gene]++
        }
    }
    
    # Report statistics
    print "" > "/dev/stderr"
    if (skipped_missing_gene > 0) {
        print "  Filtered: " skipped_missing_gene " variants (missing gene symbol)" > "/dev/stderr"
        if (skipped_missing_gene > 10) {
            print "  (only first 10 warnings shown)" > "/dev/stderr"
        }
    }
    if (skipped_all_non_can > 0) {
        print "  Filtered: " skipped_all_non_can " LoF annotations (all NON_CAN_SPLICE)" > "/dev/stderr"
    }
    if (skipped_count_error > 0) {
        print "  Filtered: " skipped_count_error " variants (NON_CAN_SPLICE_Count > Total_Anno)" > "/dev/stderr"
    }
    if (skipped_filtered > 0) {
        print "  Filtered: " skipped_filtered " variants (annotation filter: " anno_filter ")" > "/dev/stderr"
    }
    if (exact_duplicates > 0) {
        print "  Removed: " exact_duplicates " exact duplicates" > "/dev/stderr"
    }
    
    if (keep_all_mode) {
        print "  Kept: " kept_all_count " variant-annotation combinations (keepall mode)" > "/dev/stderr"
    } else {
        if (conflicts_resolved > 0) {
            print "  Resolved: " conflicts_resolved " conflicts using priority: " priority_str > "/dev/stderr"
            if (conflicts_resolved > 10) {
                print "  (only first 10 resolutions shown)" > "/dev/stderr"
            }
        }
        if (conflicts_unresolved > 0) {
            print "  Conflicts: " conflicts_unresolved " unresolved (same priority) → see " conflict_file > "/dev/stderr"
        }
    }
    
    # Output each gene with its variants and annotations (space-delimited)
    for (i = 1; i <= gene_count; i++) {
        gene = genes[i]
        if (gene in gene_vars && gene_vars[gene] != "") {
            print gene " var " gene_vars[gene]
            print gene " anno " gene_annos[gene]
        }
    }
}
' "$INPUT" > "$OUTPUT"

#==========================================
# SUMMARY
#==========================================

if [ -f "$OUTPUT" ]; then
    TOTAL_LINES=$(wc -l < "$OUTPUT")
    GENE_COUNT=$((TOTAL_LINES / 2))
    FILE_SIZE=$(du -h "$OUTPUT" | cut -f1)
    
    echo "" >&2
    echo "[$(date)] Complete!" >&2
    echo "  Genes processed: $GENE_COUNT" >&2
    echo "  Output lines: $TOTAL_LINES" >&2
    echo "  File size: $FILE_SIZE" >&2
    echo "  Output file: $OUTPUT" >&2
    
    # Check conflict file (only relevant in priority mode)
    if [ "$PRIORITY_STR" != "keepall" ]; then
        if [ -f "$CONFLICT_FILE" ] && [ -s "$CONFLICT_FILE" ]; then
            CONFLICT_COUNT=$(wc -l < "$CONFLICT_FILE")
            echo "  Conflict file: $CONFLICT_FILE ($CONFLICT_COUNT variants)" >&2
            echo "" >&2
            echo "  WARNING: $CONFLICT_COUNT variants have conflicting annotations with same priority" >&2
            echo "  Please review: $CONFLICT_FILE" >&2
            echo "" >&2
            echo "  First 5 conflicts:" >&2
            head -5 "$CONFLICT_FILE" 2>/dev/null | column -t -s $'\t' | sed 's/^/    /' >&2
        fi
    else
        echo "" >&2
        echo "  Mode: KEEP ALL - all variant-annotation combinations retained" >&2
        echo "  Note: Same variant may appear multiple times with different annotations" >&2
    fi
    
    # Show sample
    echo "" >&2
    echo "Sample output (first 2 genes):" >&2
    head -4 "$OUTPUT" | sed 's/^/  /' >&2
else
    echo "[$(date)] Error: Output file was not created!" >&2
    exit 1
fi

exit 0

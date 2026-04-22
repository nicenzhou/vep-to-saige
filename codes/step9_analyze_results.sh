#!/bin/bash

# step9_analyze_results.sh
# Interactive analysis script for SAIGE-GENE results

set -euo pipefail

OUTPUT_DIR="$(dirname "\$0")"

#==========================================
# Display Banner
#==========================================

echo "=========================================="
echo "SAIGE-GENE Results Analysis Tool"
echo "=========================================="
echo "Started: $(date)"
echo "Results directory: $OUTPUT_DIR"
echo ""

#==========================================
# Parse Command Line Arguments
#==========================================

if [ $# -eq 0 ]; then
    # Interactive mode
    INTERACTIVE=true
else
    # Command mode
    INTERACTIVE=false
    OPERATIONS="\$1"
fi

#==========================================
# Available Operations
#==========================================

show_menu() {
    echo "Available Operations:"
    echo "===================="
    echo ""
    echo "Data Combination:"
    echo "  mergechrom     - Merge chunked results by chromosome"
    echo "  mergeall       - Combine all chromosomes into single file"
    echo ""
    echo "Significance Filtering:"
    echo "  findsig        - Extract significant results (multiple thresholds)"
    echo "  findgws        - Extract genome-wide significant only (p < 5e-8)"
    echo "  findsug        - Extract suggestive results (p < 1e-5)"
    echo "  findnom        - Extract nominal results (p < 0.05)"
    echo ""
    echo "Gene Ranking:"
    echo "  toptop10       - Extract top 10 genes"
    echo "  top50          - Extract top 50 genes"
    echo "  top100         - Extract top 100 genes"
    echo ""
    echo "Summary Statistics:"
    echo "  chromsum       - Per-chromosome summary table"
    echo "  fullsum        - Complete summary report"
    echo ""
    echo "Visualization Data:"
    echo "  qqdata         - Generate QQ plot data"
    echo "  mandata        - Generate Manhattan plot data"
    echo "  plotdata       - Generate both QQ and Manhattan data"
    echo ""
    echo "Comprehensive Analysis:"
    echo "  standard       - Standard analysis (mergeall+findsig+top50+fullsum)"
    echo "  full           - Full analysis (all operations)"
    echo "  quick          - Quick analysis (mergeall+top50+fullsum)"
    echo ""
    echo "Usage Examples:"
    echo "  ./step9_analyze_results.sh                           # Interactive mode"
    echo "  ./step9_analyze_results.sh standard                  # Run standard analysis"
    echo "  ./step9_analyze_results.sh mergeall+findsig+top50    # Custom combination"
    echo "  ./step9_analyze_results.sh full                      # Complete analysis"
    echo ""
}

#==========================================
# Detect Result Files
#==========================================

detect_files() {
    echo "Detecting result files..."
    
    if ls "$OUTPUT_DIR"/chr*_combined_results.txt 1> /dev/null 2>&1; then
        RESULT_TYPE="combined"
        RESULT_PATTERN="chr*_combined_results.txt"
        echo "  Found combined results (chunks were merged)"
    elif ls "$OUTPUT_DIR"/chr*_chunk*_results.txt 1> /dev/null 2>&1; then
        RESULT_TYPE="chunked"
        RESULT_PATTERN="chr*_chunk*_results.txt"
        echo "  Found chunked results"
    elif ls "$OUTPUT_DIR"/chr*_all_results.txt 1> /dev/null 2>&1; then
        RESULT_TYPE="all"
        RESULT_PATTERN="chr*_all_results.txt"
        echo "  Found chromosome-level results"
    else
        echo "ERROR: No result files found in $OUTPUT_DIR" >&2
        exit 1
    fi
    
    echo ""
}

#==========================================
# Determine P-value Column
#==========================================

get_pcol() {
    local file="\$1"
    
    if [ ! -f "$file" ]; then
        echo "ERROR: File $file not found" >&2
        return 1
    fi
    
    HEADER=$$(head -1 "$$file")
    
    # Try to find p-value column
    PCOL=$$(echo "$$HEADER" | tr '\t' '\n' | grep -n -i -E '^p$$|^p\.value$$|^pvalue$$|^p_value$$|^Pvalue_Burden$' | cut -d: -f1 | head -1)
    
    if [ -z "$PCOL" ]; then
        # If not found, assume last column
        PCOL=$$(echo "$$HEADER" | tr '\t' '\n' | wc -l)
        echo "  WARNING: Could not identify p-value column, using last column ($PCOL)" >&2
    else
        PNAME=$$(echo "$$HEADER" | cut -f"$PCOL")
        echo "  P-value column: $$PCOL ($$PNAME)"
    fi
}

#==========================================
# Operation: Merge Chromosome Chunks
#==========================================

op_mergechrom() {
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "Merging chunked results by chromosome..."
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    
    if [ "$RESULT_TYPE" != "chunked" ]; then
        echo "  Skipping: Results are not chunked"
        echo ""
        return
    fi
    
    for chr in {1..22} X Y; do
        if ls "$$OUTPUT_DIR"/chr$${chr}_chunk*_results.txt 1> /dev/null 2>&1; then
            echo "  Processing chr${chr}..."
            
            FIRST_CHUNK=$$(ls "$$OUTPUT_DIR"/chr${chr}_chunk*_results.txt | sort -V | head -1)
            
            head -1 "$$FIRST_CHUNK" > "$$OUTPUT_DIR/chr${chr}_combined.txt"
            tail -n +2 -q "$$OUTPUT_DIR"/chr$${chr}_chunk*_results.txt >> "$$OUTPUT_DIR/chr$${chr}_combined.txt"
            
            COUNT=$$(tail -n +2 "$$OUTPUT_DIR/chr${chr}_combined.txt" | wc -l)
            echo "    Combined $COUNT genes"
        fi
    done
    
    RESULT_PATTERN="chr*_combined.txt"
    RESULT_TYPE="combined"
    echo "  ✓ Chromosome merging complete"
    echo ""
}

#==========================================
# Operation: Merge All Chromosomes
#==========================================

op_mergeall() {
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "Combining all chromosomes..."
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    
    # Determine pattern based on current state
    if [ -f "$OUTPUT_DIR/all_results.txt" ]; then
        echo "  all_results.txt already exists"
        TOTAL_GENES=$$(tail -n +2 "$$OUTPUT_DIR/all_results.txt" | wc -l)
        echo "  Total genes: $TOTAL_GENES"
        echo ""
        return
    fi
    
    FIRST_FILE=$$(ls "$$OUTPUT_DIR"/$RESULT_PATTERN 2>/dev/null | sort -V | head -1)
    
    if [ -z "$FIRST_FILE" ]; then
        echo "ERROR: No result files found" >&2
        exit 1
    fi
    
    head -1 "$$FIRST_FILE" > "$$OUTPUT_DIR/all_results.txt"
    tail -n +2 -q "$$OUTPUT_DIR"/$$RESULT_PATTERN >> "$OUTPUT_DIR/all_results.txt"
    
    TOTAL_GENES=$$(tail -n +2 "$$OUTPUT_DIR/all_results.txt" | wc -l)
    echo "  Total genes tested: $TOTAL_GENES"
    echo "  ✓ Saved to: all_results.txt"
    echo ""
}

#==========================================
# Operation: Extract Significant Results
#==========================================

op_findsig() {
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "Extracting significant results..."
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    
    if [ ! -f "$OUTPUT_DIR/all_results.txt" ]; then
        echo "  ERROR: all_results.txt not found. Run mergeall first."
        echo ""
        return 1
    fi
    
    get_pcol "$OUTPUT_DIR/all_results.txt"
    
    # Genome-wide significant
    head -1 "$$OUTPUT_DIR/all_results.txt" > "$$OUTPUT_DIR/genome_wide_sig.txt"
    tail -n +2 "$$OUTPUT_DIR/all_results.txt" | awk -v pcol="$$PCOL" -F'\t' '$$pcol != "NA" && $$pcol < 5e-8' >> "$OUTPUT_DIR/genome_wide_sig.txt"
    GWS_COUNT=$$(tail -n +2 "$$OUTPUT_DIR/genome_wide_sig.txt" | wc -l)
    echo "  Genome-wide (p < 5e-8):  $GWS_COUNT genes"
    
    # Suggestive
    head -1 "$$OUTPUT_DIR/all_results.txt" > "$$OUTPUT_DIR/suggestive_sig.txt"
    tail -n +2 "$$OUTPUT_DIR/all_results.txt" | awk -v pcol="$$PCOL" -F'\t' '$$pcol != "NA" && $$pcol < 1e-5' >> "$OUTPUT_DIR/suggestive_sig.txt"
    SUG_COUNT=$$(tail -n +2 "$$OUTPUT_DIR/suggestive_sig.txt" | wc -l)
    echo "  Suggestive (p < 1e-5):   $SUG_COUNT genes"
    
    # Nominal
    head -1 "$$OUTPUT_DIR/all_results.txt" > "$$OUTPUT_DIR/nominal_sig.txt"
    tail -n +2 "$$OUTPUT_DIR/all_results.txt" | awk -v pcol="$$PCOL" -F'\t' '$$pcol != "NA" && $$pcol < 0.05' >> "$OUTPUT_DIR/nominal_sig.txt"
    NOM_COUNT=$$(tail -n +2 "$$OUTPUT_DIR/nominal_sig.txt" | wc -l)
    echo "  Nominal (p < 0.05):      $NOM_COUNT genes"
    
    # P < 0.01
    head -1 "$$OUTPUT_DIR/all_results.txt" > "$$OUTPUT_DIR/sig_p001.txt"
    tail -n +2 "$$OUTPUT_DIR/all_results.txt" | awk -v pcol="$$PCOL" -F'\t' '$$pcol != "NA" && $$pcol < 0.01' >> "$OUTPUT_DIR/sig_p001.txt"
    P001_COUNT=$$(tail -n +2 "$$OUTPUT_DIR/sig_p001.txt" | wc -l)
    echo "  P < 0.01:                $P001_COUNT genes"
    
    echo "  ✓ Significance filtering complete"
    echo ""
}

#==========================================
# Operation: Extract Top Genes
#==========================================

op_topgenes() {
    local n=\$1
    local outfile="$$OUTPUT_DIR/top$${n}_genes.txt"
    
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "Extracting top $n genes..."
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    
    if [ ! -f "$OUTPUT_DIR/all_results.txt" ]; then
        echo "  ERROR: all_results.txt not found. Run mergeall first."
        echo ""
        return 1
    fi
    
    get_pcol "$OUTPUT_DIR/all_results.txt"
    
    head -1 "$$OUTPUT_DIR/all_results.txt" > "$$outfile"
    tail -n +2 "$OUTPUT_DIR/all_results.txt" | \
        awk -v pcol="$$PCOL" -F'\t' '$$pcol != "NA"' | \
        sort -t$$'\t' -k"$$PCOL","$PCOL"g | \
        head -$$n >> "$$outfile"
    
    echo "  ✓ Saved to: top${n}_genes.txt"
    
    # Display top 20
    if [ $n -ge 20 ]; then
        echo ""
        echo "  Top 20 genes:"
        if command -v column &> /dev/null; then
            head -21 "$$outfile" | column -t -s$$'\t' | head -21
        else
            head -21 "$outfile"
        fi
    fi
    
    echo ""
}

#==========================================
# Operation: Chromosome Summary
#==========================================

op_chromsum() {
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "Generating per-chromosome summary..."
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    
    {
        printf "%-5s %-12s %-12s %-12s %-12s %-12s\n" \
            "Chr" "Total_Genes" "GWS_p<5e-8" "Sug_p<1e-5" "Nom_p<0.05" "P<0.01"
        echo "------------------------------------------------------------------------"
        
        for chr in {1..22} X Y; do
            CHR_FILE=""
            if [ -f "$$OUTPUT_DIR/chr$${chr}_combined_results.txt" ]; then
                CHR_FILE="$$OUTPUT_DIR/chr$${chr}_combined_results.txt"
            elif [ -f "$$OUTPUT_DIR/chr$${chr}_combined.txt" ]; then
                CHR_FILE="$$OUTPUT_DIR/chr$${chr}_combined.txt"
            elif [ -f "$$OUTPUT_DIR/chr$${chr}_all_results.txt" ]; then
                CHR_FILE="$$OUTPUT_DIR/chr$${chr}_all_results.txt"
            fi
            
            if [ -n "$$CHR_FILE" ] && [ -f "$$CHR_FILE" ]; then
                get_pcol "$CHR_FILE" > /dev/null
                
                TOTAL=$$(tail -n +2 "$$CHR_FILE" | wc -l)
                GWS=$$(tail -n +2 "$$CHR_FILE" | awk -v pcol="$$PCOL" -F'\t' '$$pcol != "NA" && $pcol < 5e-8' | wc -l)
                SUG=$$(tail -n +2 "$$CHR_FILE" | awk -v pcol="$$PCOL" -F'\t' '$$pcol != "NA" && $pcol < 1e-5' | wc -l)
                NOM=$$(tail -n +2 "$$CHR_FILE" | awk -v pcol="$$PCOL" -F'\t' '$$pcol != "NA" && $pcol < 0.05' | wc -l)
                P001=$$(tail -n +2 "$$CHR_FILE" | awk -v pcol="$PCOL" -F'\t' '$$pcol != "NA" && $$pcol < 0.01' | wc -l)
                
                printf "%-5s %-12s %-12s %-12s %-12s %-12s\n" \
                    "$$chr" "$$TOTAL" "$$GWS" "$$SUG" "$$NOM" "$$P001"
            fi
        done
    } > "$OUTPUT_DIR/chromosome_summary.txt"
    
    cat "$OUTPUT_DIR/chromosome_summary.txt"
    echo ""
    echo "  ✓ Saved to: chromosome_summary.txt"
    echo ""
}

#==========================================
# Operation: Full Summary Report
#==========================================

op_fullsum() {
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "Generating full summary report..."
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    
    if [ ! -f "$OUTPUT_DIR/all_results.txt" ]; then
        echo "  ERROR: all_results.txt not found. Run mergeall first."
        echo ""
        return 1
    fi
    
    get_pcol "$OUTPUT_DIR/all_results.txt"
    
    TOTAL_GENES=$$(tail -n +2 "$$OUTPUT_DIR/all_results.txt" | wc -l)
    
    # Count significant results
    GWS_COUNT=$$(tail -n +2 "$$OUTPUT_DIR/all_results.txt" | awk -v pcol="$$PCOL" -F'\t' '$$pcol != "NA" && $pcol < 5e-8' | wc -l)
    SUG_COUNT=$$(tail -n +2 "$$OUTPUT_DIR/all_results.txt" | awk -v pcol="$$PCOL" -F'\t' '$$pcol != "NA" && $pcol < 1e-5' | wc -l)
    NOM_COUNT=$$(tail -n +2 "$$OUTPUT_DIR/all_results.txt" | awk -v pcol="$$PCOL" -F'\t' '$$pcol != "NA" && $pcol < 0.05' | wc -l)
    P001_COUNT=$$(tail -n +2 "$$OUTPUT_DIR/all_results.txt" | awk -v pcol="$$PCOL" -F'\t' '$$pcol != "NA" && $pcol < 0.01' | wc -l)
    
    {
        echo "=========================================="
        echo "SAIGE-GENE Results Summary"
        echo "=========================================="
        echo "Analysis Date: $(date)"
        echo "Results Directory: $OUTPUT_DIR"
        echo ""
        echo "Total Genes Tested: $TOTAL_GENES"
        echo ""
        echo "Significance Thresholds:"
        echo "----------------------------------------"
        printf "  %-25s %8s %8s\n" "Threshold" "Count" "Percent"
        echo "----------------------------------------"
        printf "  %-25s %8s %7.2f%%\n" "Genome-wide (p < 5e-8)" "$$GWS_COUNT" "$$(awk -v g="$$GWS_COUNT" -v t="$$TOTAL_GENES" 'BEGIN {printf "%.2f", (g/t)*100}')"
        printf "  %-25s %8s %7.2f%%\n" "Suggestive (p < 1e-5)" "$$SUG_COUNT" "$$(awk -v g="$$SUG_COUNT" -v t="$$TOTAL_GENES" 'BEGIN {printf "%.2f", (g/t)*100}')"
        printf "  %-25s %8s %7.2f%%\n" "P < 0.01" "$$P001_COUNT" "$$(awk -v g="$$P001_COUNT" -v t="$$TOTAL_GENES" 'BEGIN {printf "%.2f", (g/t)*100}')"
        printf "  %-25s %8s %7.2f%%\n" "Nominal (p < 0.05)" "$$NOM_COUNT" "$$(awk -v g="$$NOM_COUNT" -v t="$$TOTAL_GENES" 'BEGIN {printf "%.2f", (g/t)*100}')"
        echo ""
        echo "Output Files Generated:"
        echo "----------------------------------------"
        
        if [ -f "$OUTPUT_DIR/all_results.txt" ]; then
            echo "  ✓ all_results.txt"
        fi
        if [ -f "$OUTPUT_DIR/genome_wide_sig.txt" ]; then
            echo "  ✓ genome_wide_sig.txt"
        fi
        if [ -f "$OUTPUT_DIR/suggestive_sig.txt" ]; then
            echo "  ✓ suggestive_sig.txt"
        fi
        if [ -f "$OUTPUT_DIR/nominal_sig.txt" ]; then
            echo "  ✓ nominal_sig.txt"
        fi
        if [ -f "$OUTPUT_DIR/sig_p001.txt" ]; then
            echo "  ✓ sig_p001.txt"
        fi
        if [ -f "$OUTPUT_DIR/top50_genes.txt" ]; then
            echo "  ✓ top50_genes.txt"
        fi
        if [ -f "$OUTPUT_DIR/chromosome_summary.txt" ]; then
            echo "  ✓ chromosome_summary.txt"
        fi
        if [ -f "$OUTPUT_DIR/qq_plot_data.txt" ]; then
            echo "  ✓ qq_plot_data.txt"
        fi
        if [ -f "$OUTPUT_DIR/manhattan_plot_data.txt" ]; then
            echo "  ✓ manhattan_plot_data.txt"
        fi
        echo ""
        echo "=========================================="
    } | tee "$OUTPUT_DIR/analysis_summary.txt"
    
    echo ""
    echo "  ✓ Summary saved to: analysis_summary.txt"
    echo ""
}

#==========================================
# Operation: QQ Plot Data
#==========================================

op_qqdata() {
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "Generating QQ plot data..."
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    
    if [ ! -f "$OUTPUT_DIR/all_results.txt" ]; then
        echo "  ERROR: all_results.txt not found. Run mergeall first."
        echo ""
        return 1
    fi
    
    get_pcol "$OUTPUT_DIR/all_results.txt"
    
    TOTAL_GENES=$$(tail -n +2 "$$OUTPUT_DIR/all_results.txt" | wc -l)
    
    tail -n +2 "$OUTPUT_DIR/all_results.txt" | \
        awk -v pcol="$$PCOL" -F'\t' '$$pcol != "NA" && $$pcol > 0 {print $$pcol}' | \
        sort -g | \
        awk -v total="$TOTAL_GENES" '{
            obs = -log(\$1)/log(10)
            exp = -log((NR-0.5)/total)/log(10)
            print exp "\t" obs
        }' > "$OUTPUT_DIR/qq_plot_data.txt"
    
    POINTS=$$(wc -l < "$$OUTPUT_DIR/qq_plot_data.txt")
    echo "  Generated $POINTS data points"
    echo "  ✓ Saved to: qq_plot_data.txt"
    echo "  Columns: Expected_log10P  Observed_log10P"
    echo ""
}

#==========================================
# Operation: Manhattan Plot Data
#==========================================

op_mandata() {
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "Generating Manhattan plot data..."
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    
    if [ ! -f "$OUTPUT_DIR/all_results.txt" ]; then
        echo "  ERROR: all_results.txt not found. Run mergeall first."
        echo ""
        return 1
    fi
    
    get_pcol "$OUTPUT_DIR/all_results.txt"
    
    echo "Gene CHR Pvalue NegLog10P" | tr ' ' '\t' > "$OUTPUT_DIR/manhattan_plot_data.txt"
    
    tail -n +2 "$OUTPUT_DIR/all_results.txt" | \
        awk -v pcol="$PCOL" -F'\t' '
        $$pcol != "NA" && $$pcol > 0 {
            gene = \$1
            pval = $pcol
            
            # Try to extract chromosome from gene name or other columns
            chr = "NA"
            if (\$1 ~ /^chr/) {
                split(\$1, arr, /[_:]/)
                chr = arr[1]
                gsub("chr", "", chr)
            }
            
            # Calculate -log10(p)
            log10p = -log(pval)/log(10)
            
            print gene "\t" chr "\t" pval "\t" log10p
        }' >> "$OUTPUT_DIR/manhattan_plot_data.txt"
    
    POINTS=$$(tail -n +2 "$$OUTPUT_DIR/manhattan_plot_data.txt" | wc -l)
    echo "  Generated $POINTS data points"
    echo "  ✓ Saved to: manhattan_plot_data.txt"
    echo "  Columns: Gene  CHR  Pvalue  NegLog10P"
    echo ""
}

#==========================================
# Operation Dispatcher
#==========================================

run_operation() {
    local op=\$1
    
    case $op in
        mergechrom)
            op_mergechrom
            ;;
        mergeall)
            op_mergeall
            ;;
        findsig)
            op_findsig
            ;;
        findgws)
            echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
            echo "Extracting genome-wide significant results..."
            echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
            if [ ! -f "$OUTPUT_DIR/all_results.txt" ]; then
                echo "  ERROR: all_results.txt not found. Run mergeall first."
                echo ""
                return 1
            fi
            get_pcol "$OUTPUT_DIR/all_results.txt"
            head -1 "$$OUTPUT_DIR/all_results.txt" > "$$OUTPUT_DIR/genome_wide_sig.txt"
            tail -n +2 "$$OUTPUT_DIR/all_results.txt" | awk -v pcol="$$PCOL" -F'\t' '$$pcol != "NA" && $$pcol < 5e-8' >> "$OUTPUT_DIR/genome_wide_sig.txt"
            GWS_COUNT=$$(tail -n +2 "$$OUTPUT_DIR/genome_wide_sig.txt" | wc -l)
            echo "  Genome-wide significant: $GWS_COUNT genes"
            echo "  ✓ Saved to: genome_wide_sig.txt"
            echo ""
            ;;
        findsug)
            echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
            echo "Extracting suggestive results..."
            echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
            if [ ! -f "$OUTPUT_DIR/all_results.txt" ]; then
                echo "  ERROR: all_results.txt not found. Run mergeall first."
                echo ""
                return 1
            fi
            get_pcol "$OUTPUT_DIR/all_results.txt"
            head -1 "$$OUTPUT_DIR/all_results.txt" > "$$OUTPUT_DIR/suggestive_sig.txt"
            tail -n +2 "$$OUTPUT_DIR/all_results.txt" | awk -v pcol="$$PCOL" -F'\t' '$$pcol != "NA" && $$pcol < 1e-5' >> "$OUTPUT_DIR/suggestive_sig.txt"
            SUG_COUNT=$$(tail -n +2 "$$OUTPUT_DIR/suggestive_sig.txt" | wc -l)
            echo "  Suggestive results: $SUG_COUNT genes"
            echo "  ✓ Saved to: suggestive_sig.txt"
            echo ""
            ;;
        findnom)
            echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
            echo "Extracting nominal significant results..."
            echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
            if [ ! -f "$OUTPUT_DIR/all_results.txt" ]; then
                echo "  ERROR: all_results.txt not found. Run mergeall first."
                echo ""
                return 1
            fi
            get_pcol "$OUTPUT_DIR/all_results.txt"
            head -1 "$$OUTPUT_DIR/all_results.txt" > "$$OUTPUT_DIR/nominal_sig.txt"
            tail -n +2 "$$OUTPUT_DIR/all_results.txt" | awk -v pcol="$$PCOL" -F'\t' '$$pcol != "NA" && $$pcol < 0.05' >> "$OUTPUT_DIR/nominal_sig.txt"
            NOM_COUNT=$$(tail -n +2 "$$OUTPUT_DIR/nominal_sig.txt" | wc -l)
            echo "  Nominal significant: $NOM_COUNT genes"
            echo "  ✓ Saved to: nominal_sig.txt"
            echo ""
            ;;
        top10)
            op_topgenes 10
            ;;
        top50)
            op_topgenes 50
            ;;
        top100)
            op_topgenes 100
            ;;
        chromsum)
            op_chromsum
            ;;
        fullsum)
            op_fullsum
            ;;
        qqdata)
            op_qqdata
            ;;
        mandata)
            op_mandata
            ;;
        plotdata)
            op_qqdata
            op_mandata
            ;;
        standard)
            detect_files
            op_mergeall
            op_findsig
            op_topgenes 50
            op_fullsum
            ;;
        full)
            detect_files
            op_mergechrom
            op_mergeall
            op_findsig
            op_topgenes 10
            op_topgenes 50
            op_topgenes 100
            op_chromsum
            op_qqdata
            op_mandata
            op_fullsum
            ;;
        quick)
            detect_files
            op_mergeall
            op_topgenes 50
            op_fullsum
            ;;
        *)
            echo "ERROR: Unknown operation '$op'" >&2
            echo ""
            show_menu
            exit 1
            ;;
    esac
}

#==========================================
# Main Execution
#==========================================

if [ "$INTERACTIVE" = true ]; then
    show_menu
    echo ""
    read -p "Enter operations (e.g., standard or mergeall+findsig+top50): " OPERATIONS if [ -z "$OPERATIONS" ]; then
        echo "No operations specified. Exiting."
        exit 0
    fi
fi

# Detect files first
detect_files

# Parse and execute operations
IFS='+' read -ra OPS <<< "$OPERATIONS"

echo "Executing operations: ${OPS[*]}"
echo ""

for op in "${OPS[@]}"; do
    run_operation "$op"
done

#==========================================
# Final Summary
#==========================================

echo "=========================================="
echo "Analysis Complete"
echo "=========================================="
echo "Completed: $(date)"
echo ""

if [ -f "$OUTPUT_DIR/analysis_summary.txt" ]; then
    echo "Summary report available at:"
    echo "  $OUTPUT_DIR/analysis_summary.txt"
    echo ""
fi

echo "Generated files in $OUTPUT_DIR:"
ls -lh "$OUTPUT_DIR"/*.txt 2>/dev/null | awk '{print "  " $9 " (" $5 ")"}'
echo ""

#==========================================
# Next Steps Suggestions
#==========================================

echo "=========================================="
echo "Next Steps"
echo "=========================================="
echo ""
echo "Visualization:"
echo "  - Use qq_plot_data.txt for QQ plots"
echo "  - Use manhattan_plot_data.txt for Manhattan plots"
echo ""
echo "Further Analysis:"
echo "  - Review genome_wide_sig.txt for significant genes"
echo "  - Check top50_genes.txt for strongest associations"
echo "  - Examine chromosome_summary.txt for regional patterns"
echo ""
echo "Example R commands for plotting:"
echo ""
cat << 'RCODE'
# QQ Plot
data <- read.table("qq_plot_data.txt", header=FALSE)
png("qq_plot.png", width=1000, height=1000, res=150)
plot(data$V1, data$V2,
     xlab=expression(-log[10](Expected~P)),
     ylab=expression(-log[10](Observed~P)),
     main="QQ Plot - SAIGE-GENE Results",
     pch=20, cex=0.8, col=rgb(0,0,0,0.5))
abline(0, 1, col="red", lwd=2)
dev.off()

# Manhattan Plot (simple version)
data <- read.table("manhattan_plot_data.txt", header=TRUE)
data <- data[!is.na(data$NegLog10P),]
png("manhattan_plot.png", width=1400, height=800, res=150)
plot(1:nrow(data), data$NegLog10P,
     xlab="Gene Index", ylab=expression(-log[10](P)),
     main="Manhattan Plot - SAIGE-GENE Results",
     pch=20, cex=0.6, col=rgb(0,0,1,0.5))
abline(h=-log10(5e-8), col="red", lwd=2, lty=2)
abline(h=-log10(1e-5), col="blue", lwd=1, lty=2)
dev.off()
RCODE
echo ""

echo "Advanced filtering examples:"
echo ""
echo "# Extract loss-of-function variants only:"
echo "  grep -i 'lof' genome_wide_sig.txt > lof_significant.txt"
echo ""
echo "# Find specific gene:"
echo "  grep -i 'GENE_NAME' all_results.txt"
echo ""
echo "# Extract genes with MAC > 100:"
echo "  awk -F'\t' '\$MAC_COLUMN > 100' all_results.txt > high_mac_results.txt"
echo ""
echo "# Sort by beta (effect size) - adjust column number:"
echo "  (head -1 all_results.txt; tail -n +2 all_results.txt | sort -t\$'\\t' -k3,3gr) > sorted_by_beta.txt"
echo ""

echo "=========================================="
echo "Analysis script completed successfully"
echo "=========================================="
echo ""

exit 0

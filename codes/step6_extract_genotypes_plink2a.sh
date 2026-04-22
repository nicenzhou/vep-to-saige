#!/bin/bash

# Script: step6_extract_genotypes_plink2.sh
# Description: Extract genotypes for gene regions using PLINK2a with optional chunking
# Usage: step6_extract_genotypes_plink2.sh <input_dir> <regions_dir> <output_dir> [threads] [output_format] [input_format] [input_prefix] [chunk_size]

set -euo pipefail

INPUT_DIR="${1:-}"
REGIONS_DIR="${2:-}"
OUTPUT_DIR="${3:-gene_bfiles}"
THREADS="${4:-4}"
OUTPUT_FORMAT="${5:-bfile}"     # Options: bfile, pgen, vcf, bgen
INPUT_FORMAT="${6:-auto}"       # Options: auto, bfile, pgen, vcf, bgen
INPUT_PREFIX="${7:-chr}"        # Input file prefix pattern (default: chr)
CHUNK_SIZE="${8:-0}"            # Chunk size (0 = no chunking, extract all in one file)

if [ -z "$INPUT_DIR" ] || [ -z "$REGIONS_DIR" ]; then
    cat << 'EOF'
Usage: step6_extract_genotypes_plink2.sh <input_dir> <regions_dir> <output_dir> [threads] [output_format] [input_format] [input_prefix] [chunk_size]

Description:
  Extract genotypes for genes using PLINK2a.
  Supports chunking to split large extractions into smaller files.
  Creates files per chromosome (and per chunk if chunk_size > 0).

Arguments:
  input_dir      Directory with chromosome genotype files
  regions_dir    Directory with chr*_regions.txt files (from step5_match_genes_to_groups.sh)
  output_dir     Output directory for extracted files (default: gene_bfiles)
  threads        Number of threads for PLINK2a (default: 4)
  output_format  Output format: bfile, pgen, vcf, bgen (default: bfile)
  input_format   Input format: auto, bfile, pgen, vcf, bgen (default: auto)
  input_prefix   Input file prefix pattern (default: chr)
  chunk_size     Number of genes per chunk (default: 0 = no chunking)

Chunking Feature:
  chunk_size = 0   : Extract all genes in one file per chromosome (default)
  chunk_size = 20  : Extract 20 genes per chunk, multiple files per chromosome
  chunk_size = 50  : Extract 50 genes per chunk
  
  Output naming with chunks:
    chr1_genes_chunk1.bed/bim/fam   (genes 1-20)
    chr1_genes_chunk2.bed/bim/fam   (genes 21-40)
    chr1_genes_chunk3.bed/bim/fam   (genes 41-60)
    ...

Input File Prefix Examples:
  chr          -> chr1.bed, chr2.bed, ... (default)
  wgs_chr      -> wgs_chr1.bed, wgs_chr2.bed, ...
  chromosome   -> chromosome1.bed, chromosome2.bed, ...
  imputed_chr  -> imputed_chr1.bed, imputed_chr2.bed, ...
  ""           -> 1.bed, 2.bed, ... (empty string for no prefix)

Input Formats (auto-detected if input_format=auto):
  bfile  - PLINK1 binary (.bed/.bim/.fam)
  pgen   - PLINK2 binary (.pgen/.pvar/.psam)
  vcf    - VCF format (.vcf.gz or .vcf)
  bgen   - BGEN format (.bgen with .sample)

Output Formats:
  bfile  - PLINK1 binary format (.bed/.bim/.fam)
  pgen   - PLINK2 binary format (.pgen/.pvar/.psam)
  vcf    - VCF format (.vcf.gz with .tbi index)
  bgen   - BGEN v1.2 format (.bgen with .sample)

Output Files (no chunking, chunk_size=0):
  - chr1_genes.bed/bim/fam   All genes on chr1
  - chr2_genes.bed/bim/fam   All genes on chr2
  - extraction_summary.txt   Summary of extraction
  - extraction.log           Detailed log

Output Files (with chunking, chunk_size=20):
  - chr1_genes_chunk1.bed/bim/fam   Genes 1-20 on chr1
  - chr1_genes_chunk2.bed/bim/fam   Genes 21-40 on chr1
  - chr2_genes_chunk1.bed/bim/fam   Genes 1-20 on chr2
  - extraction_summary.txt          Summary with chunk info
  - extraction.log                  Detailed log

Examples:
  # Default: all genes in one file per chromosome
  ./step6_extract_genotypes_plink2.sh /data/genotypes plink_regions gene_bfiles

  # Extract 20 genes per chunk
  ./step6_extract_genotypes_plink2.sh /data/genotypes plink_regions gene_bfiles 16 bfile auto chr 20

  # Custom prefix with chunking
  ./step6_extract_genotypes_plink2.sh /data/genotypes plink_regions gene_bfiles 16 bfile auto wgs_chr 50

  # BGEN format with 30 genes per chunk
  ./step6_extract_genotypes_plink2.sh /data/bgen plink_regions gene_bfiles 16 bgen bgen imputed_chr 30

  # No chunking (extract all genes together)
  ./step6_extract_genotypes_plink2.sh /data/genotypes plink_regions gene_bfiles 8 bfile auto chr 0

Use Cases for Chunking:
  - Memory constraints: Smaller chunks use less RAM
  - Parallel processing: Process chunks independently
  - Testing: Extract small chunks first for QC
  - Large cohorts: Manage file sizes for very large datasets

Use in SAIGE-GENE (with chunks):
  # Process each chunk separately
  for chunk in gene_bfiles/chr1_genes_chunk*.bgen; do
    Rscript step2_SPAtests.R \
      --bgenFile=$chunk \
      --bgenFileIndex=${chunk}.bgi \
      --sampleFile=${chunk%.bgen}.sample \
      --groupFile=chr1_groups_$(basename $chunk .bgen).txt \
      ...
  done

EOF
    exit 1
fi

#==========================================
# Validate Inputs
#==========================================

[ ! -d "$INPUT_DIR" ] && { echo "ERROR: Input directory not found: $INPUT_DIR" >&2; exit 1; }
[ ! -d "$REGIONS_DIR" ] && { echo "ERROR: Regions directory not found: $REGIONS_DIR" >&2; exit 1; }

# Validate output format
OUTPUT_FORMAT=$(echo "$OUTPUT_FORMAT" | tr '[:upper:]' '[:lower:]')
if [[ ! "$OUTPUT_FORMAT" =~ ^(bfile|pgen|vcf|bgen)$ ]]; then
    echo "ERROR: Invalid output_format. Must be: bfile, pgen, vcf, or bgen" >&2
    exit 1
fi

# Validate input format
INPUT_FORMAT=$(echo "$INPUT_FORMAT" | tr '[:upper:]' '[:lower:]')
if [[ ! "$INPUT_FORMAT" =~ ^(auto|bfile|pgen|vcf|bgen)$ ]]; then
    echo "ERROR: Invalid input_format. Must be: auto, bfile, pgen, vcf, or bgen" >&2
    exit 1
fi

# Validate threads
if ! [[ "$THREADS" =~ ^[0-9]+$ ]]; then
    echo "ERROR: Threads must be a positive integer (got: $THREADS)" >&2
    exit 1
fi

# Validate chunk_size
if ! [[ "$CHUNK_SIZE" =~ ^[0-9]+$ ]]; then
    echo "ERROR: Chunk size must be a non-negative integer (got: $CHUNK_SIZE)" >&2
    exit 1
fi

#==========================================
# Load PLINK2a Module
#==========================================

echo "Loading PLINK2a module..." >&2

# Try to load plink2a module
if command -v module &> /dev/null; then
    if module avail plink2a 2>&1 | grep -q plink2a; then
        module load plink2a
        echo "  ✓ Module plink2a loaded" >&2
    elif module avail plink/2a 2>&1 | grep -q "plink/2a"; then
        module load plink/2a
        echo "  ✓ Module plink/2a loaded" >&2
    elif module avail PLINK/2a 2>&1 | grep -q "PLINK/2a"; then
        module load PLINK/2a
        echo "  ✓ Module PLINK/2a loaded" >&2
    else
        echo "  ⚠ No plink2a module found, using system plink2a" >&2
    fi
else
    echo "  ⚠ Module system not available, using system plink2a" >&2
fi

# Check for plink2a
if ! command -v plink2a &> /dev/null; then
    echo "ERROR: plink2a not found in PATH" >&2
    echo "Please install PLINK2 alpha or load the appropriate module" >&2
    echo "Download: https://www.cog-genomics.org/plink/2.0/" >&2
    exit 1
fi

# Verify plink2a version
PLINK_VERSION=$(plink2a --version 2>&1 | head -n1 || echo "Unknown")
echo "  PLINK2a version: $PLINK_VERSION" >&2
echo "" >&2

mkdir -p "$OUTPUT_DIR"

#==========================================
# Function to Auto-Detect Input Format
#==========================================

detect_input_format() {
    local input_dir="$1"
    local chr="$2"
    local prefix="$3"
    
    for pattern in "${prefix}${chr}" "${prefix}${chr}_merged"; do
        if [ -f "$input_dir/${pattern}.bed" ] && [ -f "$input_dir/${pattern}.bim" ] && [ -f "$input_dir/${pattern}.fam" ]; then
            echo "bfile:$input_dir/${pattern}"
            return 0
        fi
        
        if [ -f "$input_dir/${pattern}.pgen" ] && [ -f "$input_dir/${pattern}.pvar" ] && [ -f "$input_dir/${pattern}.psam" ]; then
            echo "pgen:$input_dir/${pattern}"
            return 0
        fi
        
        if [ -f "$input_dir/${pattern}.bgen" ]; then
            echo "bgen:$input_dir/${pattern}.bgen"
            return 0
        fi
        
        if [ -f "$input_dir/${pattern}.vcf.gz" ]; then
            echo "vcf:$input_dir/${pattern}.vcf.gz"
            return 0
        fi
        
        if [ -f "$input_dir/${pattern}.vcf" ]; then
            echo "vcf:$input_dir/${pattern}.vcf"
            return 0
        fi
    done
    
    echo "notfound"
    return 1
}

#==========================================
# Function to Find Input File
#==========================================

find_input_file() {
    local input_dir="$1"
    local chr="$2"
    local format="$3"
    local prefix="$4"
    
    for pattern in "${prefix}${chr}" "${prefix}${chr}_merged"; do
        case "$format" in
            bfile)
                if [ -f "$input_dir/${pattern}.bed" ]; then
                    echo "$input_dir/${pattern}"
                    return 0
                fi
                ;;
            pgen)
                if [ -f "$input_dir/${pattern}.pgen" ]; then
                    echo "$input_dir/${pattern}"
                    return 0
                fi
                ;;
            vcf)
                if [ -f "$input_dir/${pattern}.vcf.gz" ]; then
                    echo "$input_dir/${pattern}.vcf.gz"
                    return 0
                elif [ -f "$input_dir/${pattern}.vcf" ]; then
                    echo "$input_dir/${pattern}.vcf"
                    return 0
                fi
                ;;
            bgen)
                if [ -f "$input_dir/${pattern}.bgen" ]; then
                    echo "$input_dir/${pattern}.bgen"
                    return 0
                fi
                ;;
        esac
    done
    
    echo ""
    return 1
}

#==========================================
# Set Output Options Based on Format
#==========================================

case "$OUTPUT_FORMAT" in
    bfile)
        OUTPUT_FLAG="--make-bed"
        OUTPUT_EXT="bed"
        OUTPUT_DESC="BFILE (bed/bim/fam)"
        ;;
    pgen)
        OUTPUT_FLAG="--make-pgen"
        OUTPUT_EXT="pgen"
        OUTPUT_DESC="PGEN (pgen/pvar/psam)"
        ;;
    vcf)
        OUTPUT_FLAG="--export vcf bgz"
        OUTPUT_EXT="vcf.gz"
        OUTPUT_DESC="VCF (vcf.gz)"
        ;;
    bgen)
        OUTPUT_FLAG="--export bgen-1.2 'bits=8'"
        OUTPUT_EXT="bgen"
        OUTPUT_DESC="BGEN v1.2 (bgen/sample)"
        ;;
esac

#==========================================
# Log Header
#==========================================

LOG_FILE="$OUTPUT_DIR/extraction.log"
SUMMARY="$OUTPUT_DIR/extraction_summary.txt"

exec > >(tee -a "$LOG_FILE") 2>&1

echo "=========================================="
echo "Genotype Extraction Pipeline"
echo "=========================================="
echo "Started: $(date)"
echo ""
echo "Configuration:"
echo "  Input directory: $INPUT_DIR"
echo "  Input prefix: '${INPUT_PREFIX}' (looking for ${INPUT_PREFIX}1, ${INPUT_PREFIX}2, ...)"
echo "  Regions directory: $REGIONS_DIR"
echo "  Output directory: $OUTPUT_DIR"
echo "  Input format: $INPUT_FORMAT"
echo "  Output format: $OUTPUT_FORMAT ($OUTPUT_DESC)"
echo "  Threads: $THREADS"
echo "  Chunk size: $CHUNK_SIZE $([ "$CHUNK_SIZE" -eq 0 ] && echo "(no chunking - extract all genes together)" || echo "genes per chunk")"
echo "  PLINK2a version: $PLINK_VERSION"
echo ""

# Summary header
{
    echo "Chromosome Extraction Summary"
    echo "Generated: $(date)"
    echo "Input prefix: ${INPUT_PREFIX}"
    echo "Input format: $INPUT_FORMAT"
    echo "Output format: $OUTPUT_FORMAT"
    echo "Chunk size: $CHUNK_SIZE $([ "$CHUNK_SIZE" -eq 0 ] && echo "(no chunking)" || echo "genes per chunk")"
    echo ""
    if [ "$CHUNK_SIZE" -eq 0 ]; then
        printf "%-5s %-25s %-10s %-8s %-15s %-15s %-10s %-12s\n" "Chr" "Input_File" "Format" "Genes" "Vars_Input" "Vars_Output" "Samples" "Status"
    else
        printf "%-5s %-8s %-25s %-10s %-8s %-8s %-15s %-15s %-10s %-12s\n" "Chr" "Chunk" "Input_File" "Format" "Genes" "Regions" "Vars_Input" "Vars_Output" "Samples" "Status"
    fi
    echo "------------------------------------------------------------------------------------------------------------------------"
} > "$SUMMARY"

TOTAL_EXTRACTED=0
TOTAL_FAILED=0
TOTAL_GENES=0
TOTAL_VARIANTS=0
TOTAL_CHUNKS=0

#==========================================
# Process Each Chromosome
#==========================================

echo "Processing chromosomes..."
echo ""

for chr in {1..22} X Y; do
    REGION_FILE="$REGIONS_DIR/chr${chr}_regions.txt"
    
    # Skip if no region file
    if [ ! -f "$REGION_FILE" ]; then
        echo "  Skipping chr${chr}: no region file"
        continue
    fi
    
    GENE_COUNT=$(wc -l < "$REGION_FILE")
    TOTAL_GENES=$((TOTAL_GENES + GENE_COUNT))
    
    echo "  Processing chr${chr}..."
    echo "    Total genes: $GENE_COUNT"
    echo "    Looking for: ${INPUT_PREFIX}${chr}.*"
    
    # Determine input file and format
    if [ "$INPUT_FORMAT" = "auto" ]; then
        DETECTED=$(detect_input_format "$INPUT_DIR" "$chr" "$INPUT_PREFIX")
        if [ "$DETECTED" = "notfound" ]; then
            echo "    ✗ Input file not found (tried all formats with prefix '${INPUT_PREFIX}')"
            if [ "$CHUNK_SIZE" -eq 0 ]; then
                printf "%-5s %-25s %-10s %-8s %-15s %-15s %-10s %-12s\n" "$chr" "${INPUT_PREFIX}${chr}" "UNKNOWN" "$GENE_COUNT" "NA" "NA" "NA" "NOT_FOUND" >> "$SUMMARY"
            else
                printf "%-5s %-8s %-25s %-10s %-8s %-8s %-15s %-15s %-10s %-12s\n" "$chr" "N/A" "${INPUT_PREFIX}${chr}" "UNKNOWN" "$GENE_COUNT" "NA" "NA" "NA" "NA" "NOT_FOUND" >> "$SUMMARY"
            fi
            TOTAL_FAILED=$((TOTAL_FAILED + 1))
            continue
        fi
        DETECTED_FORMAT="${DETECTED%%:*}"
        INPUT_FILE="${DETECTED#*:}"
        echo "    Detected format: $DETECTED_FORMAT"
    else
        INPUT_FILE=$(find_input_file "$INPUT_DIR" "$chr" "$INPUT_FORMAT" "$INPUT_PREFIX")
        if [ -z "$INPUT_FILE" ]; then
            echo "    ✗ Input file not found for format: $INPUT_FORMAT with prefix '${INPUT_PREFIX}'"
            if [ "$CHUNK_SIZE" -eq 0 ]; then
                printf "%-5s %-25s %-10s %-8s %-15s %-15s %-10s %-12s\n" "$chr" "${INPUT_PREFIX}${chr}" "$INPUT_FORMAT" "$GENE_COUNT" "NA" "NA" "NA" "NOT_FOUND" >> "$SUMMARY"
            else
                printf "%-5s %-8s %-25s %-10s %-8s %-8s %-15s %-15s %-10s %-12s\n" "$chr" "N/A" "${INPUT_PREFIX}${chr}" "$INPUT_FORMAT" "$GENE_COUNT" "NA" "NA" "NA" "NA" "NOT_FOUND" >> "$SUMMARY"
            fi
            TOTAL_FAILED=$((TOTAL_FAILED + 1))
            continue
        fi
        DETECTED_FORMAT="$INPUT_FORMAT"
    fi
    
    INPUT_BASENAME=$(basename "$INPUT_FILE")
    echo "    Input: $INPUT_BASENAME"
    
    # Get input variant count
    VARS_INPUT="NA"
    case "$DETECTED_FORMAT" in
        bfile)
            if [ -f "${INPUT_FILE}.bim" ]; then
                VARS_INPUT=$(wc -l < "${INPUT_FILE}.bim")
            fi
            ;;
        pgen)
            if [ -f "${INPUT_FILE}.pvar" ]; then
                VARS_INPUT=$(grep -v "^#" "${INPUT_FILE}.pvar" | wc -l)
            fi
            ;;
        vcf)
            if [[ "$INPUT_FILE" == *.gz ]]; then
                VARS_INPUT=$(zgrep -v "^#" "$INPUT_FILE" | wc -l)
            else
                VARS_INPUT=$(grep -v "^#" "$INPUT_FILE" | wc -l)
            fi
            ;;
        bgen)
            if [ -f "${INPUT_FILE}.bgi" ]; then
                VARS_INPUT=$(bgenix -g "$INPUT_FILE" -list 2>/dev/null | wc -l || echo "NA")
            fi
            ;;
    esac
    
    #==========================================
    # Determine Chunking Strategy
    #==========================================
    
    if [ "$CHUNK_SIZE" -eq 0 ] || [ "$GENE_COUNT" -le "$CHUNK_SIZE" ]; then
        # No chunking or genes fit in one chunk
        NUM_CHUNKS=1
        echo "    Processing all $GENE_COUNT genes in one file"
    else
        # Calculate number of chunks needed
        NUM_CHUNKS=$(( (GENE_COUNT + CHUNK_SIZE - 1) / CHUNK_SIZE ))
        echo "    Splitting into $NUM_CHUNKS chunks ($CHUNK_SIZE genes per chunk)"
    fi
    
    TOTAL_CHUNKS=$((TOTAL_CHUNKS + NUM_CHUNKS))
    
    #==========================================
    # Process Each Chunk
    #==========================================
    
    for ((chunk=1; chunk<=NUM_CHUNKS; chunk++)); do
        
        # Calculate start and end lines for this chunk
        START_LINE=$(( (chunk - 1) * CHUNK_SIZE + 1 ))
        
        if [ "$CHUNK_SIZE" -eq 0 ]; then
            # No chunking - process all lines
            END_LINE=$GENE_COUNT
            CHUNK_REGION_FILE="$REGION_FILE"
            if [ "$NUM_CHUNKS" -eq 1 ]; then
                OUTPUT_PREFIX="$OUTPUT_DIR/chr${chr}_genes"
                CHUNK_LABEL=""
            else
                OUTPUT_PREFIX="$OUTPUT_DIR/chr${chr}_genes_chunk${chunk}"
                CHUNK_LABEL="_chunk${chunk}"
            fi
        else
            # Chunking enabled
            END_LINE=$(( chunk * CHUNK_SIZE ))
            if [ "$END_LINE" -gt "$GENE_COUNT" ]; then
                END_LINE=$GENE_COUNT
            fi
            
            # Create temporary region file for this chunk
            CHUNK_REGION_FILE="$OUTPUT_DIR/.tmp_chr${chr}_chunk${chunk}_regions.txt"
            sed -n "${START_LINE},${END_LINE}p" "$REGION_FILE" > "$CHUNK_REGION_FILE"
            
            OUTPUT_PREFIX="$OUTPUT_DIR/chr${chr}_genes_chunk${chunk}"
            CHUNK_LABEL="_chunk${chunk}"
        fi
        
        CHUNK_GENE_COUNT=$(wc -l < "$CHUNK_REGION_FILE")
        
        if [ "$NUM_CHUNKS" -gt 1 ]; then
            echo ""
            echo "    === Chunk $chunk/$NUM_CHUNKS ==="
            echo "    Genes in chunk: $CHUNK_GENE_COUNT (lines $START_LINE-$END_LINE)"
        fi
        
        # Find sample file for BGEN if needed
        SAMPLE_FILE=""
        if [ "$DETECTED_FORMAT" = "bgen" ]; then
            BGEN_DIR=$(dirname "$INPUT_FILE")
            BGEN_BASE=$(basename "$INPUT_FILE" .bgen)
            
            if [ -f "${INPUT_FILE}.sample" ]; then
                SAMPLE_FILE="${INPUT_FILE}.sample"
            elif [ -f "${BGEN_DIR}/${BGEN_BASE}.sample" ]; then
                SAMPLE_FILE="${BGEN_DIR}/${BGEN_BASE}.sample"
            elif [ -f "${BGEN_DIR}/${INPUT_PREFIX}${chr}.sample" ]; then
                SAMPLE_FILE="${BGEN_DIR}/${INPUT_PREFIX}${chr}.sample"
            elif [ -f "${BGEN_DIR}/chr${chr}.sample" ]; then
                SAMPLE_FILE="${BGEN_DIR}/chr${chr}.sample"
            fi
        fi
        
        # Build PLINK2a command
        echo "    Running PLINK2a extraction..."
        
        PLINK_CMD="plink2a"
        
        case "$DETECTED_FORMAT" in
            bfile)
                PLINK_CMD="$PLINK_CMD --bfile ${INPUT_FILE}"
                ;;
            pgen)
                PLINK_CMD="$PLINK_CMD --pfile ${INPUT_FILE}"
                ;;
            vcf)
                PLINK_CMD="$PLINK_CMD --vcf ${INPUT_FILE}"
                if [[ "$INPUT_FILE" == *"dosage"* ]] || [[ "$INPUT_FILE" == *"imputed"* ]]; then
                    PLINK_CMD="$PLINK_CMD --vcf-dosage DS"
                fi
                ;;
            bgen)
                PLINK_CMD="$PLINK_CMD --bgen ${INPUT_FILE}"
                if [ -n "$SAMPLE_FILE" ]; then
                    PLINK_CMD="$PLINK_CMD --sample ${SAMPLE_FILE}"
                fi
                ;;
        esac
        
        # Add extraction and output options
        PLINK_CMD="$PLINK_CMD --extract bed1 $CHUNK_REGION_FILE"
        PLINK_CMD="$PLINK_CMD $OUTPUT_FLAG"
        PLINK_CMD="$PLINK_CMD --threads $THREADS"
        PLINK_CMD="$PLINK_CMD --out $OUTPUT_PREFIX"
        PLINK_CMD="$PLINK_CMD --memory 8000"
        
        # Execute PLINK2a
        eval $PLINK_CMD 2>&1 | grep -E "(variants|samples|written|Loading|Writing)" || true
        
        # Check if extraction succeeded
        if [ -f "${OUTPUT_PREFIX}.${OUTPUT_EXT}" ]; then
            # Count variants and samples based on output format
            VARS_OUTPUT="NA"
            SAMPLES="NA"
            
            case "$OUTPUT_FORMAT" in
                bfile)
                    if [ -f "${OUTPUT_PREFIX}.bim" ]; then
                        VARS_OUTPUT=$(wc -l < "${OUTPUT_PREFIX}.bim")
                    fi
                    if [ -f "${OUTPUT_PREFIX}.fam" ]; then
                        SAMPLES=$(wc -l < "${OUTPUT_PREFIX}.fam")
                    fi
                    ;;
                pgen)
                    if [ -f "${OUTPUT_PREFIX}.pvar" ]; then
                        VARS_OUTPUT=$(grep -v "^#" "${OUTPUT_PREFIX}.pvar" | wc -l)
                    fi
                    if [ -f "${OUTPUT_PREFIX}.psam" ]; then
                        SAMPLES=$(tail -n +2 "${OUTPUT_PREFIX}.psam" | wc -l)
                    fi
                    ;;
                vcf)
                    if [ -f "${OUTPUT_PREFIX}.vcf.gz" ]; then
                        VARS_OUTPUT=$(zgrep -v "^#" "${OUTPUT_PREFIX}.vcf.gz" | wc -l)
                        SAMPLES=$(zgrep "^#CHROM" "${OUTPUT_PREFIX}.vcf.gz" | awk '{print NF-9}')
                        
                        if command -v tabix &> /dev/null; then
                            tabix -p vcf "${OUTPUT_PREFIX}.vcf.gz" 2>/dev/null
                        fi
                    fi
                    ;;
                bgen)
                    if [ -f "${OUTPUT_PREFIX}.sample" ]; then
                        SAMPLES=$(tail -n +3 "${OUTPUT_PREFIX}.sample" | wc -l)
                    fi
                    
                    if [ -f "${OUTPUT_PREFIX}.log" ]; then
                        VARS_OUTPUT=$(grep -i "variants written" "${OUTPUT_PREFIX}.log" | grep -oP '\d+' | head -1 || echo "NA")
                    fi
                    
                    if command -v bgenix &> /dev/null; then
                        bgenix -g "${OUTPUT_PREFIX}.bgen" -index 2>/dev/null
                    fi
                    ;;
            esac
            
            # Update totals
            if [[ "$VARS_OUTPUT" =~ ^[0-9]+$ ]]; then
                TOTAL_VARIANTS=$((TOTAL_VARIANTS + VARS_OUTPUT))
            fi
            
            echo "    ✓ Extracted: $VARS_OUTPUT variants, $SAMPLES samples"
            
            if [ "$CHUNK_SIZE" -eq 0 ]; then
                printf "%-5s %-25s %-10s %-8s %-15s %-15s %-10s %-12s\n" "$chr" "$INPUT_BASENAME" "$DETECTED_FORMAT" "$CHUNK_GENE_COUNT" "$VARS_INPUT" "$VARS_OUTPUT" "$SAMPLES" "SUCCESS" >> "$SUMMARY"
            else
                printf "%-5s %-8s %-25s %-10s %-8s %-8s %-15s %-15s %-10s %-12s\n" "$chr" "$chunk" "$INPUT_BASENAME" "$DETECTED_FORMAT" "$CHUNK_GENE_COUNT" "$CHUNK_GENE_COUNT" "$VARS_INPUT" "$VARS_OUTPUT" "$SAMPLES" "SUCCESS" >> "$SUMMARY"
            fi
            
            TOTAL_EXTRACTED=$((TOTAL_EXTRACTED + 1))
        else
            echo "    ✗ Extraction failed"
            
            if [ "$CHUNK_SIZE" -eq 0 ]; then
                printf "%-5s %-25s %-10s %-8s %-15s %-15s %-10s %-12s\n" "$chr" "$INPUT_BASENAME" "$DETECTED_FORMAT" "$CHUNK_GENE_COUNT" "$VARS_INPUT" "0" "0" "FAILED" >> "$SUMMARY"
            else
                printf "%-5s %-8s %-25s %-10s %-8s %-8s %-15s %-15s %-10s %-12s\n" "$chr" "$chunk" "$INPUT_BASENAME" "$DETECTED_FORMAT" "$CHUNK_GENE_COUNT" "$CHUNK_GENE_COUNT" "$VARS_INPUT" "0" "0" "FAILED" >> "$SUMMARY"
            fi
            
            TOTAL_FAILED=$((TOTAL_FAILED + 1))
        fi
        
        # Cleanup temporary chunk region file
        if [ "$CHUNK_SIZE" -gt 0 ] && [ -f "$CHUNK_REGION_FILE" ]; then
            rm -f "$CHUNK_REGION_FILE"
        fi
        
    done  # End chunk loop
    
    echo ""
    
done  # End chromosome loop

#==========================================
# Final Summary
#==========================================

{
    echo ""
    echo "=========================================="
    echo "Overall Summary"
    echo "=========================================="
    echo "Input prefix pattern: ${INPUT_PREFIX}"
    echo "Chunk size: $CHUNK_SIZE $([ "$CHUNK_SIZE" -eq 0 ] && echo "(no chunking)" || echo "genes per chunk")"
    echo "Chromosomes processed: $((TOTAL_EXTRACTED + TOTAL_FAILED > 0 ? TOTAL_EXTRACTED + TOTAL_FAILED : 0))"
    if [ "$CHUNK_SIZE" -gt 0 ]; then
        echo "Total chunks created: $TOTAL_CHUNKS"
    fi
    echo "Successfully extracted: $TOTAL_EXTRACTED"
    echo "Failed: $TOTAL_FAILED"
    echo "Total genes: $TOTAL_GENES"
    echo "Total variants: $TOTAL_VARIANTS"
    echo ""
} | tee -a "$SUMMARY"

echo "=========================================="
echo "Genotype Extraction Complete"
echo "=========================================="
echo "Input prefix pattern: ${INPUT_PREFIX}"
echo "Chunk size: $CHUNK_SIZE $([ "$CHUNK_SIZE" -eq 0 ] && echo "(no chunking)" || echo "genes per chunk")"
echo "Chromosomes processed: $((TOTAL_EXTRACTED + TOTAL_FAILED > 0 ? TOTAL_EXTRACTED + TOTAL_FAILED : 0))"
if [ "$CHUNK_SIZE" -gt 0 ]; then
    echo "Total chunks created: $TOTAL_CHUNKS"
fi
echo "Successfully extracted: $TOTAL_EXTRACTED"
echo "Failed: $TOTAL_FAILED"
echo "Total genes: $TOTAL_GENES"
echo "Total variants: $TOTAL_VARIANTS"
echo ""

echo "Output Files ($OUTPUT_FORMAT format):"
echo "----------------------------------------"

case "$OUTPUT_FORMAT" in
    bfile)
        if [ "$CHUNK_SIZE" -eq 0 ]; then
            ls -lh "$OUTPUT_DIR"/chr*_genes.{bed,bim,fam} 2>/dev/null | \
                awk '{printf "  %-40s %10s\n", $9, $5}' || echo "  No files generated"
        else
            ls -lh "$OUTPUT_DIR"/chr*_genes_chunk*.{bed,bim,fam} 2>/dev/null | \
                awk '{printf "  %-40s %10s\n", $9, $5}' || echo "  No files generated"
        fi
        ;;
    pgen)
        if [ "$CHUNK_SIZE" -eq 0 ]; then
            ls -lh "$OUTPUT_DIR"/chr*_genes.{pgen,pvar,psam} 2>/dev/null | \
                awk '{printf "  %-40s %10s\n", $9, $5}' || echo "  No files generated"
        else
            ls -lh "$OUTPUT_DIR"/chr*_genes_chunk*.{pgen,pvar,psam} 2>/dev/null | \
                awk '{printf "  %-40s %10s\n", $9, $5}' || echo "  No files generated"
        fi
        ;;
    vcf)
        if [ "$CHUNK_SIZE" -eq 0 ]; then
            ls -lh "$OUTPUT_DIR"/chr*_genes.vcf.gz* 2>/dev/null | \
                awk '{printf "  %-40s %10s\n", $9, $5}' || echo "  No files generated"
        else
            ls -lh "$OUTPUT_DIR"/chr*_genes_chunk*.vcf.gz* 2>/dev/null | \
                awk '{printf "  %-40s %10s\n", $9, $5}' || echo "  No files generated"
        fi
        ;;
    bgen)
        if [ "$CHUNK_SIZE" -eq 0 ]; then
            ls -lh "$OUTPUT_DIR"/chr*_genes.{bgen,sample} 2>/dev/null | \
                awk '{printf "  %-40s %10s\n", $9, $5}' || echo "  No files generated"
            if ls "$OUTPUT_DIR"/chr*_genes.bgen.bgi &>/dev/null; then
                echo ""
                echo "  BGEN index files (.bgi) also created"
            fi
        else
            ls -lh "$OUTPUT_DIR"/chr*_genes_chunk*.{bgen,sample} 2>/dev/null | \
                awk '{printf "  %-40s %10s\n", $9, $5}' || echo "  No files generated"
            if ls "$OUTPUT_DIR"/chr*_genes_chunk*.bgen.bgi &>/dev/null; then
                echo ""
                echo "  BGEN index files (.bgi) also created"
            fi
        fi
        ;;
esac

echo ""
echo "Summary Table:"
echo "----------------------------------------"
column -t -s $'\t' "$SUMMARY" | sed 's/^/  /'

echo ""
echo "=========================================="
echo "Next Steps"
echo "=========================================="
echo ""

if [ "$CHUNK_SIZE" -gt 0 ]; then
    echo "Note: Output is split into chunks of $CHUNK_SIZE genes each"
    echo ""
fi

case "$OUTPUT_FORMAT" in
    bfile)
        echo "1. Verify extraction:"
        if [ "$CHUNK_SIZE" -eq 0 ]; then
            echo "   plink2a --bfile $OUTPUT_DIR/chr1_genes --freq --out $OUTPUT_DIR/chr1_freq"
        else
            echo "   plink2a --bfile $OUTPUT_DIR/chr1_genes_chunk1 --freq --out $OUTPUT_DIR/chr1_chunk1_freq"
        fi
        echo ""
        echo "2. List all output files:"
        if [ "$CHUNK_SIZE" -eq 0 ]; then
            echo "   ls -lh $OUTPUT_DIR/chr*_genes.bed"
        else
            echo "   ls -lh $OUTPUT_DIR/chr*_genes_chunk*.bed"
        fi
        echo ""
        if [ "$CHUNK_SIZE" -gt 0 ]; then
            echo "3. Process all chunks in a loop:"
            echo "   for chunk in $OUTPUT_DIR/chr1_genes_chunk*.bed; do"
            echo "     base=${chunk%.bed}"
            echo "     plink2a --bfile $base --freq --out ${base}_freq"
            echo "   done"
            echo ""
            echo "4. Merge chunks if needed (for same chromosome):"
            echo "   # Create merge list"
            echo "   ls $OUTPUT_DIR/chr1_genes_chunk*.bed | sed 's/.bed//' > merge_list.txt"
            echo "   # Merge using PLINK"
            echo "   plink2a --pmerge-list merge_list.txt --make-bed --out $OUTPUT_DIR/chr1_genes_merged"
            echo ""
        fi
        ;;
    
    pgen)
        echo "1. Verify extraction:"
        if [ "$CHUNK_SIZE" -eq 0 ]; then
            echo "   plink2a --pfile $OUTPUT_DIR/chr1_genes --freq --out $OUTPUT_DIR/chr1_freq"
        else
            echo "   plink2a --pfile $OUTPUT_DIR/chr1_genes_chunk1 --freq --out $OUTPUT_DIR/chr1_chunk1_freq"
        fi
        echo ""
        if [ "$CHUNK_SIZE" -gt 0 ]; then
            echo "2. Process all chunks:"
            echo "   for chunk in $OUTPUT_DIR/chr1_genes_chunk*.pgen; do"
            echo "     base=${chunk%.pgen}"
            echo "     plink2a --pfile $base --freq --out ${base}_freq"
            echo "   done"
            echo ""
        fi
        ;;
    
    vcf)
        echo "1. Verify extraction:"
        if [ "$CHUNK_SIZE" -eq 0 ]; then
            echo "   bcftools stats $OUTPUT_DIR/chr1_genes.vcf.gz > $OUTPUT_DIR/chr1_stats.txt"
        else
            echo "   bcftools stats $OUTPUT_DIR/chr1_genes_chunk1.vcf.gz > $OUTPUT_DIR/chr1_chunk1_stats.txt"
        fi
        echo ""
        if [ "$CHUNK_SIZE" -gt 0 ]; then
            echo "2. Process all chunks:"
            echo "   for chunk in $OUTPUT_DIR/chr1_genes_chunk*.vcf.gz; do"
            echo "     bcftools stats $chunk > ${chunk%.vcf.gz}_stats.txt"
            echo "   done"
            echo ""
            echo "3. Merge chunks if needed:"
            echo "   bcftools concat $OUTPUT_DIR/chr1_genes_chunk*.vcf.gz -Oz -o $OUTPUT_DIR/chr1_genes_merged.vcf.gz"
            echo "   tabix -p vcf $OUTPUT_DIR/chr1_genes_merged.vcf.gz"
            echo ""
        fi
        ;;
    
    bgen)
        echo "1. Verify extraction:"
        if [ "$CHUNK_SIZE" -eq 0 ]; then
            echo "   bgenix -g $OUTPUT_DIR/chr1_genes.bgen -list | head"
        else
            echo "   bgenix -g $OUTPUT_DIR/chr1_genes_chunk1.bgen -list | head"
        fi
        echo ""
        echo "2. Use in SAIGE-GENE:"
        if [ "$CHUNK_SIZE" -eq 0 ]; then
            echo "   Rscript step2_SPAtests.R \\"
            echo "     --bgenFile=$OUTPUT_DIR/chr1_genes.bgen \\"
            echo "     --bgenFileIndex=$OUTPUT_DIR/chr1_genes.bgen.bgi \\"
            echo "     --sampleFile=$OUTPUT_DIR/chr1_genes.sample \\"
            echo "     --groupFile=chr1_groups.txt \\"
            echo "     --GMMATmodelFile=step1_null_model.rda \\"
            echo "     --varianceRatioFile=step1_variance_ratio.txt \\"
            echo "     --SAIGEOutputFile=chr1_gene_results.txt \\"
            echo "     --is_Firth_beta=TRUE \\"
            echo "     --pCutoffforFirth=0.01 \\"
            echo "     --LOCO=FALSE"
        else
            echo "   # Process each chunk separately"
            echo "   for i in {1..$TOTAL_CHUNKS}; do"
            echo "     chunk_file=$OUTPUT_DIR/chr1_genes_chunk${i}.bgen"
            echo "     if [ -f $chunk_file ]; then"
            echo "       Rscript step2_SPAtests.R \\"
            echo "         --bgenFile=$chunk_file \\"
            echo "         --bgenFileIndex=${chunk_file}.bgi \\"
            echo "         --sampleFile=${chunk_file%.bgen}.sample \\"
            echo "         --groupFile=chr1_groups_chunk${i}.txt \\"
            echo "         --GMMATmodelFile=step1_null_model.rda \\"
            echo "         --varianceRatioFile=step1_variance_ratio.txt \\"
            echo "         --SAIGEOutputFile=chr1_chunk${i}_results.txt \\"
            echo "         --is_Firth_beta=TRUE \\"
            echo "         --pCutoffforFirth=0.01 \\"
            echo "         --LOCO=FALSE"
            echo "     fi"
            echo "   done"
            echo ""
            echo "   # Combine results from all chunks"
            echo "   head -1 chr1_chunk1_results.txt > chr1_all_results.txt"
            echo "   tail -n +2 -q chr1_chunk*_results.txt >> chr1_all_results.txt"
        fi
        echo ""
        if [ "$CHUNK_SIZE" -gt 0 ]; then
            echo "3. Note: You may need to split your group file to match chunks:"
            echo "   # This ensures each chunk processes only its genes"
            echo "   # Extract gene lists from each chunk and create corresponding group files"
            echo ""
        fi
        ;;
esac

if [ "$CHUNK_SIZE" -gt 0 ]; then
    echo ""
    echo "=========================================="
    echo "Chunking Information"
    echo "=========================================="
    echo "Total genes processed: $TOTAL_GENES"
    echo "Genes per chunk: $CHUNK_SIZE"
    echo "Total chunks created: $TOTAL_CHUNKS"
    echo ""
    echo "Chunk file naming pattern:"
    echo "  chr{N}_genes_chunk{M}.{ext}"
    echo "  Where N = chromosome number, M = chunk number"
    echo ""
    echo "Benefits of chunking:"
    echo "  ✓ Reduced memory usage per job"
    echo "  ✓ Parallel processing capability"
    echo "  ✓ Easier error recovery (re-run failed chunks)"
    echo "  ✓ Better for HPC job schedulers"
    echo ""
fi

echo "=========================================="
echo "File Summary"
echo "=========================================="
if [ "$CHUNK_SIZE" -eq 0 ]; then
    echo "Output pattern: chr{1-22,X,Y}_genes.{ext}"
else
    echo "Output pattern: chr{1-22,X,Y}_genes_chunk{1-N}.{ext}"
fi
echo "Input prefix used: ${INPUT_PREFIX}"
echo ""
echo "Log Files:"
echo "  Main log: $LOG_FILE"
echo "  Summary:  $SUMMARY"
echo ""
echo "Completed: $(date)"
echo "=========================================="

exit 0

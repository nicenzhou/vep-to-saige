#!/bin/bash

# Script: step8_run_saige_gene_tests.sh
# Description: Run SAIGE-GENE Step 2 tests on extracted genotype chunks
# Usage: step8_run_saige_gene_tests.sh <config_file>

set -euo pipefail

#==========================================
# Configuration File Format
#==========================================

if [ $# -eq 0 ]; then
    cat << 'EOF'
Usage: step8_run_saige_gene_tests.sh <config_file>

Description:
  Run SAIGE-GENE Step 2 association tests on extracted genotype files.
  Supports both chunked and non-chunked files, multiple input formats.
  Handles chromosome naming: chr1/chr01, 1/01 formats.
  Optional automatic merging of chunk results by chromosome.
  Uses a configuration file for all parameters.

Configuration File Format (key=value pairs):
  # Required parameters
  GENOTYPE_DIR=/path/to/gene_bfiles
  OUTPUT_DIR=/path/to/saige_results
  GMMAT_MODEL=/path/to/step1_null_model.rda
  VARIANCE_RATIO=/path/to/step1_variance_ratio.txt
  
  # Group file options (choose one)
  GROUP_FILE=/path/to/all_genes_groups.txt     # Single file for all chromosomes
  GROUP_FILE_BY_CHR=yes                         # Use chr{N}_groups.txt files
  GROUP_DIR=/path/to/group_files                # Directory with group files (if BY_CHR=yes)
  
  # Input format
  INPUT_FORMAT=bgen          # bgen, bfile, pgen, vcf
  
  # Optional: File naming format
  CHR_PREFIX=chr             # Prefix for chromosome in filenames (default: chr)
  CHR_PADDING=no             # Whether chromosomes are zero-padded (01, 02) or not (1, 2). Options: yes/no/auto (default: auto)
  CHUNKED_INPUT=yes          # Are genotype files chunked? (default: no)
  CHUNK_PATTERN=chunk        # Pattern to identify chunks (default: chunk)
  
  # Optional: Processing control
  THREADS=4
  CHROMOSOMES=1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22  # Comma-separated
  
  # Optional: Result merging
  MERGE_CHUNKS=yes           # Merge chunk results by chromosome (default: no)
  KEEP_CHUNK_FILES=no        # Keep individual chunk files after merging (default: yes)
  
  # Optional: Allele configuration
  ALLELE_ORDER=alt-first     # alt-first or ref-first (default: alt-first)
  
  # Optional: SAIGE parameters (any SAIGE step2 parameter)
  LOCO=TRUE
  chrom=                     # Will be set automatically per chromosome
  minMAF=0.0001
  minMAC=1
  maxMAF_in_groupTest=0.0001,0.001,0.01
  annotation_in_groupTest=lof,missense;lof,missense;lof;synonymous
  is_Firth_beta=TRUE
  pCutoffforFirth=0.01
  is_output_markerList_in_groupTest=TRUE
  r_corr=0                   # NOTE: Use r_corr in config (undercore), will be converted to r.corr for SAIGE
  SPAcutoff=2
  
  # Add any other SAIGE step2_SPAtests.R parameters as needed

Example config file (saige_config.txt):
  GENOTYPE_DIR=gene_bfiles
  OUTPUT_DIR=saige_results
  GMMAT_MODEL=step1_null_model.rda
  VARIANCE_RATIO=step1_variance_ratio.txt
  GROUP_FILE_BY_CHR=yes
  GROUP_DIR=group_files
  INPUT_FORMAT=bgen
  CHR_PREFIX=chr
  CHR_PADDING=auto
  CHUNKED_INPUT=yes
  CHUNK_PATTERN=chunk
  ALLELE_ORDER=alt-first
  THREADS=8
  CHROMOSOMES=1,2,3,21,22
  MERGE_CHUNKS=yes
  KEEP_CHUNK_FILES=no
  LOCO=TRUE
  minMAF=0.0001
  is_Firth_beta=TRUE
  pCutoffforFirth=0.01
  r_corr=0

Usage:
  ./step8_run_saige_gene_tests.sh saige_config.txt

Chromosome Naming:
  The script automatically detects chromosome naming formats:
  - chr1, chr2, ... chr22, chrX, chrY
  - chr01, chr02, ... chr22, chrX, chrY
  - 1, 2, ... 22, X, Y
  - 01, 02, ... 22, X, Y

Result Merging:
  MERGE_CHUNKS=yes : Combine chunk results into chr{N}_combined_results.txt
  MERGE_CHUNKS=no  : Keep individual chunk result files
  
  KEEP_CHUNK_FILES=yes : Keep original chunk files after merging
  KEEP_CHUNK_FILES=no  : Delete chunk files after successful merge

EOF
    exit 1
fi

CONFIG_FILE="\$1"

if [ ! -f "$CONFIG_FILE" ]; then
    echo "ERROR: Configuration file not found: $CONFIG_FILE" >&2
    exit 1
fi

#==========================================
# Helper Functions
#==========================================

# Function to format chromosome number based on detected padding
format_chr() {
    local chr="\$1"
    local padding="\$2"
    
    # Handle X and Y chromosomes
    if [[ "$$chr" =~ ^[XYxy]$$ ]]; then
        echo "$$(echo $$chr | tr '[:lower:]' '[:upper:]')"
        return
    fi
    
    # Remove leading zeros for comparison
    local chr_num=$$(echo "$$chr" | sed 's/^0*//')
    
    if [ "$padding" = "yes" ]; then
        # Zero-pad to 2 digits
        printf "%02d" "$chr_num"
    else
        # No padding
        echo "$chr_num"
    fi
}

# Function to detect chromosome padding in files
detect_chr_padding() {
    local dir="\$1"
    local ext="\$2"
    local prefix="\$3"
    
    # Check for padded format (chr01, chr02, etc.)
    if ls "$${dir}/$${prefix}01"* 2>/dev/null | grep -q "\.$${ext}$$"; then
        echo "yes"
        return
    fi
    
    # Check for padded format without prefix (01, 02, etc.)
    if ls "$${dir}/01"* 2>/dev/null | grep -q "\.$${ext}$"; then
        echo "yes"
        return
    fi
    
    # Check for non-padded format (chr1, chr2, etc.)
    if ls "$${dir}/$${prefix}1"* 2>/dev/null | grep -q "\.$${ext}$$"; then
        echo "no"
        return
    fi
    
    # Check for non-padded format without prefix (1, 2, etc.)
    if ls "$${dir}/1"* 2>/dev/null | grep -q "\.$${ext}$"; then
        echo "no"
        return
    fi
    
    # Default to no padding
    echo "no"
}

# Function to find genotype file with flexible naming
find_genotype_file() {
    local dir="\$1"
    local chr="\$2"
    local ext="\$3"
    local prefix="\$4"
    local padding="\$5"
    local pattern="\$6"
    
    local chr_formatted=$$(format_chr "$$chr" "$padding")
    local found_files=()
    
    # Try different naming patterns
    local patterns=(
        "$${dir}/$${prefix}$${chr_formatted}_genes_$${pattern}*.${ext}"
        "$${dir}/$${prefix}$${chr_formatted}_genes.$${ext}"
        "$${dir}/$${chr_formatted}_genes_$${pattern}*.$${ext}"
        "$${dir}/$${chr_formatted}_genes.${ext}"
    )
    
    for pat in "${patterns[@]}"; do
        for file in $pat; do
            if [ -f "$file" ]; then
                found_files+=("$file")
            fi
        done
    done
    
    # Return unique files
    printf '%s\n' "${found_files[@]}" | sort -u
}

# Function to find group file with flexible naming
find_group_file() {
    local dir="\$1"
    local chr="\$2"
    local prefix="\$3"
    local padding="\$4"
    
    local chr_formatted=$$(format_chr "$$chr" "$padding")
    
    # Try different naming patterns
    local patterns=(
        "$${dir}/$${prefix}${chr_formatted}_groups.txt"
        "$${dir}/$${prefix}${chr_formatted}_group.txt"
        "$${dir}/$${chr_formatted}_groups.txt"
        "$${dir}/$${chr_formatted}_group.txt"
        "$${dir}/chr$${chr_formatted}_groups.txt"
        "$${dir}/chr$${chr_formatted}_group.txt"
    )
    
    for pat in "${patterns[@]}"; do
        if [ -f "$pat" ]; then
            echo "$pat"
            return 0
        fi
    done
    
    return 1
}

# Function to merge chunk results by chromosome
merge_chunk_results() {
    local chr="\$1"
    local output_dir="\$2"
    local keep_chunks="\$3"
    
    local chunk_files=("$$output_dir"/chr$${chr}_chunk*_results.txt)
    
    # Check if any chunk files exist
    if [ ! -f "${chunk_files[0]}" ]; then
        echo "  No chunk files to merge for chr${chr}"
        return 1
    fi
    
    local combined_file="$$output_dir/chr$${chr}_combined_results.txt"
    
    echo "  Merging $${#chunk_files[@]} chunk files for chr$${chr}..."
    
    # Sort chunk files by chunk number to maintain order
    local sorted_chunks=($$(printf '%s\n' "$${chunk_files[@]}" | sort -V))
    
    # Get header from first chunk
    if [ -f "${sorted_chunks[0]}" ]; then
        head -1 "$${sorted_chunks[0]}" > "$$combined_file"
        
        # Append data from all chunks (skip headers)
        for chunk_file in "${sorted_chunks[@]}"; do
            if [ -f "$chunk_file" ]; then
                tail -n +2 "$$chunk_file" >> "$$combined_file"
            fi
        done
        
        # Count total results
        local total_lines=$$(($$(wc -l < "$combined_file") - 1))
        echo "    ✓ Combined file created: $$(basename $$combined_file)"
        echo "    Total results: $total_lines"
        
        # Remove chunk files if requested
        if [ "$keep_chunks" = "no" ]; then
            echo "    Removing individual chunk files..."
            for chunk_file in "${sorted_chunks[@]}"; do
                rm -f "$chunk_file"
            done
            echo "    ✓ Chunk files removed"
        fi
        
        return 0
    else
        echo "    ✗ Failed to merge: First chunk file not found"
        return 1
    fi
}

#==========================================
# Load Configuration
#==========================================

echo "=========================================="
echo "SAIGE-GENE Step 2 Test Runner"
echo "=========================================="
echo "Started: $(date)"
echo ""
echo "Loading configuration from: $CONFIG_FILE"
echo ""

# Source configuration file
source "$CONFIG_FILE"

# Validate required parameters
REQUIRED_PARAMS=("GENOTYPE_DIR" "OUTPUT_DIR" "GMMAT_MODEL" "VARIANCE_RATIO")
for param in "${REQUIRED_PARAMS[@]}"; do
    if [ -z "${!param:-}" ]; then
        echo "ERROR: Required parameter $param not set in config file" >&2
        exit 1
    fi
done

# Set defaults
THREADS="${THREADS:-4}"
INPUT_FORMAT="${INPUT_FORMAT:-bgen}"
CHUNK_PATTERN="${CHUNK_PATTERN:-chunk}"
CHUNKED_INPUT="${CHUNKED_INPUT:-no}"
CHROMOSOMES="${CHROMOSOMES:-1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}"
GROUP_FILE_BY_CHR="${GROUP_FILE_BY_CHR:-no}"
CHR_PREFIX="${CHR_PREFIX:-chr}"
CHR_PADDING="${CHR_PADDING:-auto}"
MERGE_CHUNKS="${MERGE_CHUNKS:-no}"
KEEP_CHUNK_FILES="${KEEP_CHUNK_FILES:-yes}"
ALLELE_ORDER="${ALLELE_ORDER:-alt-first}"

# Validate allele order
if [[ ! "$$ALLELE_ORDER" =~ ^(alt-first|ref-first)$$ ]]; then
    echo "ERROR: ALLELE_ORDER must be 'alt-first' or 'ref-first', got: $ALLELE_ORDER" >&2
    exit 1
fi

# Validate directories
[ ! -d "$$GENOTYPE_DIR" ] && { echo "ERROR: GENOTYPE_DIR not found: $$GENOTYPE_DIR" >&2; exit 1; }
[ ! -f "$$GMMAT_MODEL" ] && { echo "ERROR: GMMAT_MODEL not found: $$GMMAT_MODEL" >&2; exit 1; }
[ ! -f "$$VARIANCE_RATIO" ] && { echo "ERROR: VARIANCE_RATIO not found: $$VARIANCE_RATIO" >&2; exit 1; }

mkdir -p "$OUTPUT_DIR"

# Validate group file configuration
if [ "$GROUP_FILE_BY_CHR" = "yes" ]; then
    GROUP_DIR="$${GROUP_DIR:-$$GENOTYPE_DIR}"
    if [ ! -d "$GROUP_DIR" ]; then
        echo "ERROR: GROUP_DIR not found: $GROUP_DIR" >&2
        exit 1
    fi
    echo "Using chromosome-specific group files from: $GROUP_DIR"
else
    if [ -z "${GROUP_FILE:-}" ]; then
        echo "ERROR: Either GROUP_FILE or GROUP_FILE_BY_CHR=yes must be specified" >&2
        exit 1
    fi
    if [ ! -f "$GROUP_FILE" ]; then
        echo "ERROR: GROUP_FILE not found: $GROUP_FILE" >&2
        exit 1
    fi
    echo "Using single group file: $GROUP_FILE"
fi

#==========================================
# Determine File Extensions
#==========================================

INPUT_FORMAT=$$(echo "$$INPUT_FORMAT" | tr '[:upper:]' '[:lower:]')

case "$INPUT_FORMAT" in
    bgen)
        GENO_EXT="bgen"
        SAMPLE_EXT="sample"
        INDEX_EXT="bgen.bgi"
        SAIGE_FLAG="--bgenFile"
        SAMPLE_FLAG="--sampleFile"
        INDEX_FLAG="--bgenFileIndex"
        ;;
    bfile|bed)
        GENO_EXT="bed"
        SAIGE_FLAG="--bedFile"
        BIM_FLAG="--bimFile"
        FAM_FLAG="--famFile"
        ;;
    pgen)
        GENO_EXT="pgen"
        PVAR_EXT="pvar"
        PSAM_EXT="psam"
        SAIGE_FLAG="--pgenFile"
        PVAR_FLAG="--pvarFile"
        PSAM_FLAG="--psamFile"
        ;;
    vcf)
        GENO_EXT="vcf.gz"
        INDEX_EXT="vcf.gz.tbi"
        SAIGE_FLAG="--vcfFile"
        INDEX_FLAG="--vcfFileIndex"
        ;;
    *)
        echo "ERROR: Unsupported INPUT_FORMAT: $INPUT_FORMAT" >&2
        echo "Supported formats: bgen, bfile, pgen, vcf" >&2
        exit 1
        ;;
esac

#==========================================
# Auto-detect Chromosome Padding
#==========================================

if [ "$CHR_PADDING" = "auto" ]; then
    echo "Auto-detecting chromosome padding format..."
    CHR_PADDING=$$(detect_chr_padding "$$GENOTYPE_DIR" "$$GENO_EXT" "$$CHR_PREFIX")
    echo "  Detected padding: $CHR_PADDING"
    
    if [ "$GROUP_FILE_BY_CHR" = "yes" ]; then
        GROUP_PADDING=$$(detect_chr_padding "$$GROUP_DIR" "txt" "$CHR_PREFIX")
        echo "  Group file padding: $GROUP_PADDING"
    fi
else
    GROUP_PADDING="$CHR_PADDING"
fi

echo ""

#==========================================
# Load SAIGE Module
#==========================================

echo "Loading SAIGE module..."
if command -v module &> /dev/null; then
    if module avail SAIGE 2>&1 | grep -q SAIGE; then
        module load SAIGE
        echo "  ✓ SAIGE module loaded"
    elif module avail saige 2>&1 | grep -q saige; then
        module load saige
        echo "  ✓ saige module loaded"
    else
        echo "  ⚠ No SAIGE module found, using system SAIGE"
    fi
else
    echo "  ⚠ Module system not available, using system SAIGE"
fi

# Check for SAIGE step2 script
if ! command -v step2_SPAtests.R &> /dev/null; then
    echo "ERROR: step2_SPAtests.R not found in PATH" >&2
    echo "Please load SAIGE module or add to PATH" >&2
    exit 1
fi

echo ""

#==========================================
# Build SAIGE Parameter String
#==========================================

build_saige_params() {
    local chr="\$1"
    local params=""
    
    # Add chromosome-specific parameter if LOCO is enabled
    if [ "${LOCO:-FALSE}" = "TRUE" ]; then
        # Use original chromosome number (no prefix, no padding) for SAIGE
        local chr_num=$$(echo "$$chr" | sed 's/^0*//')
        params="$$params --chrom=$$chr_num"
    fi
    
    # Add AlleleOrder parameter
    params="$$params --AlleleOrder=$$ALLELE_ORDER"
    
    # Add optional SAIGE parameters from config
    # NOTE: r_corr in config becomes r.corr for SAIGE
    for param in minMAF minMAC minInfo maxMissing \
                 maxMAF_in_groupTest maxMAC_in_groupTest \
                 annotation_in_groupTest is_Firth_beta pCutoffforFirth \
                 is_output_markerList_in_groupTest SPAcutoff \
                 markers_per_chunk groups_per_chunk \
                 is_output_moreDetails is_overwrite_output \
                 MACCutoff_to_CollapseUltraRare weights_beta \
                 markers_per_chunk_in_groupTest is_single_in_groupTest \
                 is_no_weight_in_groupTest is_fastTest \
                 max_MAC_for_ER dosage_zerod_cutoff dosage_zerod_MAC_cutoff \
                 minGroupMAC_in_BurdenTest impute_method is_imputed_data \
                 cateVarRatioMinMACVecExclude cateVarRatioMaxMACVecInclude; do
        
        if [ -n "${!param:-}" ]; then
            params="$$params --$${param}=${!param}"
        fi
    done
    
    # Handle r_corr specially (config uses r_corr, SAIGE expects r.corr)
    if [ -n "${r_corr:-}" ]; then
        params="$$params --r.corr=$${r_corr}"
    fi
    
    # Add LOCO flag
    if [ "${LOCO:-FALSE}" = "TRUE" ]; then
        params="$params --LOCO=TRUE"
    fi
    
    echo "$params"
}

#==========================================
# Process Each Chromosome
#==========================================

LOG_FILE="$OUTPUT_DIR/saige_run.log"
SUMMARY_FILE="$OUTPUT_DIR/saige_run_summary.txt"

exec > >(tee -a "$LOG_FILE") 2>&1

echo "=========================================="
echo "Configuration Summary"
echo "=========================================="
echo "Genotype directory: $GENOTYPE_DIR"
echo "Output directory: $OUTPUT_DIR"
echo "Input format: $INPUT_FORMAT"
echo "Chromosome prefix: $CHR_PREFIX"
echo "Chromosome padding: $CHR_PADDING"
echo "Chunked input: $CHUNKED_INPUT"
if [ "$CHUNKED_INPUT" = "yes" ]; then
    echo "Chunk pattern: $CHUNK_PATTERN"
fi
echo "Allele order: $ALLELE_ORDER"
echo "GMMAT model: $GMMAT_MODEL"
echo "Variance ratio: $VARIANCE_RATIO"
echo "Threads: $THREADS"
echo "Chromosomes: $CHROMOSOMES"
echo "Merge chunks: $MERGE_CHUNKS"
if [ "$MERGE_CHUNKS" = "yes" ]; then
    echo "Keep chunk files: $KEEP_CHUNK_FILES"
fi
echo ""
echo "Additional SAIGE parameters:"
echo "$(build_saige_params 1)" | tr ' ' '\n' | grep '^--' | sed 's/^/  /'
echo ""
echo "=========================================="
echo ""

# Convert chromosome list to array
IFS=',' read -ra CHR_ARRAY <<< "$CHROMOSOMES"

TOTAL_JOBS=0
SUCCESSFUL_JOBS=0
FAILED_JOBS=0

# Track chromosomes with chunks for merging
declare -A CHR_HAS_CHUNKS

{
    echo "SAIGE-GENE Run Summary"
    echo "Generated: $(date)"
    echo ""
    printf "%-5s %-10s %-45s %-40s %-10s\n" "Chr" "Chunk" "Genotype_File" "Output_File" "Status"
    echo "---------------------------------------------------------------------------------------------------------------"
} > "$SUMMARY_FILE"

for chr in "${CHR_ARRAY[@]}"; do
    
    echo "=========================================="
    echo "Processing Chromosome $chr"
    echo "=========================================="
    
    # Format chromosome for file names
    CHR_FORMATTED=$$(format_chr "$$chr" "$CHR_PADDING")
    
    # Determine group file
    if [ "$GROUP_FILE_BY_CHR" = "yes" ]; then
        CHR_GROUP_FILE=$$(find_group_file "$$GROUP_DIR" "$$chr" "$$CHR_PREFIX" "$GROUP_PADDING")
        if [ -z "$CHR_GROUP_FILE" ]; then
            echo "  ✗ Group file not found for chromosome $chr"
            echo "    Tried patterns with chr format: $${CHR_PREFIX}$${CHR_FORMATTED}"
            echo "  Skipping chromosome $chr"
            echo ""
            continue
        fi
    else
        CHR_GROUP_FILE="$GROUP_FILE"
    fi
    
    echo "  Group file: $CHR_GROUP_FILE"
    
    # Find all genotype files for this chromosome
    GENO_FILES=()
    
    # Use the helper function to find files
    while IFS= read -r geno_file; do
        GENO_FILES+=("$geno_file")
    done < <(find_genotype_file "$$GENOTYPE_DIR" "$$chr" "$$GENO_EXT" "$$CHR_PREFIX" "$$CHR_PADDING" "$$CHUNK_PATTERN")
    
    if [ ${#GENO_FILES[@]} -eq 0 ]; then
        echo "  ✗ No genotype files found for chromosome $chr"
        echo "    Tried patterns with chr format: $${CHR_PREFIX}$${CHR_FORMATTED}"
        echo ""
        continue
    fi
    
    echo "  Found ${#GENO_FILES[@]} genotype file(s)"
    
    # Check if this chromosome has chunks
    if [ ${#GENO_FILES[@]} -gt 1 ]; then
        CHR_HAS_CHUNKS[$chr]=1
        echo "  Chunks detected: will merge after processing" 
    fi
    
    echo ""
    
    # Build chromosome-specific SAIGE parameters
    SAIGE_CHR_PARAMS=$$(build_saige_params "$$chr")
    
    # Process each genotype file
    for geno_file in "${GENO_FILES[@]}"; do
        
        TOTAL_JOBS=$((TOTAL_JOBS + 1))
        
        GENO_BASE=$$(basename "$$geno_file" .${GENO_EXT})
        
        # Determine chunk number if applicable
        if [[ "$$GENO_BASE" =~ $${CHUNK_PATTERN}([0-9]+) ]]; then
            CHUNK_NUM="${BASH_REMATCH[1]}"
            CHUNK_LABEL="chunk${CHUNK_NUM}"
        else
            CHUNK_NUM="all"
            CHUNK_LABEL="all"
        fi
        
        OUTPUT_FILE="$$OUTPUT_DIR/chr$${chr}_${CHUNK_LABEL}_results.txt"
        
        echo "  ========================================"
        echo "  Job $$TOTAL_JOBS: Chr$$chr $CHUNK_LABEL"
        echo "  Genotype: $GENO_BASE"
        echo "  Output: $$(basename $$OUTPUT_FILE)"
        echo ""
        
        # Build SAIGE command
        SAIGE_CMD="step2_SPAtests.R"
        
        # Add genotype file(s)
        case "$INPUT_FORMAT" in
            bgen)
                GENO_PREFIX="$${geno_file%.$${GENO_EXT}}"
                SAIGE_CMD="$$SAIGE_CMD $$SAIGE_FLAG=$geno_file"
                
                # Sample file
                if [ -f "$${GENO_PREFIX}.$${SAMPLE_EXT}" ]; then
                    SAIGE_CMD="$$SAIGE_CMD $$SAMPLE_FLAG=$${GENO_PREFIX}.$${SAMPLE_EXT}"
                elif [ -f "$${GENOTYPE_DIR}/$${CHR_PREFIX}$${CHR_FORMATTED}.$${SAMPLE_EXT}" ]; then
                    SAIGE_CMD="$$SAIGE_CMD $$SAMPLE_FLAG=$${GENOTYPE_DIR}/$${CHR_PREFIX}$${CHR_FORMATTED}.$${SAMPLE_EXT}"
                fi
                
                # Index file
                if [ -f "$${geno_file}.$${INDEX_EXT}" ]; then
                    SAIGE_CMD="$$SAIGE_CMD $$INDEX_FLAG=$${geno_file}.$${INDEX_EXT}"
                fi
                ;;
            
            bfile)
                GENO_PREFIX="$${geno_file%.$${GENO_EXT}}"
                SAIGE_CMD="$$SAIGE_CMD $$SAIGE_FLAG=$geno_file"
                SAIGE_CMD="$$SAIGE_CMD $$BIM_FLAG=${GENO_PREFIX}.bim"
                SAIGE_CMD="$$SAIGE_CMD $$FAM_FLAG=${GENO_PREFIX}.fam"
                ;;
            
            pgen)
                GENO_PREFIX="$${geno_file%.$${GENO_EXT}}"
                SAIGE_CMD="$$SAIGE_CMD $$SAIGE_FLAG=$geno_file"
                SAIGE_CMD="$$SAIGE_CMD $$PVAR_FLAG=$${GENO_PREFIX}.$${PVAR_EXT}"
                SAIGE_CMD="$$SAIGE_CMD $$PSAM_FLAG=$${GENO_PREFIX}.$${PSAM_EXT}"
                ;;
            
            vcf)
                SAIGE_CMD="$$SAIGE_CMD $$SAIGE_FLAG=$geno_file"
                
                # Index file
                if [ -f "$${geno_file}.$${INDEX_EXT}" ]; then
                    SAIGE_CMD="$$SAIGE_CMD $$INDEX_FLAG=$${geno_file}.$${INDEX_EXT}"
                fi
                ;;
        esac
        
        # Add required parameters
        SAIGE_CMD="$$SAIGE_CMD --GMMATmodelFile=$$GMMAT_MODEL"
        SAIGE_CMD="$$SAIGE_CMD --varianceRatioFile=$$VARIANCE_RATIO"
        SAIGE_CMD="$$SAIGE_CMD --groupFile=$$CHR_GROUP_FILE"
        SAIGE_CMD="$$SAIGE_CMD --SAIGEOutputFile=$$OUTPUT_FILE"
        
        # Add chromosome-specific and optional parameters
        SAIGE_CMD="$$SAIGE_CMD $$SAIGE_CHR_PARAMS"
        
        # Display command
        echo "  Command:"
        echo "$SAIGE_CMD" | sed 's/ --/\n    --/g' | sed 's/^/    /'
        echo ""
        
        # Execute SAIGE
        echo "  Running SAIGE..."
        if eval "$SAIGE_CMD"; then
            echo "  ✓ Job completed successfully"
            SUCCESSFUL_JOBS=$((SUCCESSFUL_JOBS + 1))
            printf "%-5s %-10s %-45s %-40s %-10s\n" "$$chr" "$$CHUNK_LABEL" "$$GENO_BASE" "$$(basename $$OUTPUT_FILE)" "SUCCESS" >> "$$SUMMARY_FILE"
        else
            echo "  ✗ Job failed"
            FAILED_JOBS=$((FAILED_JOBS + 1))
            printf "%-5s %-10s %-45s %-40s %-10s\n" "$$chr" "$$CHUNK_LABEL" "$$GENO_BASE" "$$(basename $$OUTPUT_FILE)" "FAILED" >> "$$SUMMARY_FILE"
        fi
        
        echo ""
        
    done  # End genotype file loop
    
    echo ""
    
done  # End chromosome loop

#==========================================
# Merge Chunks by Chromosome
#==========================================

if [ "$MERGE_CHUNKS" = "yes" ]; then
    echo "=========================================="
    echo "Merging Chunk Results by Chromosome"
    echo "=========================================="
    echo ""
    
    MERGED_COUNT=0
    MERGE_FAILED=0
    
    for chr in "${CHR_ARRAY[@]}"; do
        # Only merge if this chromosome had chunks
        if [ "$${CHR_HAS_CHUNKS[$$chr]:-0}" -eq 1 ]; then
            echo "Chromosome $chr:"
            if merge_chunk_results "$$chr" "$$OUTPUT_DIR" "$KEEP_CHUNK_FILES"; then
                MERGED_COUNT=$((MERGED_COUNT + 1))
            else
                MERGE_FAILED=$((MERGE_FAILED + 1))
            fi
            echo ""
        fi
    done
    
    if [ $MERGED_COUNT -gt 0 ]; then
        echo "Summary:"
        echo "  Chromosomes merged: $MERGED_COUNT"
        if [ $MERGE_FAILED -gt 0 ]; then
            echo "  Merge failures: $MERGE_FAILED"
        fi
        echo ""
        
        # Update summary file
        {
            echo ""
            echo "=========================================="
            echo "Chunk Merging Summary"
            echo "=========================================="
            echo "Chromosomes merged: $MERGED_COUNT"
            echo "Merge failures: $MERGE_FAILED"
            echo "Keep chunk files: $KEEP_CHUNK_FILES"
            echo ""
        } >> "$SUMMARY_FILE"
    else
        echo "No chunks to merge (all chromosomes have single files)"
        echo ""
    fi
fi

#==========================================
# Final Summary
#==========================================

{
    echo ""
    echo "=========================================="
    echo "Overall Summary"
    echo "=========================================="
    echo "Total jobs: $TOTAL_JOBS"
    echo "Successful: $SUCCESSFUL_JOBS"
    echo "Failed: $FAILED_JOBS"
    echo ""
    if [ "$MERGE_CHUNKS" = "yes" ]; then
        echo "Chunk merging: Enabled"
        echo "  Chromosomes merged: ${MERGED_COUNT:-0}"
        echo "  Keep chunk files: $KEEP_CHUNK_FILES"
        echo ""
    fi
    echo "Results directory: $OUTPUT_DIR"
    echo ""
} | tee -a "$SUMMARY_FILE"

echo "=========================================="
echo "SAIGE-GENE Run Complete"
echo "=========================================="
echo "Total jobs: $TOTAL_JOBS"
echo "Successful: $SUCCESSFUL_JOBS"
echo "Failed: $FAILED_JOBS"
echo ""
if [ "$MERGE_CHUNKS" = "yes" ] && [ "${MERGED_COUNT:-0}" -gt 0 ]; then
    echo "Merged results available:"
    ls -1 "$OUTPUT_DIR"/chr*_combined_results.txt 2>/dev/null | while read f; do
        lines=$(($(wc -l < "$f") - 1))
        echo "  $(basename $f): $lines results"
    done
    echo ""
fi
echo "Log file: $LOG_FILE"
echo "Summary file: $SUMMARY_FILE"
echo ""

if [ $FAILED_JOBS -gt 0 ]; then
    echo "⚠ WARNING: $FAILED_JOBS job(s) failed"
    echo "Check log file for details: $LOG_FILE"
    echo ""
fi

echo "=========================================="
echo "Completed: $(date)"
echo "=========================================="

exit 0

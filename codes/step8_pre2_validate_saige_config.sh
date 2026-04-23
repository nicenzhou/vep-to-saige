#!/bin/bash
# step8_pre2_validate_saige_config.sh
# Validates SAIGE configuration file before running analysis

CONFIG_FILE="$1"

if [ -z "$CONFIG_FILE" ]; then
    echo "Usage: ./step8_pre2_validate_saige_config.sh <config_file>"
    exit 1
fi

if [ ! -f "$CONFIG_FILE" ]; then
    echo "ERROR: Configuration file not found: $CONFIG_FILE"
    exit 1
fi

echo "=========================================="
echo "SAIGE Configuration Validator"
echo "=========================================="
echo ""

# Source the config
source "$CONFIG_FILE"

ERRORS=0
WARNINGS=0

echo "Checking required parameters..."
echo ""

# Check required parameters
check_required() {
    local param_name="$1"
    local param_value="${!param_name:-}"
    
    if [ -z "$param_value" ]; then
        echo "  ✗ ERROR: $param_name is not set"
        ERRORS=$((ERRORS + 1))
        return 1
    else
        echo "  ✓ $param_name is set"
        return 0
    fi
}

# Check file/directory exists
check_exists() {
    local param_name="$1"
    local param_value="${!param_name:-}"
    local type="$2"  # file or directory
    
    if [ -z "$param_value" ]; then
        return 1
    fi
    
    if [ "$type" = "directory" ]; then
        if [ ! -d "$param_value" ]; then
            echo "    ⚠ WARNING: Directory not found: $param_value"
            WARNINGS=$((WARNINGS + 1))
            return 1
        fi
    else
        if [ ! -f "$param_value" ]; then
            echo "    ⚠ WARNING: File not found: $param_value"
            WARNINGS=$((WARNINGS + 1))
            return 1
        fi
    fi
    
    echo "    ✓ Exists: $param_value"
    return 0
}

# Required parameters
check_required GENOTYPE_DIR && check_exists GENOTYPE_DIR directory
check_required OUTPUT_DIR
check_required GMMAT_MODEL && check_exists GMMAT_MODEL file
check_required VARIANCE_RATIO && check_exists VARIANCE_RATIO file

echo ""
echo "Checking group file configuration..."
echo ""

if [ "${GROUP_FILE_BY_CHR:-no}" = "yes" ]; then
    echo "  Using per-chromosome group files"
    GROUP_DIR="${GROUP_DIR:-$GENOTYPE_DIR}"
    check_exists GROUP_DIR directory
    
    # Check if at least one chr group file exists
    if ls "$GROUP_DIR"/chr*_group*.txt 1> /dev/null 2>&1; then
        echo "    ✓ Found group files in directory"
    else
        echo "    ✗ ERROR: No chr*_group*.txt files found in $GROUP_DIR"
        ERRORS=$((ERRORS + 1))
    fi
else
    echo "  Using single group file"
    check_required GROUP_FILE && check_exists GROUP_FILE file
fi

echo ""
echo "Checking module and command configuration..."
echo ""

USE_MODULE="${USE_MODULE:-no}"
echo "  Module system: $USE_MODULE"

if [ "$USE_MODULE" = "yes" ]; then
    # Validate module configuration
    MODULE_NAME="${MODULE_NAME:-}"
    SAIGE_CMD="${SAIGE_CMD:-}"
    
    if [ -z "$MODULE_NAME" ]; then
        echo "    ✗ ERROR: USE_MODULE=yes but MODULE_NAME is not set"
        ERRORS=$((ERRORS + 1))
    else
        echo "  Module name: $MODULE_NAME"
        
        # Check if module command exists
        if command -v module &> /dev/null; then
            echo "    ✓ 'module' command is available"
            
            # Try to check if module exists (different systems have different commands)
            echo "    Checking if module '$MODULE_NAME' is available..."
            
            # Method 1: Try module avail
            if module avail "$MODULE_NAME" 2>&1 | grep -qi "$MODULE_NAME"; then
                echo "    ✓ Module '$MODULE_NAME' is available"
                
                # Try to actually load the module in a subshell
                echo "    Testing module load..."
                if (module load "$MODULE_NAME" 2>&1) > /dev/null 2>&1; then
                    echo "    ✓ Module '$MODULE_NAME' can be loaded"
                    
                    # Test if command becomes available after loading
                    if [ -n "$SAIGE_CMD" ]; then
                        TEST_OUTPUT=$(module load "$MODULE_NAME" 2>&1 && command -v "$SAIGE_CMD" 2>&1)
                        if echo "$TEST_OUTPUT" | grep -q "$SAIGE_CMD"; then
                            echo "    ✓ Command '$SAIGE_CMD' available after loading module"
                        else
                            echo "    ⚠ WARNING: Command '$SAIGE_CMD' not found after loading module"
                            echo "      Module may not provide this command"
                            WARNINGS=$((WARNINGS + 1))
                        fi
                    fi
                else
                    echo "    ✗ ERROR: Failed to load module '$MODULE_NAME'"
                    echo "      Please check module name and availability"
                    ERRORS=$((ERRORS + 1))
                fi
            else
                echo "    ⚠ WARNING: Module '$MODULE_NAME' not found in module avail"
                echo "      Attempting to load anyway (some systems don't list all modules)"
                WARNINGS=$((WARNINGS + 1))
                
                # Still try to load it
                if (module load "$MODULE_NAME" 2>&1) > /dev/null 2>&1; then
                    echo "    ✓ Module '$MODULE_NAME' loaded successfully despite not appearing in avail"
                else
                    echo "    ✗ ERROR: Cannot load module '$MODULE_NAME'"
                    ERRORS=$((ERRORS + 1))
                fi
            fi
        else
            echo "    ✗ ERROR: 'module' command not found"
            echo "      Environment module system is not available on this system"
            echo "      Either install environment modules or set USE_MODULE=no"
            ERRORS=$((ERRORS + 1))
        fi
    fi
    
    if [ -z "$SAIGE_CMD" ]; then
        echo "    ✗ ERROR: USE_MODULE=yes but SAIGE_CMD is not set"
        ERRORS=$((ERRORS + 1))
    else
        echo "  SAIGE command: $SAIGE_CMD"
    fi
    
else
    # Not using modules - validate SAIGE_CMD directly
    SAIGE_CMD="${SAIGE_CMD:-Rscript}"
    echo "  SAIGE command: $SAIGE_CMD"
    
    # Check if it's a path or command name
    if [[ "$SAIGE_CMD" == *"/"* ]]; then
        # It's a path
        if [ -f "$SAIGE_CMD" ]; then
            if [ -x "$SAIGE_CMD" ]; then
                echo "    ✓ Custom SAIGE executable found and is executable"
            else
                echo "    ✗ ERROR: File exists but is not executable: $SAIGE_CMD"
                ERRORS=$((ERRORS + 1))
            fi
        else
            echo "    ✗ ERROR: SAIGE command path not found: $SAIGE_CMD"
            ERRORS=$((ERRORS + 1))
        fi
    else
        # It's a command name - check if it exists in PATH
        if command -v "$SAIGE_CMD" &> /dev/null; then
            SAIGE_PATH=$(command -v "$SAIGE_CMD")
            echo "    ✓ Command '$SAIGE_CMD' found: $SAIGE_PATH"
            
            # For Rscript, also check R version
            if [ "$SAIGE_CMD" = "Rscript" ]; then
                R_VERSION=$(Rscript --version 2>&1 | head -n1)
                echo "    R version: $R_VERSION"
            fi
        else
            echo "    ✗ ERROR: Command '$SAIGE_CMD' not found in PATH"
            echo "      Either install it, provide full path, or use module system"
            ERRORS=$((ERRORS + 1))
        fi
    fi
fi

# Additional check: If using Rscript, verify SAIGE package can be loaded
if [ "$ERRORS" -eq 0 ]; then
    ACTUAL_CMD="$SAIGE_CMD"
    
    # Determine the actual command to test
    if [ "$USE_MODULE" = "yes" ]; then
        # Will need to load module first
        TEST_CMD="module load $MODULE_NAME 2>&1 && $SAIGE_CMD"
    else
        TEST_CMD="$SAIGE_CMD"
    fi
    
    # Only test if it's Rscript-based
    if [[ "$ACTUAL_CMD" =~ [Rr]script ]] || [ "$ACTUAL_CMD" = "Rscript" ]; then
        echo ""
        echo "  Testing SAIGE package availability..."
        
        TEST_SCRIPT='library(SAIGE); cat("SAIGE package loaded successfully\n")'
        
        if [ "$USE_MODULE" = "yes" ]; then
            TEST_OUTPUT=$(module load "$MODULE_NAME" 2>&1 && echo "$TEST_SCRIPT" | $SAIGE_CMD --vanilla - 2>&1)
        else
            TEST_OUTPUT=$(echo "$TEST_SCRIPT" | $SAIGE_CMD --vanilla - 2>&1)
        fi
        
        if echo "$TEST_OUTPUT" | grep -q "SAIGE package loaded successfully"; then
            echo "    ✓ SAIGE R package is available"
        else
            echo "    ⚠ WARNING: Could not verify SAIGE R package"
            echo "      Make sure SAIGE is installed in R"
            if echo "$TEST_OUTPUT" | grep -qi "error"; then
                echo "      Error: $(echo "$TEST_OUTPUT" | grep -i error | head -n1)"
            fi
            WARNINGS=$((WARNINGS + 1))
        fi
    fi
fi

echo ""
echo "Checking format and naming..."
echo ""

INPUT_FORMAT="${INPUT_FORMAT:-bgen}"
echo "  Input format: $INPUT_FORMAT"

case "$INPUT_FORMAT" in
    bgen|bfile|pgen|vcf)
        echo "    ✓ Valid format"
        ;;
    *)
        echo "    ✗ ERROR: Invalid INPUT_FORMAT: $INPUT_FORMAT"
        echo "      Must be: bgen, bfile, pgen, or vcf"
        ERRORS=$((ERRORS + 1))
        ;;
esac

CHR_PREFIX="${CHR_PREFIX:-chr}"
CHR_PADDING="${CHR_PADDING:-auto}"
CHUNK_PATTERN="${CHUNK_PATTERN:-chunk}"
CHUNKED_INPUT="${CHUNKED_INPUT:-no}"

echo "  Chromosome prefix: '$CHR_PREFIX'"
echo "  Chromosome padding: $CHR_PADDING"
echo "  Chunked input: $CHUNKED_INPUT"

if [ "$CHUNKED_INPUT" = "yes" ]; then
    echo "  Chunk pattern: $CHUNK_PATTERN"
fi

if [[ ! "$CHR_PADDING" =~ ^(auto|yes|no)$ ]]; then
    echo "    ⚠ WARNING: CHR_PADDING should be 'auto', 'yes', or 'no'"
    WARNINGS=$((WARNINGS + 1))
fi

echo ""
echo "Checking processing parameters..."
echo ""

THREADS="${THREADS:-4}"
echo "  Threads: $THREADS"
if ! [[ "$THREADS" =~ ^[0-9]+$ ]]; then
    echo "    ✗ ERROR: THREADS must be a positive integer"
    ERRORS=$((ERRORS + 1))
fi

CHROMOSOMES="${CHROMOSOMES:-1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}"
echo "  Chromosomes: $CHROMOSOMES"

# Check if CHROMOSOMES is properly quoted in config
if grep -E "^CHROMOSOMES=[^\"'].*," "$CONFIG_FILE" 2>/dev/null | grep -qv '^#'; then
    echo "    ⚠ WARNING: CHROMOSOMES may not be properly quoted in config"
    echo "      Values with commas should be quoted"
    echo "      Correct: CHROMOSOMES=\"1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22\""
    WARNINGS=$((WARNINGS + 1))
fi

MERGE_CHUNKS="${MERGE_CHUNKS:-no}"
KEEP_CHUNK_FILES="${KEEP_CHUNK_FILES:-yes}"
echo "  Merge chunks: $MERGE_CHUNKS"
echo "  Keep chunk files: $KEEP_CHUNK_FILES"

ALLELE_ORDER="${ALLELE_ORDER:-alt-first}"
echo "  Allele order: $ALLELE_ORDER"

if [[ ! "$ALLELE_ORDER" =~ ^(alt-first|ref-first)$ ]]; then
    echo "    ⚠ WARNING: ALLELE_ORDER should be 'alt-first' or 'ref-first'"
    WARNINGS=$((WARNINGS + 1))
fi

echo ""
echo "Checking SAIGE parameters..."
echo ""

LOCO="${LOCO:-FALSE}"
echo "  LOCO: $LOCO"
if [[ ! "$LOCO" =~ ^(TRUE|FALSE)$ ]]; then
    echo "    ⚠ WARNING: LOCO should be TRUE or FALSE"
    WARNINGS=$((WARNINGS + 1))
fi

echo "  minMAF: ${minMAF:-0}"
echo "  minMAC: ${minMAC:-1}"
echo "  maxMAF_in_groupTest: ${maxMAF_in_groupTest:-0.0001,0.001,0.01}"
echo "  annotation_in_groupTest: ${annotation_in_groupTest:-lof,missense:lof:missense:synonymous}"

# Check for problematic separators in annotation_in_groupTest
if [ -n "${annotation_in_groupTest:-}" ]; then
    # Check if it contains semicolons (old format)
    if [[ "$annotation_in_groupTest" == *";"* ]]; then
        echo "    ⚠ WARNING: annotation_in_groupTest contains semicolons (;)"
        echo "      SAIGE expects colons (:) as separators"
        echo "      Current value: $annotation_in_groupTest"
        echo "      Expected format: lof,missense:lof:missense:synonymous"
        WARNINGS=$((WARNINGS + 1))
    fi
    
    # Check if value is properly quoted in config file
    if grep -E "^annotation_in_groupTest=[^\"'].*[:;]" "$CONFIG_FILE" 2>/dev/null | grep -qv '^#'; then
        echo "    ⚠ WARNING: annotation_in_groupTest may not be properly quoted in config"
        echo "      Values with colons or semicolons should be quoted"
        echo "      Correct: annotation_in_groupTest=\"lof,missense:lof:missense:synonymous\""
        WARNINGS=$((WARNINGS + 1))
    fi
fi

# Check for problematic separators in maxMAF_in_groupTest
if [ -n "${maxMAF_in_groupTest:-}" ]; then
    # Check if value is properly quoted in config file (contains commas)
    if grep -E "^maxMAF_in_groupTest=[^\"'].*," "$CONFIG_FILE" 2>/dev/null | grep -qv '^#'; then
        echo "    ⚠ WARNING: maxMAF_in_groupTest may not be properly quoted in config"
        echo "      Values with commas should be quoted"
        echo "      Correct: maxMAF_in_groupTest=\"0.0001,0.001,0.01\""
        WARNINGS=$((WARNINGS + 1))
    fi
fi

# Check r_corr (stored with underscore, converted to r.corr for SAIGE)
r_corr="${r_corr:-0}"
echo "  r_corr: $r_corr (converted to --r.corr for SAIGE)"
if [[ ! "$r_corr" =~ ^[01]$ ]]; then
    echo "    ⚠ WARNING: r_corr should be 0 (SKAT-O) or 1 (Burden)"
    WARNINGS=$((WARNINGS + 1))
fi

# Check for old r.corr format and give helpful error
if grep -E "^r\.corr=" "$CONFIG_FILE" 2>/dev/null | grep -qv '^#'; then
    echo "    ✗ ERROR: Found 'r.corr=' in config file"
    echo "      Please use 'r_corr=' instead (bash variables cannot contain dots)"
    echo "      The script will automatically convert r_corr to --r.corr for SAIGE"
    ERRORS=$((ERRORS + 1))
fi

is_Firth_beta="${is_Firth_beta:-FALSE}"
echo "  is_Firth_beta: $is_Firth_beta"
if [[ ! "$is_Firth_beta" =~ ^(TRUE|FALSE)$ ]]; then
    echo "    ⚠ WARNING: is_Firth_beta should be TRUE or FALSE"
    WARNINGS=$((WARNINGS + 1))
fi

if [ "$is_Firth_beta" = "TRUE" ]; then
    pCutoffforFirth="${pCutoffforFirth:-0.01}"
    echo "  pCutoffforFirth: $pCutoffforFirth"
fi

# Check other boolean parameters
is_output_markerList_in_groupTest="${is_output_markerList_in_groupTest:-TRUE}"
echo "  is_output_markerList_in_groupTest: $is_output_markerList_in_groupTest"
if [[ ! "$is_output_markerList_in_groupTest" =~ ^(TRUE|FALSE)$ ]]; then
    echo "    ⚠ WARNING: is_output_markerList_in_groupTest should be TRUE or FALSE"
    WARNINGS=$((WARNINGS + 1))
fi

# Check imputed data parameters
is_imputed_data="${is_imputed_data:-FALSE}"
if [ "$is_imputed_data" = "TRUE" ]; then
    echo ""
    echo "Imputed data parameters:"
    echo "  is_imputed_data: TRUE"
    echo "  minInfo: ${minInfo:-not set}"
    echo "  dosage_zerod_cutoff: ${dosage_zerod_cutoff:-not set}"
    echo "  dosage_zerod_MAC_cutoff: ${dosage_zerod_MAC_cutoff:-not set}"
    echo "  impute_method: ${impute_method:-not set}"
fi

echo ""
echo "Checking genotype files..."
echo ""

# Count genotype files
if [ -d "$GENOTYPE_DIR" ]; then
    
    # Determine file extension based on format
    case "$INPUT_FORMAT" in
        bgen)
            FILE_EXT="bgen"
            ;;
        bfile)
            FILE_EXT="bed"
            ;;
        pgen)
            FILE_EXT="pgen"
            ;;
        vcf)
            FILE_EXT="vcf.gz"
            ;;
    esac
    
    if [ "$CHUNKED_INPUT" = "yes" ]; then
        # Check for chunked files
        GENO_COUNT=$(ls -1 "$GENOTYPE_DIR"/${CHR_PREFIX}*_genes_${CHUNK_PATTERN}*.${FILE_EXT} 2>/dev/null | wc -l)
        echo "  Found $GENO_COUNT chunked genotype files"
        echo "  Pattern: ${CHR_PREFIX}*_genes_${CHUNK_PATTERN}*.${FILE_EXT}"
        
        if [ "$GENO_COUNT" -eq 0 ]; then
            echo "    ⚠ WARNING: No chunked genotype files found"
            WARNINGS=$((WARNINGS + 1))
        else
            # Show sample files
            echo ""
            echo "  Sample chunked files found:"
            ls -1 "$GENOTYPE_DIR"/${CHR_PREFIX}*_genes_${CHUNK_PATTERN}*.${FILE_EXT} 2>/dev/null | head -3 | while read file; do
                echo "    - $(basename "$file")"
            done
            
            # Count chunks per chromosome - FIXED ORDERING
            echo ""
            echo "  Chunks per chromosome:"
            
            # Parse CHROMOSOMES into array, maintaining order
            IFS=',' read -ra CHR_ARRAY <<< "$CHROMOSOMES"
            
            # Sort the array numerically to ensure proper order (1,2,3...10,11...22)
            # Use sort with version sort to handle numbers properly
            SORTED_CHRS=$(printf '%s\n' "${CHR_ARRAY[@]}" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//' | sort -V)
            
            while IFS= read -r chr; do
                # Determine chromosome number format
                if [ "$CHR_PADDING" = "yes" ] || ([ "$CHR_PADDING" = "auto" ] && [ ${#chr} -eq 1 ] && [[ "$chr" =~ ^[0-9]$ ]]); then
                    CHR_NUM=$(printf "%02d" "$chr" 2>/dev/null || echo "$chr")
                else
                    CHR_NUM="$chr"
                fi
                
                # Count chunks for this chromosome
                CHUNK_COUNT=$(ls -1 "$GENOTYPE_DIR"/${CHR_PREFIX}${CHR_NUM}_genes_${CHUNK_PATTERN}*.${FILE_EXT} 2>/dev/null | wc -l)
                
                if [ "$CHUNK_COUNT" -gt 0 ]; then
                    printf "    Chr%-3s: %3d chunks\n" "$chr" "$CHUNK_COUNT"
                fi
            done <<< "$SORTED_CHRS"
            
            # Also check for any chromosomes not in the list
            echo ""
            echo "  Checking for additional chromosomes not in CHROMOSOMES list..."
            ALL_CHRS=$(ls -1 "$GENOTYPE_DIR"/${CHR_PREFIX}*_genes_${CHUNK_PATTERN}*.${FILE_EXT} 2>/dev/null | \
                       sed -n "s/.*${CHR_PREFIX}\([^_]*\)_genes_${CHUNK_PATTERN}.*/\1/p" | \
                       sort -u | sort -V)
            
            EXTRA_CHRS=""
            while IFS= read -r found_chr; do
                # Check if this chr is in CHROMOSOMES
                if ! echo "$CHROMOSOMES" | grep -qw "$found_chr"; then
                    # Also check without leading zero
                    found_chr_no_zero=$(echo "$found_chr" | sed 's/^0*//')
                    if ! echo "$CHROMOSOMES" | grep -qw "$found_chr_no_zero"; then
                        EXTRA_CHRS="$EXTRA_CHRS $found_chr"
                    fi
                fi
            done <<< "$ALL_CHRS"
            
            if [ -n "$EXTRA_CHRS" ]; then
                echo "    ℹ INFO: Found files for chromosomes not in CHROMOSOMES list:$EXTRA_CHRS"
                echo "      These will be skipped during analysis"
            else
                echo "    ✓ All found chromosomes are in CHROMOSOMES list"
            fi
        fi
    else
        # Check for non-chunked files
        GENO_COUNT=$(ls -1 "$GENOTYPE_DIR"/${CHR_PREFIX}*_genes.${FILE_EXT} 2>/dev/null | wc -l)
        echo "  Found $GENO_COUNT non-chunked genotype files"
        echo "  Pattern: ${CHR_PREFIX}*_genes.${FILE_EXT}"
        
        if [ "$GENO_COUNT" -eq 0 ]; then
            # Also check for generic pattern
            GENO_COUNT=$(ls -1 "$GENOTYPE_DIR"/*.${FILE_EXT} 2>/dev/null | wc -l)
            echo "  Found $GENO_COUNT total .$FILE_EXT files in directory"
            
            if [ "$GENO_COUNT" -eq 0 ]; then
                echo "    ⚠ WARNING: No genotype files found"
                WARNINGS=$((WARNINGS + 1))
            fi
        else
            # Show sample files
            echo ""
            echo "  Sample files found:"
            ls -1 "$GENOTYPE_DIR"/${CHR_PREFIX}*_genes.${FILE_EXT} 2>/dev/null | head -5 | while read file; do
                echo "    - $(basename "$file")"
            done
            
            # List all chromosomes found - FIXED ORDERING
            echo ""
            echo "  Chromosomes with genotype files:"
            ALL_CHR_FILES=$(ls -1 "$GENOTYPE_DIR"/${CHR_PREFIX}*_genes.${FILE_EXT} 2>/dev/null | \
                           sed -n "s/.*${CHR_PREFIX}\([^_]*\)_genes\.${FILE_EXT}/\1/p" | \
                           sort -V)
            
            while IFS= read -r found_chr; do
                # Check if in CHROMOSOMES list
                if echo "$CHROMOSOMES" | grep -qw "$found_chr"; then
                    STATUS="✓"
                else
                    # Check without leading zero
                    found_chr_no_zero=$(echo "$found_chr" | sed 's/^0*//')
                    if echo "$CHROMOSOMES" | grep -qw "$found_chr_no_zero"; then
                        STATUS="✓"
                    else
                        STATUS="⊗ (not in CHROMOSOMES list)"
                    fi
                fi
                printf "    Chr%-3s %s\n" "$found_chr" "$STATUS"
            done <<< "$ALL_CHR_FILES"
        fi
    fi
else
    echo "  ⚠ WARNING: Cannot check genotype files (directory not found)"
fi

# Check for companion files (bgen needs .bgi, vcf needs .tbi, etc.)
if [ -d "$GENOTYPE_DIR" ] && [ "$GENO_COUNT" -gt 0 ]; then
    echo ""
    echo "Checking companion index files..."
    
    case "$INPUT_FORMAT" in
        bgen)
            INDEX_EXT="bgen.bgi"
            INDEX_COUNT=$(ls -1 "$GENOTYPE_DIR"/*.bgen.bgi 2>/dev/null | wc -l)
            INDEX_NAME="BGEN index (.bgi)"
            ;;
        vcf)
            INDEX_EXT="vcf.gz.tbi"
            INDEX_COUNT=$(ls -1 "$GENOTYPE_DIR"/*.vcf.gz.tbi 2>/dev/null | wc -l)
            INDEX_NAME="Tabix index (.tbi)"
            ;;
        bfile)
            INDEX_EXT="bim"
            INDEX_COUNT=$(ls -1 "$GENOTYPE_DIR"/*.bim 2>/dev/null | wc -l)
            INDEX_NAME="PLINK .bim"
            
            # Also check for .fam files
            FAM_COUNT=$(ls -1 "$GENOTYPE_DIR"/*.fam 2>/dev/null | wc -l)
            echo "  Found $FAM_COUNT .fam files"
            if [ "$FAM_COUNT" -lt "$GENO_COUNT" ]; then
                echo "    ⚠ WARNING: Fewer .fam files than .bed files"
                WARNINGS=$((WARNINGS + 1))
            fi
            ;;
        pgen)
            INDEX_EXT="pvar"
            INDEX_COUNT=$(ls -1 "$GENOTYPE_DIR"/*.pvar 2>/dev/null | wc -l)
            INDEX_NAME="PLINK2 .pvar"
            
            # Also check for .psam files
            PSAM_COUNT=$(ls -1 "$GENOTYPE_DIR"/*.psam 2>/dev/null | wc -l)
            echo "  Found $PSAM_COUNT .psam files"
            if [ "$PSAM_COUNT" -lt "$GENO_COUNT" ]; then
                echo "    ⚠ WARNING: Fewer .psam files than .pgen files"
                WARNINGS=$((WARNINGS + 1))
            fi
            ;;
    esac
    
    if [ -n "$INDEX_EXT" ]; then
        echo "  Found $INDEX_COUNT $INDEX_NAME files"
        
        if [ "$INDEX_COUNT" -lt "$GENO_COUNT" ]; then
            echo "    ⚠ WARNING: Fewer index files than genotype files"
            echo "      Expected: $GENO_COUNT, Found: $INDEX_COUNT"
            WARNINGS=$((WARNINGS + 1))
        elif [ "$INDEX_COUNT" -gt "$GENO_COUNT" ]; then
            echo "    ℹ INFO: More index files than genotype files (may include backups)"
        else
            echo "    ✓ Index file count matches genotype file count"
        fi
    fi
fi

# Check sample file if bgen format
if [ "$INPUT_FORMAT" = "bgen" ] && [ -d "$GENOTYPE_DIR" ]; then
    echo ""
    echo "Checking for BGEN sample files..."
    SAMPLE_COUNT=$(ls -1 "$GENOTYPE_DIR"/*.sample 2>/dev/null | wc -l)
    if [ "$SAMPLE_COUNT" -gt 0 ]; then
        echo "  Found $SAMPLE_COUNT .sample files"
        echo "    ✓ Sample files present"
    else
        echo "  No .sample files found"
        echo "    ℹ INFO: BGEN files may contain sample information internally"
    fi
fi

echo ""
echo "=========================================="
echo "Validation Summary"
echo "=========================================="
echo ""

# Configuration summary
echo "Configuration Overview:"
echo "  Input: $INPUT_FORMAT files from $GENOTYPE_DIR"
echo "  Output: $OUTPUT_DIR"
echo "  Chromosomes: $(echo "$CHROMOSOMES" | tr ',' ' ')"
echo "  Module system: $USE_MODULE"
if [ "$USE_MODULE" = "yes" ]; then
    echo "  Module: $MODULE_NAME"
fi
echo "  Command: $SAIGE_CMD"
echo "  Test type: $([ "$r_corr" -eq 0 ] && echo 'SKAT-O (r.corr=0)' || echo 'Burden (r.corr=1)')"
echo "  Chunked: $CHUNKED_INPUT"
if [ "$CHUNKED_INPUT" = "yes" ]; then
    echo "  Merge chunks: $MERGE_CHUNKS"
fi
echo ""

echo "Results:"
echo "  Errors: $ERRORS"
echo "  Warnings: $WARNINGS"
echo ""

if [ $ERRORS -eq 0 ]; then
    if [ $WARNINGS -eq 0 ]; then
        echo "✓✓✓ Configuration is valid with no warnings ✓✓✓"
    else
        echo "✓ Configuration is valid (with $WARNINGS warning(s))"
        echo ""
        echo "Note: Warnings indicate potential issues but won't prevent running."
        echo "      Review warnings above and proceed if acceptable."
    fi
    echo ""
    echo "Ready to run:"
    echo "  ./step8_run_saige_gene_tests.sh $CONFIG_FILE"
    echo ""
    exit 0
else
    echo "✗✗✗ Configuration has $ERRORS error(s) that must be fixed ✗✗✗"
    if [ $WARNINGS -gt 0 ]; then
        echo "    Also has $WARNINGS warning(s)"
    fi
    echo ""
    exit 1
fi

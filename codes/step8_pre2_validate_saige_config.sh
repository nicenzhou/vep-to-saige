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
echo "  annotation_in_groupTest: ${annotation_in_groupTest:-lof,missense;lof,missense;lof;synonymous}"

# Check r_corr (stored with underscore, converted to r.corr for SAIGE)
r_corr="${r_corr:-0}"
echo "  r_corr: $r_corr (converted to --r.corr for SAIGE)"
if [[ ! "$r_corr" =~ ^[01]$ ]]; then
    echo "    ⚠ WARNING: r_corr should be 0 (SKAT-O) or 1 (Burden)"
    WARNINGS=$((WARNINGS + 1))
fi

# Check for old r.corr format and give helpful error
if grep -q "^r\.corr=" "$CONFIG_FILE" 2>/dev/null; then
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
            
            # Count chunks per chromosome
            echo ""
            echo "  Chunks per chromosome:"
            for chr in {1..22} X Y; do
                if [ "$CHR_PADDING" = "yes" ] || ([ "$CHR_PADDING" = "auto" ] && [[ "$chr" =~ ^[0-9]$ ]]); then
                    CHR_NUM=$(printf "%02d" "$chr" 2>/dev/null || echo "$chr")
                else
                    CHR_NUM="$chr"
                fi
                
                CHUNK_COUNT=$(ls -1 "$GENOTYPE_DIR"/${CHR_PREFIX}${CHR_NUM}_genes_${CHUNK_PATTERN}*.${FILE_EXT} 2>/dev/null | wc -l)
                
                if [ "$CHUNK_COUNT" -gt 0 ]; then
                    echo "    Chr${chr}: $CHUNK_COUNT chunks"
                fi
            done
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
            ;;
        vcf)
            INDEX_EXT="vcf.gz.tbi"
            INDEX_COUNT=$(ls -1 "$GENOTYPE_DIR"/*.vcf.gz.tbi 2>/dev/null | wc -l)
            ;;
        bfile)
            INDEX_EXT="bim"
            INDEX_COUNT=$(ls -1 "$GENOTYPE_DIR"/*.bim 2>/dev/null | wc -l)
            ;;
        pgen)
            INDEX_EXT="pvar"
            INDEX_COUNT=$(ls -1 "$GENOTYPE_DIR"/*.pvar 2>/dev/null | wc -l)
            ;;
    esac
    
    if [ -n "$INDEX_EXT" ]; then
        echo "  Found $INDEX_COUNT .$INDEX_EXT files"
        
        if [ "$INDEX_COUNT" -lt "$GENO_COUNT" ]; then
            echo "    ⚠ WARNING: Fewer index files than genotype files"
            WARNINGS=$((WARNINGS + 1))
        fi
    fi
fi

echo ""
echo "=========================================="
echo "Validation Summary"
echo "=========================================="
echo "Errors: $ERRORS"
echo "Warnings: $WARNINGS"
echo ""

if [ $ERRORS -eq 0 ]; then
    if [ $WARNINGS -eq 0 ]; then
        echo "✓ Configuration is valid with no warnings"
    else
        echo "✓ Configuration is valid (with $WARNINGS warning(s))"
    fi
    echo ""
    echo "Ready to run:"
    echo "  ./step8_run_saige_gene_tests.sh $CONFIG_FILE"
    echo ""
    exit 0
else
    echo "✗ Configuration has $ERRORS error(s) that must be fixed"
    if [ $WARNINGS -gt 0 ]; then
        echo "  Also has $WARNINGS warning(s)"
    fi
    echo ""
    exit 1
fi

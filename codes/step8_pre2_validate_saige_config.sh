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
echo "  Chromosome prefix: '$CHR_PREFIX'"
echo "  Chromosome padding: $CHR_PADDING"

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

r_corr="${r.corr:-0}"
echo "  r.corr: $r_corr"
if [[ ! "$r_corr" =~ ^[01]$ ]]; then
    echo "    ⚠ WARNING: r.corr should be 0 or 1"
    WARNINGS=$((WARNINGS + 1))
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
    case "$INPUT_FORMAT" in
        bgen)
            GENO_COUNT=$(ls -1 "$GENOTYPE_DIR"/*.bgen 2>/dev/null | wc -l)
            ;;
        bfile)
            GENO_COUNT=$(ls -1 "$GENOTYPE_DIR"/*.bed 2>/dev/null | wc -l)
            ;;
        pgen)
            GENO_COUNT=$(ls -1 "$GENOTYPE_DIR"/*.pgen 2>/dev/null | wc -l)
            ;;
        vcf)
            GENO_COUNT=$(ls -1 "$GENOTYPE_DIR"/*.vcf.gz 2>/dev/null | wc -l)
            ;;
    esac
    
    echo "  Found $GENO_COUNT genotype files in $GENOTYPE_DIR"
    
    if [ "$GENO_COUNT" -eq 0 ]; then
        echo "    ⚠ WARNING: No genotype files found"
        WARNINGS=$((WARNINGS + 1))
    fi
else
    echo "  ⚠ WARNING: Cannot check genotype files (directory not found)"
fi

echo ""
echo "=========================================="
echo "Validation Summary"
echo "=========================================="
echo "Errors: $ERRORS"
echo "Warnings: $WARNINGS"
echo ""

if [ $ERRORS -eq 0 ]; then
    echo "✓ Configuration appears valid"
    echo ""
    echo "Ready to run:"
    echo "  ./step8_run_saige_gene_tests.sh $CONFIG_FILE"
    echo ""
    exit 0
else
    echo "✗ Configuration has errors that must be fixed"
    echo ""
    exit 1
fi

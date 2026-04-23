#!/bin/bash
# step8_pre1_build_saige_config.sh
# Interactive SAIGE configuration file builder

echo "=========================================="
echo "SAIGE Configuration Builder"
echo "=========================================="
echo ""
echo "This wizard will help you create a SAIGE configuration file."
echo ""

read -p "Output configuration file name [saige_config.txt]: " CONFIG_FILE
CONFIG_FILE="${CONFIG_FILE:-saige_config.txt}"

if [ -f "$CONFIG_FILE" ]; then
    read -p "File exists. Overwrite? (yes/no): " OVERWRITE
    if [ "$OVERWRITE" != "yes" ]; then
        echo "Aborted."
        exit 1
    fi
fi

echo ""
echo "=== Required Parameters ==="
echo ""

read -p "Genotype directory: " GENOTYPE_DIR
read -p "Output directory: " OUTPUT_DIR
read -p "GMMAT model file: " GMMAT_MODEL
read -p "Variance ratio file: " VARIANCE_RATIO

echo ""
read -p "Use separate group file per chromosome? (yes/no) [yes]: " GROUP_BY_CHR
GROUP_BY_CHR="${GROUP_BY_CHR:-yes}"

if [ "$GROUP_BY_CHR" = "yes" ]; then
    read -p "Group files directory: " GROUP_DIR
    GROUP_FILE=""
else
    read -p "Single group file path: " GROUP_FILE
    GROUP_DIR=""
fi

echo ""
echo "=== Format and Processing ==="
echo ""

read -p "Input format (bgen/bfile/pgen/vcf) [bgen]: " INPUT_FORMAT
INPUT_FORMAT="${INPUT_FORMAT:-bgen}"

echo ""
echo "=== Module Configuration ==="
echo ""

read -p "Use environment module system? (yes/no) [no]: " USE_MODULE
USE_MODULE="${USE_MODULE:-no}"

if [ "$USE_MODULE" = "yes" ]; then
    read -p "Module name to load [SAIGE]: " MODULE_NAME
    MODULE_NAME="${MODULE_NAME:-SAIGE}"
    
    read -p "SAIGE command after loading module [saige]: " SAIGE_CMD
    SAIGE_CMD="${SAIGE_CMD:-saige}"
else
    read -p "SAIGE command or path [Rscript]: " SAIGE_CMD
    SAIGE_CMD="${SAIGE_CMD:-Rscript}"
    MODULE_NAME=""
fi

echo ""
echo "=== File Naming Configuration ==="
echo ""

read -p "Chromosome prefix [chr]: " CHR_PREFIX
CHR_PREFIX="${CHR_PREFIX:-chr}"

read -p "Chromosome padding (auto/yes/no) [auto]: " CHR_PADDING
CHR_PADDING="${CHR_PADDING:-auto}"

read -p "Number of threads [8]: " THREADS
THREADS="${THREADS:-8}"

read -p "Chromosomes to process [1-22]: " CHROMOSOMES
CHROMOSOMES="${CHROMOSOMES:-1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}"

echo ""
echo "=== Chunked File Configuration ==="
echo ""

read -p "Are genotype files chunked? (yes/no) [no]: " CHUNKED_INPUT
CHUNKED_INPUT="${CHUNKED_INPUT:-no}"

if [ "$CHUNKED_INPUT" = "yes" ]; then
    read -p "Chunk pattern in filename [chunk]: " CHUNK_PATTERN
    CHUNK_PATTERN="${CHUNK_PATTERN:-chunk}"
    echo "  Expected file naming: $${CHR_PREFIX}*_genes_$${CHUNK_PATTERN}*"
else
    CHUNK_PATTERN="chunk"
fi

echo ""
read -p "Merge chunk results? (yes/no) [yes]: " MERGE_CHUNKS
MERGE_CHUNKS="${MERGE_CHUNKS:-yes}"

if [ "$MERGE_CHUNKS" = "yes" ]; then
    read -p "Keep individual chunk files? (yes/no) [no]: " KEEP_CHUNK_FILES
    KEEP_CHUNK_FILES="${KEEP_CHUNK_FILES:-no}"
else
    KEEP_CHUNK_FILES="yes"
fi

echo ""
echo "=== Allele Configuration ==="
echo ""

read -p "Allele order (alt-first/ref-first) [alt-first]: " ALLELE_ORDER
ALLELE_ORDER="${ALLELE_ORDER:-alt-first}"

echo ""
echo "=== SAIGE Parameters ==="
echo ""

read -p "Use LOCO? (TRUE/FALSE) [TRUE]: " LOCO
LOCO="${LOCO:-TRUE}"

read -p "Minimum MAF [0]: " minMAF
minMAF="${minMAF:-0}"

read -p "Minimum MAC [1]: " minMAC
minMAC="${minMAC:-1}"

read -p "MAF cutoffs for gene tests [0.0001,0.001,0.01]: " maxMAF_in_groupTest
maxMAF_in_groupTest="${maxMAF_in_groupTest:-0.0001,0.001,0.01}"

echo ""
echo "NOTE: Annotation categories use colon (:) to separate different masks"
echo "      Example: lof,missense:lof:missense:synonymous creates 4 masks"
read -p "Annotation categories [lof,missense:lof:missense:synonymous]: " annotation_in_groupTest
annotation_in_groupTest="${annotation_in_groupTest:-lof,missense:lof:missense:synonymous}"

echo ""
echo "NOTE: For test type, enter 0 for SKAT-O or 1 for Burden test"
echo "      (internally stored as r_corr, converted to --r.corr for SAIGE)"
read -p "Test type - r_corr (0=SKAT-O, 1=Burden) [0]: " r_corr
r_corr="${r_corr:-0}"

read -p "Use Firth correction? (TRUE/FALSE) [TRUE]: " is_Firth_beta
is_Firth_beta="${is_Firth_beta:-TRUE}"

if [ "$is_Firth_beta" = "TRUE" ]; then
    read -p "P-value cutoff for Firth [0.01]: " pCutoffforFirth
    pCutoffforFirth="${pCutoffforFirth:-0.01}"
else
    pCutoffforFirth=""
fi

echo ""
echo "=== Advanced Options (press Enter to skip) ==="
echo ""

read -p "Is data imputed? (TRUE/FALSE) [FALSE]: " is_imputed_data
is_imputed_data="${is_imputed_data:-FALSE}"

if [ "$is_imputed_data" = "TRUE" ]; then
    read -p "Minimum Info score [0.3]: " minInfo
    minInfo="${minInfo:-0.3}"
    
    read -p "Dosage zero cutoff [0.2]: " dosage_zerod_cutoff
    dosage_zerod_cutoff="${dosage_zerod_cutoff:-0.2}"
    
    read -p "Dosage zero MAC cutoff [10]: " dosage_zerod_MAC_cutoff
    dosage_zerod_MAC_cutoff="${dosage_zerod_MAC_cutoff:-10}"
    
    read -p "Impute method (minor/mean/best_guess) [mean]: " impute_method
    impute_method="${impute_method:-mean}"
else
    minInfo=""
    dosage_zerod_cutoff=""
    dosage_zerod_MAC_cutoff=""
    impute_method=""
fi

read -p "Output marker lists? (TRUE/FALSE) [TRUE]: " is_output_markerList_in_groupTest
is_output_markerList_in_groupTest="${is_output_markerList_in_groupTest:-TRUE}"

read -p "Maximum missing rate [0.15]: " maxMissing
maxMissing="${maxMissing:-0.15}"

read -p "SPA cutoff [2]: " SPAcutoff
SPAcutoff="${SPAcutoff:-2}"

echo ""
echo "=== Writing Configuration File ==="
echo ""

cat > "$CONFIG_FILE" << EOF
#!/bin/bash
# SAIGE-GENE Configuration File
# Generated: $(date)
# Generated by: step8_pre1_build_saige_config.sh

#==========================================================================
# REQUIRED PARAMETERS
#==========================================================================

# Input/Output Directories
GENOTYPE_DIR="$GENOTYPE_DIR"
OUTPUT_DIR="$OUTPUT_DIR"
GMMAT_MODEL="$GMMAT_MODEL"
VARIANCE_RATIO="$VARIANCE_RATIO"

# Group File Configuration
EOF

if [ "$GROUP_BY_CHR" = "yes" ]; then
    cat >> "$CONFIG_FILE" << EOF
GROUP_FILE_BY_CHR=yes
GROUP_DIR="$GROUP_DIR"
EOF
else
    cat >> "$CONFIG_FILE" << EOF
GROUP_FILE="$GROUP_FILE"
GROUP_FILE_BY_CHR=no
EOF
fi

cat >> "$CONFIG_FILE" << EOF

#==========================================================================
# FORMAT AND PROCESSING
#==========================================================================

# Input format
INPUT_FORMAT=$INPUT_FORMAT

# Module Configuration
USE_MODULE=$USE_MODULE
EOF

if [ "$USE_MODULE" = "yes" ]; then
    cat >> "$CONFIG_FILE" << EOF
MODULE_NAME=$MODULE_NAME
SAIGE_CMD=$SAIGE_CMD
EOF
else
    cat >> "$CONFIG_FILE" << EOF
SAIGE_CMD=$SAIGE_CMD
EOF
fi

cat >> "$CONFIG_FILE" << EOF

# File naming format
CHR_PREFIX=$CHR_PREFIX
CHR_PADDING=$CHR_PADDING

# Processing control
THREADS=$THREADS
CHROMOSOMES="$CHROMOSOMES"

# Chunked file configuration
CHUNKED_INPUT=$CHUNKED_INPUT
CHUNK_PATTERN=$CHUNK_PATTERN

# Result merging
MERGE_CHUNKS=$MERGE_CHUNKS
KEEP_CHUNK_FILES=$KEEP_CHUNK_FILES

# Allele configuration
ALLELE_ORDER=$ALLELE_ORDER

#==========================================================================
# SAIGE PARAMETERS
#==========================================================================

# LOCO and filtering
LOCO=$LOCO
minMAF=$minMAF
minMAC=$minMAC
EOF

if [ -n "$maxMissing" ]; then
    echo "maxMissing=$$maxMissing" >> "$$CONFIG_FILE"
fi

cat >> "$CONFIG_FILE" << EOF

# Gene-based test parameters
maxMAF_in_groupTest="$maxMAF_in_groupTest"
annotation_in_groupTest="$annotation_in_groupTest"
r_corr=$r_corr

# Effect size estimation
is_Firth_beta=$is_Firth_beta
EOF

if [ -n "$pCutoffforFirth" ]; then
    echo "pCutoffforFirth=$$pCutoffforFirth" >> "$$CONFIG_FILE"
fi

# Add SPA cutoff
if [ -n "$SPAcutoff" ]; then
    cat >> "$CONFIG_FILE" << EOF

# SPA parameters
SPAcutoff=$SPAcutoff
EOF
fi

# Add imputed data parameters if applicable
if [ "$is_imputed_data" = "TRUE" ]; then
    cat >> "$CONFIG_FILE" << EOF

# Imputed data parameters
is_imputed_data=$is_imputed_data
EOF
    [ -n "$$minInfo" ] && echo "minInfo=$$minInfo" >> "$CONFIG_FILE"
    [ -n "$$dosage_zerod_cutoff" ] && echo "dosage_zerod_cutoff=$$dosage_zerod_cutoff" >> "$CONFIG_FILE"
    [ -n "$$dosage_zerod_MAC_cutoff" ] && echo "dosage_zerod_MAC_cutoff=$$dosage_zerod_MAC_cutoff" >> "$CONFIG_FILE"
    [ -n "$$impute_method" ] && echo "impute_method=$$impute_method" >> "$CONFIG_FILE"
fi

# Add output parameters
cat >> "$CONFIG_FILE" << EOF

# Output control
is_output_markerList_in_groupTest=$is_output_markerList_in_groupTest

#==========================================================================
# ADDITIONAL OPTIONS
#==========================================================================

# Uncomment and modify as needed:
# MACCutoff_to_CollapseUltraRare=10
# minGroupMAC_in_BurdenTest=5
# weights_beta="1,25"
# is_single_in_groupTest=TRUE
# is_no_weight_in_groupTest=FALSE
# markers_per_chunk=10000
# groups_per_chunk=100

# IMPORTANT NOTES:
# - r_corr is used in config (underscore) and converted to --r.corr for SAIGE
#   (bash variable names cannot contain dots)
# - annotation_in_groupTest uses colons (:) to separate different annotation masks
#   Example: "lof,missense:lof:missense:synonymous" creates 4 masks:
#     1. lof,missense (combined)
#     2. lof (only)
#     3. missense (only)
#     4. synonymous (only)
# - Values with special characters (colons, semicolons, commas) MUST be quoted
# - CHROMOSOMES parameter MUST be quoted when containing commas

# MODULE CONFIGURATION EXAMPLES:
# Example 1: Use module system
#   USE_MODULE=yes
#   MODULE_NAME=SAIGE
#   SAIGE_CMD=saige
#
# Example 2: Use Rscript directly (no module)
#   USE_MODULE=no
#   SAIGE_CMD=Rscript
#
# Example 3: Custom SAIGE installation
#   USE_MODULE=no
#   SAIGE_CMD=/path/to/saige/bin/saige

EOF

echo "✓ Configuration file created: $CONFIG_FILE"
echo ""

# Validate the configuration
echo "=== Validating Configuration ==="
echo ""

ERRORS=0
WARNINGS=0

# Check required paths
if [ ! -d "$GENOTYPE_DIR" ]; then
    echo "⚠ WARNING: Genotype directory not found: $GENOTYPE_DIR"
    ERRORS=$((ERRORS + 1))
fi

if [ ! -f "$GMMAT_MODEL" ]; then
    echo "⚠ WARNING: GMMAT model file not found: $GMMAT_MODEL"
    ERRORS=$((ERRORS + 1))
fi

if [ ! -f "$VARIANCE_RATIO" ]; then
    echo "⚠ WARNING: Variance ratio file not found: $VARIANCE_RATIO"
    ERRORS=$((ERRORS + 1))
fi

if [ "$GROUP_BY_CHR" = "yes" ]; then
    if [ ! -d "$GROUP_DIR" ]; then
        echo "⚠ WARNING: Group directory not found: $GROUP_DIR"
        ERRORS=$((ERRORS + 1))
    fi
else
    if [ ! -f "$GROUP_FILE" ]; then
        echo "⚠ WARNING: Group file not found: $GROUP_FILE"
        ERRORS=$((ERRORS + 1))
    fi
fi

# Validate input format
if [[ ! "$$INPUT_FORMAT" =~ ^(bgen|bfile|pgen|vcf)$$ ]]; then
    echo "⚠ WARNING: Invalid input format '$INPUT_FORMAT' (should be bgen, bfile, pgen, or vcf)"
    ERRORS=$((ERRORS + 1))
fi

# Validate allele order
if [[ ! "$$ALLELE_ORDER" =~ ^(alt-first|ref-first)$$ ]]; then
    echo "⚠ WARNING: Invalid allele order '$ALLELE_ORDER' (should be alt-first or ref-first)"
    ERRORS=$((ERRORS + 1))
fi

# Validate r_corr
if [[ ! "$$r_corr" =~ ^[01]$$ ]]; then
    echo "⚠ WARNING: r_corr should be 0 (SKAT-O) or 1(Burden), got: $r_corr"
    WARNINGS=$((WARNINGS + 1))
fi

# Validate module configuration
if [ "$USE_MODULE" = "yes" ]; then
    if [ -z "$MODULE_NAME" ]; then
        echo "⚠ WARNING: USE_MODULE=yes but MODULE_NAME is empty"
        ERRORS=$((ERRORS + 1))
    fi
    if [ -z "$SAIGE_CMD" ]; then
        echo "⚠ WARNING: USE_MODULE=yes but SAIGE_CMD is empty"
        ERRORS=$((ERRORS + 1))
    fi
    echo "ℹ Module configuration:"
    echo "  - Will load module: $MODULE_NAME"
    echo "  - Will use command: $SAIGE_CMD"
else
    if [ -z "$SAIGE_CMD" ]; then
        echo "⚠ WARNING: SAIGE_CMD is empty"
        ERRORS=$((ERRORS + 1))
    else
        echo "ℹ Module configuration:"
        echo "  - Module system: disabled"
        echo "  - Will use command: $SAIGE_CMD"
        
        # Check if command exists (only if it's not a path)
        if [[ ! "$SAIGE_CMD" =~ / ]]; then
            if ! command -v "$SAIGE_CMD" &> /dev/null; then
                echo "⚠ WARNING: Command '$SAIGE_CMD' not found in PATH"
                WARNINGS=$((WARNINGS + 1))
            fi
        elif [ ! -f "$SAIGE_CMD" ] && [ ! -x "$SAIGE_CMD" ]; then
            echo "⚠ WARNING: SAIGE_CMD path '$SAIGE_CMD' not found or not executable"
            WARNINGS=$((WARNINGS + 1))
        fi
    fi
fi

# Validate chromosome padding
if [[ ! "$CHR_PADDING" =~ ^(auto|yes|no)$ ]]; then
    echo "⚠ WARNING: Invalid CHR_PADDING '$CHR_PADDING' (should be auto, yes, or no)"
    WARNINGS=$((WARNINGS + 1))
fi

# Check for chunked files if specified
if [ "$CHUNKED_INPUT" = "yes" ] && [ -d "$GENOTYPE_DIR" ]; then
    CHUNK_COUNT=$(ls "$GENOTYPE_DIR"/${CHR_PREFIX}*_genes_${CHUNK_PATTERN}* 2>/dev/null | wc -l)
    if [ $CHUNK_COUNT -eq 0 ]; then
        echo "⚠ WARNING: No chunked files found matching pattern: ${CHR_PREFIX}*_genes_${CHUNK_PATTERN}*"
        ERRORS=$((ERRORS + 1))
    else
        echo "✓ Found $CHUNK_COUNT chunked genotype files"
    fi
fi

# Validate numeric parameters
if ! [[ "$THREADS" =~ ^[0-9]+$ ]]; then
    echo "⚠ WARNING: THREADS should be a positive integer, got: $THREADS"
    WARNINGS=$((WARNINGS + 1))
fi

if ! [[ "$minMAF" =~ ^[0-9]*\.?[0-9]+$ ]]; then
    echo "⚠ WARNING: minMAF should be numeric, got: $minMAF"
    WARNINGS=$((WARNINGS + 1))
fi

if ! [[ "$minMAC" =~ ^[0-9]+$ ]]; then
    echo "⚠ WARNING: minMAC should be a positive integer, got: $minMAC"
    WARNINGS=$((WARNINGS + 1))
fi

# Validate boolean parameters
for param in LOCO is_Firth_beta is_output_markerList_in_groupTest is_imputed_data; do
    eval value=\$$param
    if [ -n "$value" ] && [[ ! "$value" =~ ^(TRUE|FALSE)$ ]]; then
        echo "⚠ WARNING: $param should be TRUE or FALSE, got: $value"
        WARNINGS=$((WARNINGS + 1))
    fi
done

echo ""

# Summary
if [ $ERRORS -eq 0 ] && [ $WARNINGS -eq 0 ]; then
    echo "=========================================="
    echo "✓ Validation passed successfully!"
    echo "=========================================="
    echo ""
    echo "Configuration Summary:"
    echo "  Input format: $INPUT_FORMAT"
    echo "  Module system: $USE_MODULE"
    if [ "$USE_MODULE" = "yes" ]; then
        echo "  Module name: $MODULE_NAME"
    fi
    echo "  SAIGE command: $SAIGE_CMD"
    echo "  Chromosome prefix: $CHR_PREFIX"
    echo "  Chromosome padding: $CHR_PADDING"
    echo "  Chunked input: $CHUNKED_INPUT"
    echo "  Merge chunks: $MERGE_CHUNKS"
    echo "  Allele order: $ALLELE_ORDER"
    echo "  Test type: $([ "$r_corr" -eq 0 ] && echo 'SKAT-O' || echo 'Burden')"
    echo ""
    echo "Next steps:"
    echo "  1. Review configuration: cat $CONFIG_FILE"
    echo "  2. Run full validation: ./validate_saige_config.sh $CONFIG_FILE"
    echo "  3. Run SAIGE analysis: ./step8_run_saige_gene_tests.sh $CONFIG_FILE"
elif [ $ERRORS -eq 0 ]; then
    echo "=========================================="
    echo "✓ Validation passed with $WARNINGS warning(s)"
    echo "=========================================="
    echo ""
    echo "The configuration file was created but has warnings."
    echo "Please review the warnings above."
    echo ""
    echo "Next steps:"
    echo "  1. Review configuration: cat $CONFIG_FILE"
    echo "  2. Run full validation: ./validate_saige_config.sh $CONFIG_FILE"
    echo "  3. Fix warnings if needed, then run: ./step8_run_saige_gene_tests.sh $CONFIG_FILE"
else
    echo "=========================================="
    echo "⚠ Validation found issues"
    echo "=========================================="
    echo ""
    echo "Found $ERRORS error(s) and $WARNINGS warning(s)"
    echo "Please review and fix the issues above before running analysis."
    echo ""
    echo "Next steps:"
    echo "  1. Review configuration: cat $CONFIG_FILE"
    echo "  2. Fix the errors listed above"
    echo "  3. Run validation: ./validate_saige_config.sh $CONFIG_FILE"
    echo "  4. When ready: ./step8_run_saige_gene_tests.sh $CONFIG_FILE"
fi

echo ""
echo "=========================================="
echo "Configuration File Contents"
echo "=========================================="
cat "$CONFIG_FILE"
echo "=========================================="

exit 0

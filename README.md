# vep-to-saige: VEP Annotations to SAIGE Gene Groups Pipeline

Convert VEP-annotated VCF files into SAIGE-compatible gene group files for rare variant burden testing.

---

## Clone and setup

```bash
git clone https://github.com/nicenzhou/vep-to-saige.git
cd vep-to-saige
sed -i 's/\r$//' *.sh  # Fix Windows line endings if needed
chmod +x *.sh
```

---

## Quick Start: Full Pipeline

```
# Steps 1-3: VEP to SAIGE groups
./step1_vep_ann_clean.sh input.vcf.gz chr1_anno.txt 4 vep_lof
./step2_create_gene_groups.sh chr1_anno.txt chr1_groups.txt all keepall
./step3_merge_and_validate_groups.sh all_genes.txt .

# Steps 4-5: Gene coordinates and region matching
# Download gene coordinates (requires internet)
./step4_download_gene_coords.sh 115 GRCh38 gene_coords

# Or use pre-downloaded coordinates from GitHub
tar -xzf gene_coords_ensembl115.tar.gz

# Match genes to groups and create region files
./step5_match_genes_to_groups.sh gene_coords all_genes_groups.txt plink_regions 10 no yes

# Step 6: Extract genotypes (supports chunking)
./step6_extract_genotypes_plink2.sh /data/genotypes plink_regions gene_bfiles 16 bgen auto chr 0

# Step 7: Verify extraction quality
./step7_verify_extraction.sh gene_bfiles plink_regions verification_report.txt bgen

# Step 8: Run SAIGE-GENE association tests
# Create configuration interactively
./step8_pre1_build_saige_config.sh

# Or create manually
cat > saige_config.txt << EOF
GENOTYPE_DIR=gene_bfiles
OUTPUT_DIR=saige_results
GMMAT_MODEL=step1_null_model.rda
VARIANCE_RATIO=step1_variance_ratio.txt
GROUP_FILE_BY_CHR=yes
GROUP_DIR=.
INPUT_FORMAT=bgen
THREADS=16
CHROMOSOMES=1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22
MERGE_CHUNKS=yes
LOCO=TRUE
maxMAF_in_groupTest=0.0001,0.001,0.01
annotation_in_groupTest=lof,missense;lof;missense;synonymous
is_Firth_beta=TRUE
r.corr=0
EOF

# Validate configuration
./step8_pre2_validate_saige_config.sh saige_config.txt

# Run SAIGE analysis
./step8_run_saige_gene_tests.sh saige_config.txt

# Analyze results
./saige_results/step9_analyze_results.sh
```

---

## Quick Start: Alternative Workflows

```
#--------------------------------------------------------------------------
# Workflow 1: Basic gene-based analysis (no genotype extraction)
#--------------------------------------------------------------------------

# Steps 1-3: Create gene groups from VEP
./step1_vep_ann_clean.sh input.vcf.gz chr1_anno.txt 4 vep_lof
./step2_create_gene_groups.sh chr1_anno.txt chr1_groups.txt all keepall
./step3_merge_and_validate_groups.sh all_genes.txt .

# Use groups directly with your genotypes in SAIGE
# (Skip steps 4-7 if you already have organized genotype files)

#--------------------------------------------------------------------------
# Workflow 2: Large WGS dataset with chunking
#--------------------------------------------------------------------------

# Steps 1-3: Create gene groups
./step1_vep_ann_clean.sh wgs.vcf.gz wgs_anno.txt 16 vep_lof
./step2_create_gene_groups.sh wgs_anno.txt wgs_groups.txt all keepall
./step3_merge_and_validate_groups.sh all_genes.txt .

# Steps 4-5: Gene coordinates with regeneration for missing genes
tar -xzf gene_coords_ensembl115.tar.gz
./step5_match_genes_to_groups.sh gene_coords all_genes_groups.txt plink_regions 50 yes yes

# Step 6: Extract with chunking (20 genes per chunk)
./step6_extract_genotypes_plink2.sh /data/bgen plink_regions gene_bfiles 16 bgen bgen imputed_chr 20

# Step 7: Verify chunked extraction
./step7_verify_extraction.sh gene_bfiles plink_regions qc_report.txt bgen

# Step 8: Run SAIGE with chunk merging
cat > wgs_config.txt << EOF
GENOTYPE_DIR=gene_bfiles
OUTPUT_DIR=wgs_saige_results
GMMAT_MODEL=wgs_null_model.rda
VARIANCE_RATIO=wgs_variance_ratio.txt
GROUP_FILE_BY_CHR=yes
GROUP_DIR=.
INPUT_FORMAT=bgen
CHR_PREFIX=imputed_chr
CHR_PADDING=yes
THREADS=32
CHROMOSOMES=1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22
MERGE_CHUNKS=yes
KEEP_CHUNK_FILES=no
LOCO=TRUE
is_imputed_data=TRUE
minInfo=0.3
maxMAF_in_groupTest=0.0001,0.001,0.01
is_Firth_beta=TRUE
EOF

./step8_run_saige_gene_tests.sh wgs_config.txt
./saige_results/step9_analyze_results.sh

#--------------------------------------------------------------------------
# Workflow 3: Quick test on chromosome 22
#--------------------------------------------------------------------------

# Steps 1-3: Process chr22 only
./step1_vep_ann_clean.sh chr22.vcf.gz chr22_anno.txt 4 vep_lof
./step2_create_gene_groups.sh chr22_anno.txt chr22_groups.txt all keepall

# Steps 4-5: Match genes
tar -xzf gene_coords_ensembl115.tar.gz
./step5_match_genes_to_groups.sh gene_coords chr22_groups.txt test_regions 10

# Step 6: Extract genotypes
./step6_extract_genotypes_plink2.sh /data/bfiles test_regions test_output 4 bfile auto chr 0

# Step 7: Verify
./step7_verify_extraction.sh test_output test_regions test_qc.txt bfile

# Step 8: Quick SAIGE test
cat > test_config.txt << EOF
GENOTYPE_DIR=test_output
OUTPUT_DIR=test_results
GMMAT_MODEL=step1_null_model.rda
VARIANCE_RATIO=step1_variance_ratio.txt
GROUP_FILE=chr22_groups.txt
GROUP_FILE_BY_CHR=no
INPUT_FORMAT=bfile
THREADS=4
CHROMOSOMES=22
MERGE_CHUNKS=no
LOCO=FALSE
minMAC=5
maxMAF_in_groupTest=0.01
annotation_in_groupTest=lof,missense
r.corr=0
EOF

./step8_run_saige_gene_tests.sh test_config.txt
./saige_results/step9_analyze_results.sh

#--------------------------------------------------------------------------
# Workflow 4: Using interactive configuration builder
#--------------------------------------------------------------------------

# Steps 1-5: Same as above
# ...

# Step 6: Extract genotypes
./step6_extract_genotypes_plink2.sh /data/bgen plink_regions gene_bfiles 16 bgen bgen chr 0

# Step 7: Verify
./step7_verify_extraction.sh gene_bfiles plink_regions

# Step 8: Interactive configuration
./step8_pre1_build_saige_config.sh
# Follow prompts to create configuration

# Validate before running
./step8_pre2_validate_saige_config.sh saige_config.txt

# Run analysis
./step8_run_saige_gene_tests.sh saige_config.txt

# Analyze results
./saige_results/step9_analyze_results.sh
```

--- 

## Expected output structure 

```
# vep-to-saige/
# ├── chr*_anno.txt              # Step 1 output
# ├── chr*_groups.txt            # Step 2 output
# ├── all_genes_groups.txt       # Step 3 output
# ├── gene_coords/               # Step 4 output
# │   ├── chr1_genes.txt
# │   ├── chr2_genes.txt
# │   └── ...
# ├── plink_regions/             # Step 5 output
# │   ├── chr*_regions.txt
# │   ├── chr*_gene_list.txt
# │   ├── matched_genes.txt
# │   ├── missing_genes.txt
# │   ├── regenerated_genes.txt  # If force_regen=yes
# │   └── matching_summary.txt
# ├── gene_bfiles/               # Step 6 output
# │   ├── chr*_genes.bed/bim/fam # Or .bgen, .pgen, .vcf.gz
# │   ├── chr*_genes_chunk*.bed  # If chunked
# │   ├── extraction_summary.txt
# │   └── extraction.log
# ├── saige_results/             # Step 8 output
# │   ├── chr*_combined_results.txt  # If MERGE_CHUNKS=yes
# │   ├── chr*_chunk*_results.txt    # If MERGE_CHUNKS=no
# │   ├── saige_run_summary.txt
# │   ├── saige_run.log
# │   └── step9_analyze_results.sh
# └── *.sh                       # Pipeline scripts
```

---

## Requirements

- ```bcftools``` (VCF processing)
- ```awk``` (GNU awk recommended)
- ```plink2a``` (for genotype extraction, Step 6)
- Bash 4.0+

---

## Pipeline Overview

### Step 1: VEP Annotation Extraction

Extract and classify variants by functional impact.

**Usage:**
```bash
./step1_vep_ann_clean.sh <input.vcf.gz> <output.txt> [threads] [lof_mode]
```

**LoF Modes:**

| Mode | Definition | Use Case |
|------|------------|----------|
| ```loftee_only``` | LOFTEE HC only | High confidence LoF |
| ```vep_lof``` (default) | VEP consequences | Standard analysis |
| ```any``` | LOFTEE HC OR VEP LoF | Maximum sensitivity |

**VEP LoF Consequences** (```vep_lof``` mode):
- ```frameshift_variant``` - Frameshift indels
- ```stop_gained``` - Nonsense mutations
- ```stop_lost``` - Stop codon loss
- ```splice_acceptor_variant``` - Splice acceptor disruption
- ```splice_donor_variant``` - Splice donor disruption
- ```inframe_deletion``` - In-frame deletions
- ```inframe_insertion``` - In-frame insertions

**How VEP LoF Works:**
1. VEP assigns consequence terms based on variant position/impact
2. Script matches consequences against LoF patterns using regex
3. Matched variants → ```group='lof'```
4. Non-matched coding variants → ```missense``` or ```synonymous```

**Examples:**
```bash
# Default VEP LoF mode
./step1_vep_ann_clean.sh chr1.vcf.gz chr1_anno.txt 4

# LOFTEE high-confidence only
./step1_vep_ann_clean.sh chr1.vcf.gz chr1_anno.txt 4 loftee_only

# Batch processing
for chr in {1..22}; do
  ./step1_vep_ann_clean.sh chr$${chr}.vcf.gz chr$${chr}_anno.txt 4 vep_lof
done
```

**Output:** Tab-delimited file with columns:
```
Variant  Gene  Consequence  LoF  LoF_flags  Group  Total_Anno  NON_CAN_SPLICE_Count
```

---

### Step 2: Gene Group Creation

Convert annotations to SAIGE format with quality filtering and duplicate handling.

**Usage:**
```bash
./step2_create_gene_groups.sh <input_anno.txt> <output_groups.txt> [filter] [priority]
```

**Annotation Filters:**
- ```all``` - All variants (default)
- ```lof``` - LoF only
- ```missense``` - Missense only
- ```synonymous``` - Synonymous only
- ```lof+missense``` - LoF and missense

**Priority Modes:**

SAIGE cannot take one single variant with two or more different annotations. This will keep the annotation based on the priority you desire.

| Mode | Behavior |
|------|----------|
| ```lof,missense,synonymous``` | Keep highest priority; lof>missense>synonymous (default) |
| ```missense,lof,synonymous``` | Prioritize missense; missense>lof>synonymous |
| ```keepall``` | Keep all annotations separately |

**Examples:**
```bash
# Default priority (lof > missense > synonymous)
./step2_create_gene_groups.sh chr1_anno.txt chr1_groups.txt

# Keep all annotations
./step2_create_gene_groups.sh chr1_anno.txt chr1_groups.txt all keepall

# LoF variants only
./step2_create_gene_groups.sh chr1_anno.txt chr1_lof.txt lof

# Batch with keepall
for chr in {01..22}; do
  ./step2_create_gene_groups.sh chr$${chr}_anno.txt chr$${chr}_groups.txt all keepall
done
```

**Output:** Space-delimited SAIGE format:
```
GENE1 var variant1 variant2 variant3
GENE1 anno lof missense lof
```

**Quality Filters:**
- Missing gene symbols removed
- LoF with all ```NON_CAN_SPLICE``` removed
- Data errors (```NON_CAN_SPLICE_Count > Total_Anno```) removed
- Exact duplicates removed
- Priority-based duplicate resolution

---

### Step 3: Merge and Validate

Merge chromosome files into a genome-wide file.

**Usage:**
```bash
./step3_merge_and_validate_groups.sh <output_file> [input_dir] [pattern]
```

**Examples:**
```bash
# Auto-detect in current directory
./step3_merge_and_validate_groups.sh all_genes_groups.txt

# Specify input directory
./step3_merge_and_validate_groups.sh all_genes.txt ./results

# Custom pattern
./step3_merge_and_validate_groups.sh output.txt ./data "custom_chr*_groups.txt"
```

**Features:**
- Auto-detects ```chr01``` or ```chr1``` naming
- Merges in chromosomal order (1→22)
- Reports gene counts and annotation breakdown
- No validation (assumes Step 2 QC complete)

---

### Step 4: Download Gene Coordinates (Require Internet)

Download Ensembl gene coordinates to match with SAIGE group files.

**Usage:**
```bash
./step4_download_gene_coords.sh [ensembl_version] [genome_build] [output_dir]
```

**Arguments:**
- ```ensembl_version``` - Ensembl release version (default: 115)
- ```genome_build``` - GRCh38 or GRCh37 (default: GRCh38)
- ```output_dir``` - Output directory (default: gene_coords)

**Examples:**
```bash
# Download latest Ensembl coordinates
./step4_download_gene_coords.sh 115 GRCh38
```

**Output:**
- Per-chromosome gene coordinate files (```chr1_genes.txt```, ```chr2_genes.txt```, etc.)

---

### Step 5: Match Genes to Groups

Match gene coordinates with SAIGE group files and create PLINK2a-compatible region files in BED format.

**Usage:**
```
./step5_match_genes_to_groups.sh <gene_coords_dir> <group_file> <output_dir> [buffer_kb] [force_regen] [merge_regen]
```

**Arguments:**
- ```gene_coords_dir``` - Directory with ```chr*_genes.txt``` files (from gene coordinate extraction)
- ```group_file``` - SAIGE group file (from Step 3)
- ```output_dir``` - Output directory for region files (default: plink_regions)
- ```buffer_kb``` - Buffer around genes in kb (default: 10)
- ```force_regen``` - Force regenerate missing genes: yes/no (default: no)
- ```merge_regen``` - Merge regenerated genes with matched: yes/no (default: yes)

**Force Regenerate Mode:**

When ```force_regen=yes```, missing genes are regenerated using variant positions from the group file:

- **Start position**: First variant position - buffer_kb - 10kb (soft buffer)
- **End position**: Last variant position + buffer_kb + 10kb (soft buffer)
- **Total buffer**: User-defined buffer + 10kb safety margin on each side
- **Chromosome extraction**: Automatically extracts chromosome from variant format (chr:pos:ref:alt)
- Missing genes are saved to ```regenerated_genes.txt``` file

**Merge Options:**

- **merge_regen=yes** (default): 
  - Regenerated genes are merged into main ```chr*_regions.txt``` files
  - All regions sorted by genomic coordinates (start position)
  - Gene lists updated to maintain genomic order
  - Ideal for streamlined downstream analysis

- **merge_regen=no**: 
  - Regenerated genes saved to separate ```chr*_regions_recovered.txt``` files
  - Allows independent quality control of regenerated vs. matched genes
  - Enables separate statistical analyses if needed
  - Useful for comparing annotation-based vs. variant-based gene boundaries

**Examples:**
```
# Extract gene coordinates (prerequisite)
tar -xzf gene_coords_ensembl111.tar.gz

# Basic usage - match genes with 10kb buffer (default settings)
./step5_match_genes_to_groups.sh gene_coords all_genes_groups.txt plink_regions

# Custom 50kb buffer without regeneration
./step5_match_genes_to_groups.sh gene_coords all_genes_groups.txt plink_regions 50

# Force regenerate missing genes and merge with matched genes
./step5_match_genes_to_groups.sh gene_coords all_genes_groups.txt plink_regions 10 yes yes

# Force regenerate missing genes but keep them separate for QC
./step5_match_genes_to_groups.sh gene_coords all_genes_groups.txt plink_regions 10 yes no

# Custom 50kb buffer with regeneration and merge
./step5_match_genes_to_groups.sh gene_coords all_genes_groups.txt plink_regions 50 yes yes

# Large buffer for regulatory regions, keep separate
./step5_match_genes_to_groups.sh gene_coords all_genes_groups.txt plink_regions 100 yes no
```

**Output Files (merge_regen=yes):**
- ```chr*_regions.txt``` - BED format region files per chromosome (matched + regenerated, sorted by coordinates)
- ```chr*_gene_list.txt``` - Gene symbols for each chromosome (genomic order)
- ```matched_genes.txt``` - All matched genes with full coordinates
- ```regenerated_genes.txt``` - Regenerated genes with coordinates (if force_regen=yes)
- ```missing_genes.txt``` - Genes still not found after regeneration attempt
- ```matching_summary.txt``` - Detailed summary statistics with match rates
- ```matching.log``` - Complete processing log with timestamps

**Output Files (merge_regen=no):**
- ```chr*_regions.txt``` - BED format region files (matched genes only)
- ```chr*_regions_recovered.txt``` - BED format region files (regenerated genes only)
- ```chr*_gene_list.txt``` - Gene symbols (matched genes only)
- ```chr*_gene_list_recovered.txt``` - Gene symbols (regenerated genes only)
- ```matched_genes.txt``` - Matched genes with coordinates
- ```regenerated_genes.txt``` - Regenerated genes with coordinates
- ```missing_genes.txt``` - Genes still not found
- ```matching_summary.txt``` - Detailed summary statistics
- ```matching.log``` - Complete processing log

**Region File Format (BED):**

BED format with 0-based start, 1-based end:
```
# chr  start(0-based)  end(1-based)  gene_symbol  gene_id
1      11868           14409         DDX11L1      ENSG00000223972
1      14403           29570         WASH7P       ENSG00000227232
1      69090           70008         OR4F5        ENSG00000186092
```

For regenerated genes:
```
# chr  start(0-based)  end(1-based)  gene_symbol  source
2      12345           67890         GENE1        REGENERATED
3      98765           123456        GENE2        REGENERATED
```

---

### Step 6: Extract Genotypes with PLINK2a

Extract genotypes for gene regions using PLINK2a with optional chunking support.
Creates output files per chromosome (and per chunk if chunking enabled).

**Usage:**
```
./step6_extract_genotypes_plink2.sh <input_dir> <regions_dir> <output_dir> [threads] [output_format] [input_format] [input_prefix] [chunk_size]
```

**Arguments:**
- ```input_dir``` - Directory with chromosome genotype files
- ```regions_dir``` - Directory with chr*_regions.txt files (from step5_match_genes_to_groups.sh)
- ```output_dir``` - Output directory for extracted files (default: gene_bfiles)
- ```threads``` - Number of threads for PLINK2a (default: 4)
- ```output_format``` - Output format: bfile, pgen, vcf, bgen (default: bfile)
- ```input_format``` - Input format: auto, bfile, pgen, vcf, bgen (default: auto)
- ```input_prefix``` - Input file prefix pattern (default: chr)
- ```chunk_size``` - Number of genes per chunk (default: 0 = no chunking)

**Chunking Feature:**

Chunk size controls how genes are grouped into output files:

- ```chunk_size = 0``` : Extract all genes in one file per chromosome (default)
  - Output: ```chr1_genes.bed```, ```chr2_genes.bed```, etc.
  - Best for: Standard workflows, smaller datasets

- ```chunk_size = 20``` : Extract 20 genes per chunk, multiple files per chromosome
  - Output: ```chr1_genes_chunk1.bed```, ```chr1_genes_chunk2.bed```, etc.
  - Best for: Memory constraints, parallel processing

- ```chunk_size = 50``` : Extract 50 genes per chunk
  - Output: ```chr1_genes_chunk1.bed```, ```chr1_genes_chunk2.bed```, etc.
  - Best for: Balanced processing, moderate datasets

**Input File Prefix Examples:**
- ```chr``` → chr1.bed, chr2.bed, ... (default)
- ```wgs_chr``` → wgs_chr1.bed, wgs_chr2.bed, ...
- ```imputed_chr``` → imputed_chr1.bed, imputed_chr2.bed, ...
- ```""``` → 1.bed, 2.bed, ... (empty string for no prefix)

**Supported Input Formats (auto-detected):**
- **BFILE**: ```.bed/.bim/.fam``` (PLINK1 binary)
- **PGEN**: ```.pgen/.pvar/.psam``` (PLINK2 binary)
- **VCF**: ```.vcf.gz``` or ```.vcf``` (VCF format)
- **BGEN**: ```.bgen``` with ```.sample``` file (BGEN format)

**Supported Output Formats:**
- **bfile**: PLINK1 binary (```.bed/.bim/.fam```)
- **pgen**: PLINK2 binary (```.pgen/.pvar/.psam```)
- **vcf**: VCF format (```.vcf.gz``` with ```.tbi``` index)
- **bgen**: BGEN v1.2 (```.bgen``` with ```.sample```)

**Examples:**

```
# Default: no chunking, all genes in one file per chromosome
./step6_extract_genotypes_plink2.sh /data/genotypes plink_regions gene_bfiles

# Extract 20 genes per chunk for memory efficiency
./step6_extract_genotypes_plink2.sh /data/genotypes plink_regions gene_bfiles 16 bfile auto chr 20

# Extract 50 genes per chunk with custom prefix
./step6_extract_genotypes_plink2.sh /data/genotypes plink_regions gene_bfiles 16 bfile auto wgs_chr 50

# BGEN output with 30 genes per chunk (for SAIGE)
./step6_extract_genotypes_plink2.sh /data/bgen plink_regions gene_bfiles 16 bgen bgen imputed_chr 30

# VCF input, BFILE output, no chunking
./step6_extract_genotypes_plink2.sh /data/vcf plink_regions gene_bfiles 8 bfile vcf chr 0

# Auto-detect input, PGEN output, 100 genes per chunk
./step6_extract_genotypes_plink2.sh /data/genotypes plink_regions gene_bfiles 8 pgen auto chr 100

# No prefix (files named 1.bed, 2.bed, etc.), 25 genes per chunk
./step6_extract_genotypes_plink2.sh /data/genotypes plink_regions gene_bfiles 16 bfile auto "" 25
```

**Output Files (no chunking, chunk_size=0):**
- ```chr1_genes.bed/bim/fam``` - All genes on chr1
- ```chr2_genes.bed/bim/fam``` - All genes on chr2
- ```extraction_summary.txt``` - Summary of extraction
- ```extraction.log``` - Detailed processing log

**Output Files (with chunking, e.g., chunk_size=20):**
- ```chr1_genes_chunk1.bed/bim/fam``` - Genes 1-20 on chr1
- ```chr1_genes_chunk2.bed/bim/fam``` - Genes 21-40 on chr1
- ```chr1_genes_chunk3.bed/bim/fam``` - Genes 41-60 on chr1
- ```chr2_genes_chunk1.bed/bim/fam``` - Genes 1-20 on chr2
- ```extraction_summary.txt``` - Summary with chunk information
- ```extraction.log``` - Detailed processing log

---

### Step 7: Verify Extraction

Verify extraction quality and completeness.

**Usage:**
```bash
./step7_verify_extraction.sh <gene_bfiles_dir> <regions_dir> [output_report]
```

**Examples:**
```bash
# Run verification
./step7_verify_extraction.sh gene_pfiles plink_regions

# Custom report name
./step7_verify_extraction.sh gene_pfiles plink_regions my_qc_report.txt
```

**Output:**
- Verification report with:
  - Variant counts per chromosome
  - Sample counts
  - MAF distribution
  - File sizes
  - Missing genes (if any)

---

### Step 8: Run SAIGE-GENE Association Tests

Run SAIGE-GENE Step 2 association tests on extracted genotype files.
Supports chunked files, multiple formats, and automatic result merging.

**Usage:**
```bash
./step8_run_saige_gene_tests.sh <config_file>
```

**Configuration File:**

Create a configuration file with required and optional parameters.

**Option 1: Interactive Configuration Builder (Recommended)**

```
# Run the interactive configuration builder
./step8_pre1_build_saige_config.sh

# Follow the prompts - key parameters:
# - Output configuration file name [saige_config.txt]: my_config.txt
# - Genotype directory: gene_bfiles
# - Output directory: saige_results
# - GMMAT model file: step1_null_model.rda
# - Variance ratio file: step1_variance_ratio.txt
# - Use separate group file per chromosome? (yes/no) [yes]: yes
# - Group files directory: group_files
# - Input format (bgen/bfile/pgen/vcf) [bgen]: bgen
# - Use environment module system? (yes/no) [no]: yes
# - Module name to load [SAIGE]: SAIGE
# - SAIGE command after loading module [saige]: saige
# - Chromosome prefix [chr]: chr
# - Chromosome padding (auto/yes/no) [auto]: auto
# - Number of threads [8]: 16
# - Chromosomes to process [1-22]: 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22
# - Are genotype files chunked? (yes/no) [no]: yes
# - Chunk pattern in filename [chunk]: chunk
# - Merge chunk results? (yes/no) [yes]: yes
# - Keep individual chunk files? (yes/no) [no]: no
# - Allele order (alt-first/ref-first) [alt-first]: alt-first
# - Use LOCO? (TRUE/FALSE) [TRUE]: TRUE
# - Minimum MAF [0]: 0
# - Minimum MAC [1]: 1
# - MAF cutoffs for gene tests [0.0001,0.001,0.01]: 0.0001,0.001,0.01
# - Annotation categories [lof,missense:lof:missense:synonymous]: lof,missense:lof:missense:synonymous
# - Test type - r_corr (0=SKAT-O, 1=Burden) [0]: 0
# - Use Firth correction? (TRUE/FALSE) [TRUE]: TRUE
# - P-value cutoff for Firth [0.01]: 0.01
# - SPA cutoff [2]: 2
# - Is data imputed? (TRUE/FALSE) [FALSE]: FALSE
# - Output marker lists? (TRUE/FALSE) [TRUE]: TRUE
# - Maximum missing rate [0.15]: 0.15

# The wizard creates and validates the configuration file
# Then run:
./step8_run_saige_gene_tests.sh my_config.txt
```

**Option 2: Manual Configuration File**

```
# Required parameters
GENOTYPE_DIR=/path/to/gene_bfiles
OUTPUT_DIR=/path/to/saige_results
GMMAT_MODEL=/path/to/step1_null_model.rda
VARIANCE_RATIO=/path/to/step1_variance_ratio.txt

# Group file (choose one option)
GROUP_FILE=/path/to/all_genes_groups.txt
GROUP_FILE_BY_CHR=no

# OR use per-chromosome files
# GROUP_FILE_BY_CHR=yes
# GROUP_DIR=/path/to/group_files

# Format and processing
INPUT_FORMAT=bgen
CHR_PREFIX=chr
CHR_PADDING=auto
CHUNKED_INPUT=yes
CHUNK_PATTERN=chunk
THREADS=16
CHROMOSOMES="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22"

# Module configuration (choose one method)
# Method 1: Use module system (recommended for HPC)
USE_MODULE=yes
MODULE_NAME=SAIGE
SAIGE_CMD=saige

# Method 2: Use Rscript directly
# USE_MODULE=no
# SAIGE_CMD=Rscript

# Method 3: Custom SAIGE installation
# USE_MODULE=no
# SAIGE_CMD=/path/to/saige/bin/saige

# Result merging
MERGE_CHUNKS=yes
KEEP_CHUNK_FILES=no

# Allele configuration
ALLELE_ORDER=alt-first

# SAIGE parameters (NOTE: use r_corr not r.corr)
LOCO=TRUE
minMAF=0
minMAC=1
maxMAF_in_groupTest="0.0001,0.001,0.01"
annotation_in_groupTest="lof,missense:lof,missense:lof:synonymous"
r_corr=0
is_Firth_beta=TRUE
pCutoffforFirth=0.01
SPAcutoff=2
is_output_markerList_in_groupTest=TRUE
```

**Examples:**
```bash
# Step 1: Create configuration
./step8_pre1_build_saige_config.sh

# Step 2: Validate configuration
./step8_pre2_validate_saige_config.sh saige_config.txt

# Step 3: Run SAIGE analysis
./step8_run_saige_gene_tests.sh saige_config.txt
```

**Key Configuration Parameters**

**Required Parameters**
- `GENOTYPE_DIR` - Directory with extracted genotype files
- `OUTPUT_DIR` - Output directory for results
- `GMMAT_MODEL` - Null model from SAIGE Step 1
- `VARIANCE_RATIO` - Variance ratio from SAIGE Step 1
- `GROUP_FILE` or `GROUP_FILE_BY_CHR` - Gene group definitions

**Module Configuration (choose one)**
- `USE_MODULE=yes` + `MODULE_NAME=SAIGE` + `SAIGE_CMD=saige` (for HPC with modules)
- `USE_MODULE=no` + `SAIGE_CMD=Rscript` (direct Rscript usage)
- `USE_MODULE=no` + `SAIGE_CMD=/path/to/saige` (custom installation)

**Format Parameters**
- `INPUT_FORMAT` - bgen, bfile, pgen, or vcf
- `CHR_PREFIX` - Prefix in file names (chr, wgs_chr, etc.)
- `CHR_PADDING` - auto/yes/no for chr01 vs chr1
- `CHUNKED_INPUT` - yes/no for chunked genotype files
- `CHUNK_PATTERN` - Pattern to identify chunks (default: chunk)

**Processing Parameters**
- `THREADS` - Number of CPU threads
- `CHROMOSOMES` - Comma-separated list (MUST BE QUOTED)
- `MERGE_CHUNKS` - Combine chunk results (yes/no)
- `KEEP_CHUNK_FILES` - Keep chunks after merge (yes/no)
- `ALLELE_ORDER` - alt-first or ref-first

**SAIGE Parameters**
- `LOCO` - Leave-one-chromosome-out (TRUE/FALSE)
- `minMAF` - Minimum minor allele frequency
- `minMAC` - Minimum minor allele count
- `maxMAF_in_groupTest` - MAF masks (MUST BE QUOTED)
- `annotation_in_groupTest` - Variant annotations (use COLONS, MUST BE QUOTED)
- `is_Firth_beta` - Use Firth correction (TRUE/FALSE)
- `r_corr` - 0 for SKAT-O, 1 for Burden (NOTE: underscore not dot)
- `SPAcutoff` - SPA test threshold
- `pCutoffforFirth` - P-value cutoff for Firth
- `is_output_markerList_in_groupTest` - Output marker lists (TRUE/FALSE)
- `maxMissing` - Maximum missing rate

**Imputed Data Parameters (optional)**
- `is_imputed_data` - TRUE/FALSE
- `minInfo` - Minimum imputation INFO score
- `dosage_zerod_cutoff` - Dosage cutoff for zero
- `dosage_zerod_MAC_cutoff` - MAC cutoff for dosage zero
- `impute_method` - minor/mean/best_guess

**Important Notes**
1. Use `r_corr` (underscore) not `r.corr` in config file
2. Use colons `:` not semicolons `;` in `annotation_in_groupTest`
3. Quote comma-separated values: `CHROMOSOMES="1,2,3"`
4. Quote special characters: `maxMAF_in_groupTest="0.0001,0.001,0.01"`
5. Quote annotations: `annotation_in_groupTest="lof,missense:lof:missense:synonymous"`
6. Processing order: chr1 chunk1, chr1 chunk2, ..., chr2 chunk1, chr2 chunk2, ...
7. Merging order: chunk1, chunk2, ..., chunk10, chunk11 (numerical order)
8. Header only from chunk1, removed from all other chunks when merging

**Output Files**

*With Chunk Merging (MERGE_CHUNKS=yes)*
```
chr1_combined_results.txt    - All results for chr1 (header from chunk1 only)
chr2_combined_results.txt    - All results for chr2 (header from chunk1 only)
chr3_combined_results.txt    - All results for chr3 (header from chunk1 only)
...
chr22_combined_results.txt   - All results for chr22 (header from chunk1 only)
saige_run_summary.txt         - Job summary with status
saige_run.log                 - Detailed execution log
```

*Without Chunk Merging (MERGE_CHUNKS=no)*
```
chr1_chunk1_results.txt       - Results for chr1 chunk1
chr1_chunk2_results.txt       - Results for chr1 chunk2
chr1_chunk3_results.txt       - Results for chr1 chunk3
...
chr2_chunk1_results.txt       - Results for chr2 chunk1
chr2_chunk2_results.txt       - Results for chr2 chunk2
...
chr22_chunk1_results.txt      - Results for chr22 chunk1
chr22_chunk2_results.txt      - Results for chr22 chunk2
...
saige_run_summary.txt         - Job summary
saige_run.log                 - Detailed log
```

*Result File Columns*

- `Gene` - Gene identifier
- `Region` - Genomic region
- `Group` - Annotation group (e.g., lof, missense)
- `max_MAF` - Maximum MAF threshold used
- `Pvalue` - Association p-value
- `Pvalue_Burden` - Burden test p-value (if r_corr=0)
- `Pvalue_SKAT` - SKAT test p-value (if r_corr=0)
- `BETA_Burden` - Burden effect size
- `SE_Burden` - Standard error
- `N_VARIANTS` - Number of variants in gene
- `N_CASES` - Number of cases with variant
- `N_CONTROLS` - Number of controls with variant

*Log File Contents (saige_run.log)*

- Configuration summary
- SAIGE command for each job
- Processing status per chromosome/chunk
- Merge operations (if enabled)
- Error messages (if any)
- Timing information

*Summary File Contents (saige_run_summary.txt)*

- Header with date/time
- Table format: Chr | Chunk | Genotype_File | Output_File | Status
- Merge summary (if MERGE_CHUNKS=yes)
- Overall statistics (total/successful/failed jobs)

---

### Step 9: Post-Analysis

Interactive analysis tool for SAIGE-GENE test results with flexible operation modes.

**Usage:**

```bash
./step9_analyze_results.sh [operations]
```

**Arguments:**
- ```operations``` - Plus-separated list of analysis operations (optional, prompts if not provided)

**Available Operations:**

*Data Combination:*
- ```mergechrom``` - Merge chunked results by chromosome
- ```mergeall``` - Combine all chromosomes into single file

*Significance Filtering:*
- ```findsig``` - Extract significant results at multiple thresholds
- ```findgws``` - Extract genome-wide significant only (p < 5e-8)
- ```findsug``` - Extract suggestive results (p < 1e-5)
- ```findnom``` - Extract nominal results (p < 0.05)

*Gene Ranking:*
- ```top10``` - Extract top 10 genes
- ```top50``` - Extract top 50 genes
- ```top100``` - Extract top 100 genes

*Summary Statistics:*
- ```chromsum``` - Per-chromosome summary table
- ```fullsum``` - Complete summary report

*Visualization Data:*
- ```qqdata``` - Generate QQ plot data
- ```mandata``` - Generate Manhattan plot data
- ```plotdata``` - Generate both QQ and Manhattan data

*Presets:*
- ```standard``` - Standard analysis (mergeall+findsig+top50+fullsum)
- ```full``` - Complete analysis (all operations)
- ```quick``` - Quick overview (mergeall+top50+fullsum)

**Examples:**
```bash
# Interactive mode (shows menu)
./step9_analyze_results.sh

# Standard analysis
./step9_analyze_results.sh standard

# Custom combination
./step9_analyze_results.sh mergeall+findsig+top50+qqdata

# Quick check of top genes
./step9_analyze_results.sh quick

# Full comprehensive analysis
./step9_analyze_results.sh full

# Find significant genes only
./step9_analyze_results.sh findgws+top10
```

**Output:**
- ```all_results.txt``` - Combined results from all chromosomes
- ```genome_wide_sig.txt``` - Genome-wide significant genes (p < 5e-8)
- ```suggestive_sig.txt``` - Suggestive associations (p < 1e-5)
- ```nominal_sig.txt``` - Nominal significant genes (p < 0.05)
- ```top50_genes.txt``` - Top 50 genes ranked by p-value
- ```chromosome_summary.txt``` - Per-chromosome statistics table
- ```qq_plot_data.txt``` - Data for QQ plots (Expected vs Observed -log10P)
- ```manhattan_plot_data.txt``` - Data for Manhattan plots (Gene, CHR, P-value)
- ```analysis_summary.txt``` - Complete analysis summary report

**Other options**
```
# Manual result combination
head -1 saige_results/chr1_combined_results.txt > all_results.txt
tail -n +2 -q saige_results/chr*_combined_results.txt >> all_results.txt

# Extract significant results (p < 5e-8)
head -1 all_results.txt > genome_wide_sig.txt
awk 'NR>1 && $NF < 5e-8' all_results.txt >> genome_wide_sig.txt

# Extract suggestive results (p < 1e-5)
head -1 all_results.txt > suggestive_sig.txt
awk 'NR>1 && $NF < 1e-5' all_results.txt >> suggestive_sig.txt

# Count genes tested per chromosome
for chr in {1..22}; do
  count=$$(tail -n +2 saige_results/chr$${chr}_combined_results.txt | wc -l)
  echo "Chr$${chr}: $${count} genes"
done

# Top 20 genes
head -1 all_results.txt > top20.txt
tail -n +2 all_results.txt | sort -gk$(head -1 all_results.txt | tr '\t' '\n' | grep -n '^p' | cut -d: -f1) | head -20 >> top20.txt
```

---

## Important notes

```
# Gene Coordinate Files:
# - Pre-built coordinates available in gene_coords_ensembl115.tar.gz
# - No need to regenerate if using Ensembl 115 GRCh38
# - To use: tar -xzf gene_coords_ensembl115.tar.gz

# Step 5 Parameters:
# - Argument 4: buffer_kb (default: 10)
# - Argument 5: force_regen (yes/no, default: no)
# - Argument 6: merge_regen (yes/no, default: yes)

# Step 6 Parameters:
# - Argument 4: threads (default: 4)
# - Argument 5: output_format (bfile/pgen/vcf/bgen, default: bfile)
# - Argument 6: input_format (auto/bfile/pgen/vcf/bgen, default: auto)
# - Argument 7: input_prefix (default: chr)
# - Argument 8: chunk_size (0=no chunking, >0=genes per chunk, default: 0)

# Step 7 Parameters:
# - Argument 3: output_report (default: gene_bfiles/verification_report.txt)
# - Argument 4: output_format (bfile/pgen/vcf/bgen, default: bfile)

# Chunking Strategy:
# - No chunking (chunk_size=0): All genes in one file per chromosome
# - Small chunks (20-50): For memory-limited systems
# - Large chunks (100-200): For large datasets with parallel processing
```

---

## Quality Control

**After Step 1:**
```bash
# Check annotation distribution
awk -F'\t' 'NR>1 {print \$6}' chr1_anno.txt | sort | uniq -c
```

**After Step 2:**
```bash
# Verify var/anno pairing
awk '{
  if(NR%2==1) var_count=NF-2
  else {
    anno_count=NF-2
    if(var_count != anno_count) print "Mismatch at line " NR
  }
}' chr1_groups.txt

# Check for duplicates
awk '{if(NR%2==1) print \$1}' chr1_groups.txt | sort | uniq -d
```

**After Step 3:**
```bash
# Verify total gene count
wc -l all_genes_groups.txt  # Should be even

# Count unique genes
awk '{if(NR%2==1) print \$1}' all_genes_groups.txt | sort -u | wc -l
```

**After Step 5:**
```bash
# Check matched genes
wc -l plink_regions/matched_genes.txt

# Review missing genes
cat plink_regions/missing_genes.txt
```

**After Step 6:**
```bash
# Check extraction summary
column -t -s $'\t' gene_pfiles/extraction_summary.txt

# Verify variant counts
for chr in {1..22}; do
  echo "chr$${chr}: $$(wc -l < gene_pfiles/chr${chr}_genes.pvar) variants"
done
```

---

## Troubleshooting

**Windows line endings:**
```bash
sed -i 's/\r$//' *.sh
chmod +x *.sh
```

**LSF cluster submission:**
```bash
for chr in {1..22}; do
  bsub -q normal -n 4 -R "rusage[mem=8GB]" \
    "./step1_vep_ann_clean.sh chr${chr}.vcf.gz chr${chr}_anno.txt 4 vep_lof"
done
```

**Input file not found:**
```bash
# Use absolute paths
./step1_vep_ann_clean.sh /full/path/to/input.vcf.gz output.txt
```

**PLINK2 memory issues:**
```bash
# Increase memory limit
plink2 --pfile input --memory 16000 --threads 8 ...
```

**Missing genes in Step 5:**
```bash
# Check Ensembl version matches VEP version
# Review missing_genes.txt for outdated/retired gene symbols
# Consider using alternative gene names or ENSG IDs
```

**Format detection fails in Step 6:**
```bash
# Manually specify format
./step6_extract_genotypes_plink2.sh /data/geno plink_regions output vcf 8
```

**BGEN conversion for SAIGE:**
```bash
# Convert PGEN to BGEN (v1.2, 8-bit for SAIGE compatibility)
plink2 --pfile chr1_genes \
       --export bgen-1.2 bits=8 \
       --out chr1_genes_bgen

# Index BGEN file
bgenix -g chr1_genes_bgen.bgen -index
```

---

## Citation

- **VEP:** McLaren et al., Genome Biology (2016)
- **LOFTEE:** Karczewski et al., Nature (2020)
- **SAIGE:** Zhou et al., Nature Genetics (2018)
- **SAIGE-GENE+:** Zhou et al., Nature Genetics (2022)
- **Ensembl:** Cunningham et al., Nucleic Acids Research (2022)

---

## Support

- **Issues:** https://github.com/nicenzhou/vep-to-saige/issues
- **Email:** jyzhou@stanford.edu

---

*Version 0.0.1 | Last updated: 2026-04-23*

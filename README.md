# vep-to-saige: VEP Annotations to SAIGE Gene Groups Pipeline

Convert VEP-annotated VCF files into SAIGE-compatible gene group files for rare variant burden testing.

---

## Quick Start

```bash
# Clone and setup
git clone https://github.com/nicenzhou/vep-to-saige.git
cd vep-to-saige
sed -i 's/\r$//' *.sh  # Fix Windows line endings if needed
chmod +x *.sh

# Run pipeline

# Steps 1-3: VEP to SAIGE groups
./step1_vep_ann_clean.sh input.vcf.gz chr1_anno.txt 4 vep_lof
./step2_create_gene_groups.sh chr1_anno.txt chr1_groups.txt all keepall
./step3_merge_and_validate_groups.sh all_genes.txt .

# Steps 4-7: Gene coordinates and genotype extraction [Optional; for large WGS; Requires Internet for Step 4] 
./step4_download_gene_coords.sh 115 GRCh38 gene_coords
tar -xzf gene_coords_ensembl115.tar.gz
./step5_match_genes_to_groups.sh gene_coords all_genes_groups.txt plink_regions 10
./step6_extract_genotypes_plink2.sh /data/vcf plink_regions gene_pfiles vcf 16
./step7_verify_extraction.sh gene_pfiles plink_regions

# * Gene coords files are available in the gene_coords_ensembl115 folder on GitHub, where no regeneration is needed.
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

**Compatible with PLINK2a:**
```
# Extract genotypes using BED region file
plink2a --bfile chr1 \
  --extract bed1 plink_regions/chr1_regions.txt \
  --make-bed \
  --out chr1_genes
```

**Workflow Recommendations:**

1. **Initial Assessment**:
   ```
   # First run with default settings to check match rate
   ./step5_match_genes_to_groups.sh gene_coords groups.txt output 10
   # Review matching_summary.txt for match statistics
   ```

2. **High Match Rate (>95%)**:
   ```
   # Proceed with matched genes only
   ./step5_match_genes_to_groups.sh gene_coords groups.txt output 10
   ```

3. **Moderate Match Rate (80-95%)**:
   ```
   # Regenerate and merge for maximum recovery
   ./step5_match_genes_to_groups.sh gene_coords groups.txt output 10 yes yes
   ```

4. **Low Match Rate (<80%) or Quality Control**:
   ```
   # Regenerate but keep separate for review
   ./step5_match_genes_to_groups.sh gene_coords groups.txt output 10 yes no
   ```

5. **Custom Buffer Sizes**:
   - **Small buffer (10kb)**: Coding regions and immediate regulatory elements
   - **Medium buffer (50kb)**: Include proximal regulatory regions
   - **Large buffer (100kb)**: Include distal regulatory elements and TADs

**Quality Control Checks:**

```
# Check match statistics
cat plink_regions/matching_summary.txt

# View first few regions
head plink_regions/chr1_regions.txt

# Count genes per chromosome
for i in {1..22}; do 
  echo "Chr$i: $(wc -l < plink_regions/chr${i}_regions.txt) genes"
done

# Compare matched vs regenerated (if merge_regen=no)
echo "Matched: $(wc -l < plink_regions/matched_genes.txt)"
echo "Regenerated: $(wc -l < plink_regions/regenerated_genes.txt)"
echo "Missing: $(wc -l < plink_regions/missing_genes.txt)"
```

**Troubleshooting:**

| Issue | Solution |
|-------|----------|
| Many missing genes | Use ```force_regen=yes``` to recover from variant positions |
| Gene symbol mismatch | Check Ensembl release version matches your annotation |
| No variants for gene | Gene may not have qualifying variants in group file |
| Regenerated boundaries too wide | Reduce buffer_kb parameter |
| Overlapping regions | Normal; PLINK2a handles overlaps correctly |

---

### Step 6: Extract Genotypes with PLINK2

Extract genotypes for gene regions. Creates ONE BFILE/PFILE PER CHROMOSOME (efficient for SAIGE).

**Usage:**
```bash
./step6_extract_genotypes_plink2.sh <input_dir> <regions_dir> <output_dir> [format] [threads]
```

**Arguments:**
- ```input_dir``` - Directory with genotype files
- ```regions_dir``` - Directory with region files (from Step 5)
- ```output_dir``` - Output directory
- ```format``` - Input format: auto, vcf, pgen, bgen, bed (default: auto)
- ```threads``` - Number of threads (default: 4)

**Supported Input Formats:**
- **VCF**: ```chr1.vcf.gz```, ```chr2.vcf.gz```, ...
- **PGEN**: ```chr1.pgen/pvar/psam```, ```chr2.pgen/pvar/psam```, ...
- **BGEN**: ```chr1.bgen```, ```chr2.bgen```, ... (requires .sample file)
- **BED**: ```chr1.bed/bim/fam```, ```chr2.bed/bim/fam```, ...

**Examples:**
```bash
# Auto-detect format
./step6_extract_genotypes_plink2.sh /data/vcf plink_regions gene_pfiles

# Specify VCF input, 16 threads
./step6_extract_genotypes_plink2.sh /data/vcf plink_regions gene_pfiles vcf 16

# BGEN input
./step6_extract_genotypes_plink2.sh /data/bgen plink_regions gene_pfiles bgen 8

# Convert to BGEN for SAIGE
CONVERT_TO_BGEN=yes ./step6_extract_genotypes_plink2.sh /data/pgen plink_regions gene_bgen pgen 8
```

**Output:**
- ```chr*_genes.bed/bim/fam`` - Extracted genotypes per chromosome (BFILE format)
- ```extraction_summary.txt``` - Summary of extraction
- Optional: BGEN format for direct SAIGE use

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

## SAIGE Integration

```bash
# Run SAIGE-GENE+ burden test
Rscript step2_SPAtests.R        \
     --bgenFile=./input/genotype_100markers.bgen    \
     --bgenFileIndex=./input/genotype_100markers.bgen.bgi \
     --SAIGEOutputFile=./output/genotype_100markers_bgen_groupTest_out.txt \
     --chrom=1 \
     --LOCO=TRUE    \
     --AlleleOrder=ref-first \
     --minMAF=0 \
     --minMAC=0.5 \
     --sampleFile=./input/samplelist.txt \
     --GMMATmodelFile=./output/example_binary_fullGRM.rda \
     --varianceRatioFile=./output/example_binary_fullGRM.varianceRatio.txt	\
     --sparseGRMFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx   \
     --sparseGRMSampleIDFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt  \
     --groupFile=./input/all_genes_groups.txt	\
     --annotation_in_groupTest="lof,missense:lof,missense:lof:synonymous"        \
     --maxMAF_in_groupTest=0.0001,0.001,0.01	\
     --is_output_markerList_in_groupTest=TRUE \
     --is_fastTest=TRUE
```

*Full details see: https://saigegit.github.io/SAIGE-doc/docs/set_step2.html*
 
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

*Version 0.0.1 | Last updated: 2026-04-16*

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
./step1_vep_ann_clean.sh input.vcf.gz chr1_anno.txt 4 vep_lof
./step2_create_gene_groups.sh chr1_anno.txt chr1_groups.txt all keepall
./step3_merge_and_validate_groups.sh all_genes.txt .
```

---

## Requirements

- ```bcftools``` (VCF processing)
- ```awk``` (GNU awk recommended)
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
  ./step1_vep_ann_clean.sh chr${chr}.vcf.gz chr${chr}_anno.txt 4 vep_lof
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

| Mode | Behavior |
|------|----------|
| ```lof,missense,synonymous``` | Keep highest priority (default) |
| ```missense,lof,synonymous``` | Prioritize missense |
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
  ./step2_create_gene_groups.sh chr${chr}_anno.txt chr${chr}_groups.txt all keepall
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

Merge chromosome files into genome-wide file.

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

## Complete Workflow

```bash
# Step 1: Extract annotations (all chromosomes)
for chr in {1..22}; do
  ./step1_vep_ann_clean.sh \
    /path/to/chr${chr}.vcf.gz \
    chr${chr}_anno.txt \
    4 \
    vep_lof
done

# Step 2: Create gene groups
for chr in {1..22}; do
  ./step2_create_gene_groups.sh \
    chr${chr}_anno.txt \
    chr${chr}_groups.txt \
    all \
    keepall
done

# Step 3: Merge
./step3_merge_and_validate_groups.sh genome_wide_groups.txt .
```

---

## Quality Control

**After Step 1:**
```bash
# Check annotation distribution
awk -F'\t' 'NR>1 {print $6}' chr1_anno.txt | sort | uniq -c
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
awk '{if(NR%2==1) print $1}' chr1_groups.txt | sort | uniq -d
```

**After Step 3:**
```bash
# Verify total gene count
wc -l all_genes_groups.txt  # Should be even

# Count unique genes
awk '{if(NR%2==1) print $1}' all_genes_groups.txt | sort -u | wc -l
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
*Full details see: https://saigegit.github.io/SAIGE-doc/docs/set_step2.html
 
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

---

## Citation

- **VEP:** McLaren et al., Genome Biology (2016)
- **LOFTEE:** Karczewski et al., Nature (2020)
- **SAIGE:** Zhou et al., Nature Genetics (2018)
- **SAIGE-GENE+:** Zhou et al., Nature Genetics (2022)

---

## Support

- **Issues:** https://github.com/nicenzhou/vep-to-saige/issues
- **Email:** jyzhou@stanford.edu

---

*Version 0.0.1 | Last updated: 2026-04-16*

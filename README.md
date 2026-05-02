# vep-to-saige: VEP Annotations to SAIGE Gene-Based Association Pipeline

Convert VEP-annotated VCF files into SAIGE-compatible gene group files and run genome-wide set-based association analysis for rare variants.

## Overview

This pipeline processes VEP-annotated VCF files to create gene-based variant groups and performs rare variant association testing using SAIGE-GENE. It supports:

- **VEP annotation parsing** - Extract functional annotations from VEP-annotated VCFs
- **Gene group generation** - Create variant sets grouped by genes and functional categories
- **Genotype extraction** - Extract gene-specific genotypes in multiple formats (BGEN, PLINK, PGEN, VCF)
- **SAIGE-GENE testing** - Run genome-wide gene-based burden and SKAT-O tests
- **Result analysis** - Identify significant gene associations

---

## Clone and setup

```bash
git clone https://github.com/nicenzhou/vep-to-saige.git
cd vep-to-saige
# Fix Windows line endings if needed (GNU sed: sed -i 's/\r$//' codes/*.sh)
sed -i '' 's/\r$//' codes/*.sh
chmod +x codes/*.sh
```

---

## Quick Start: Full Pipeline

```
# Steps 1-3: VEP to SAIGE groups
./codes/step1_vep_ann_clean.sh input.vcf.gz chr1_anno.txt 4 vep_lof
./codes/step2_create_gene_groups.sh chr1_anno.txt chr1_groups.txt all keepall
./codes/step3_merge_and_validate_groups.sh all_genes.txt .

# Steps 4-5: Gene coordinates and region matching
./codes/step4_download_gene_coords.sh 115 GRCh38 gene_coords   # needs wget + network

# Or unpack pre-built coordinates
tar -xzf gene_coords_ensembl115.tar.gz

./codes/step5_match_genes_to_groups.sh gene_coords all_genes_groups.txt plink_regions 10 no yes

# Step 6: Extract genotypes (chunk last arg: 0 = whole chr)
./codes/step6_extract_genotypes_plink2a.sh /data/genotypes plink_regions gene_bfiles 16 bgen auto chr 0

./codes/step7_verify_extraction.sh gene_bfiles plink_regions verification_report.txt bgen

# Step 8: SAIGE-GENE — configure, validate, run
./codes/step8_pre1_build_saige_config.sh

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
annotation_in_groupTest=lof,missense:lof:missense:synonymous
is_Firth_beta=TRUE
r_corr=0
EOF

./codes/step8_pre2_validate_saige_config.sh saige_config.txt
./codes/step8_run_saige_gene_tests.sh saige_config.txt

./codes/step9_analyze_results.sh --results-dir saige_results standard
```

---

## Quick Start: Alternative Workflows

```
#--------------------------------------------------------------------------
# Workflow 1: Basic gene-based analysis (no genotype extraction)
#--------------------------------------------------------------------------

# Steps 1-3: Create gene groups from VEP
./codes/step1_vep_ann_clean.sh input.vcf.gz chr1_anno.txt 4 vep_lof
./codes/step2_create_gene_groups.sh chr1_anno.txt chr1_groups.txt all keepall
./codes/step3_merge_and_validate_groups.sh all_genes.txt .

# Use groups directly with your genotypes in SAIGE
# (Skip steps 4-7 if you already have organized genotype files)

#--------------------------------------------------------------------------
# Workflow 2: Large WGS dataset with chunking
#--------------------------------------------------------------------------

# Steps 1-3: Create gene groups
./codes/step1_vep_ann_clean.sh wgs.vcf.gz wgs_anno.txt 16 vep_lof
./codes/step2_create_gene_groups.sh wgs_anno.txt wgs_groups.txt all keepall
./codes/step3_merge_and_validate_groups.sh all_genes.txt .

# Steps 4-5: Gene coordinates with regeneration for missing genes
tar -xzf gene_coords_ensembl115.tar.gz
./codes/step5_match_genes_to_groups.sh gene_coords all_genes_groups.txt plink_regions 50 yes yes

# Step 6: Extract with chunking (20 genes per chunk)
./codes/step6_extract_genotypes_plink2a.sh /data/bgen plink_regions gene_bfiles 16 bgen bgen imputed_chr 20

# Step 7: Verify chunked extraction
./codes/step7_verify_extraction.sh gene_bfiles plink_regions qc_report.txt bgen

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

./codes/step8_run_saige_gene_tests.sh wgs_config.txt
./codes/step9_analyze_results.sh --results-dir wgs_saige_results standard

#--------------------------------------------------------------------------
# Workflow 3: Quick test on chromosome 22
#--------------------------------------------------------------------------

# Steps 1-3: Process chr22 only
./codes/step1_vep_ann_clean.sh chr22.vcf.gz chr22_anno.txt 4 vep_lof
./codes/step2_create_gene_groups.sh chr22_anno.txt chr22_groups.txt all keepall

# Steps 4-5: Match genes
tar -xzf gene_coords_ensembl115.tar.gz
./codes/step5_match_genes_to_groups.sh gene_coords chr22_groups.txt test_regions 10

# Step 6: Extract genotypes
./codes/step6_extract_genotypes_plink2a.sh /data/bfiles test_regions test_output 4 bfile auto chr 0

# Step 7: Verify
./codes/step7_verify_extraction.sh test_output test_regions test_qc.txt bfile

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
r_corr=0
EOF

./codes/step8_run_saige_gene_tests.sh test_config.txt
./codes/step9_analyze_results.sh --results-dir test_results quick

#--------------------------------------------------------------------------
# Workflow 4: Using interactive configuration builder
#--------------------------------------------------------------------------

# Steps 1-5: Same as above
# ...

# Step 6: Extract genotypes
./codes/step6_extract_genotypes_plink2a.sh /data/bgen plink_regions gene_bfiles 16 bgen bgen chr 0

# Step 7: Verify
./codes/step7_verify_extraction.sh gene_bfiles plink_regions

./codes/step8_pre1_build_saige_config.sh
./codes/step8_pre2_validate_saige_config.sh saige_config.txt
./codes/step8_run_saige_gene_tests.sh saige_config.txt
./codes/step9_analyze_results.sh --results-dir saige_results standard
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
# └── codes/*.sh                 # Pipeline scripts (incl. step9_analyze_results.sh, step9_noninteractive_example.sh)
```

---

## Requirements

- ```bcftools``` — Step 1
- ```awk``` — pipeline (GNU awk recommended)
- ```wget``` — Step 4 download
- ```plink2a``` — Step 6 extraction
- Bash 4.0+
- Optional: ```python3``` / ```certifi``` (Step 9 Ensembl), ```Rscript``` + ```ggplot2``` (Step 9 plots; ```ggrepel``` optional, improves Manhattan text labels)

---


## Pipeline scripts (steps 1–9)

Scripts are under **`codes/`**. Use **`./codes/<script>.sh`**, **`cd codes`**, or **`export PATH="$PWD/codes:$PATH"`**.

| Step | Script | Role |
|------|--------|------|
| 1 | `step1_vep_ann_clean.sh` | VEP → annotation TSV |
| 2 | `step2_create_gene_groups.sh` | TSV → SAIGE `GENE var …` / `GENE anno …` |
| 3 | `step3_merge_and_validate_groups.sh` | Merge per-chromosome group files |
| 4 | `step4_download_gene_coords.sh` | Ensembl GTF → `gene_coords/chr*_genes.txt` |
| 5 | `step5_match_genes_to_groups.sh` | Coordinates + groups → PLINK regions |
| 6 | `step6_extract_genotypes_plink2a.sh` | Extract cohort genotypes per gene region |
| 7 | `step7_verify_extraction.sh` | QC vs region lists |
| 8 | `step8_pre1_build_saige_config.sh`, `step8_pre2_validate_saige_config.sh`, `step8_run_saige_gene_tests.sh` | SAIGE-GENE config + Step 2 |
| 9 | `step9_analyze_results.sh` (optional `step9_noninteractive_example.sh`) | Merge & summarize SAIGE-GENE results, optional QQ/Manhattan |

**Performance (recent versions):** Step **4** builds per-chromosome tables in **one pass** over `all_genes_coords.txt` (not one full scan per chromosome). Step **5** does **one `awk` per chromosome** for filter + BED + gene list (no large temp “matched” files). Step **1** stderr summaries use **awk counts + one numeric sort** per table instead of `sort | uniq -c` on every row. Step **3** “unique genes” uses **var lines only** (`NR%2==1`). Step **6** BGEN log parsing avoids GNU-only `grep -P` (portable on macOS).

### Step 1 — `step1_vep_ann_clean.sh`

```bash
./codes/step1_vep_ann_clean.sh <input.vcf.gz> <output.txt> [threads] [lof_mode]
```

| | Default | Notes |
|---|---------|-------|
| `threads` | `1` | Passed to `bcftools view --threads` |
| `lof_mode` | `vep_lof` | `loftee_only` (LOFTEE HC only), `vep_lof` (VEP LoF consequences), `any` (HC or VEP LoF) |

Requires **`bcftools`**. Drops `LoF=LC`, deduplicates identical annotation rows. **Output:** tab-separated header + `Variant`, `Gene`, `Consequence`, `LoF`, `LoF_flags`, `Group`, `Total_Anno`, `NON_CAN_SPLICE_Count`.

### Step 2 — `step2_create_gene_groups.sh`

```bash
./codes/step2_create_gene_groups.sh <anno.txt> <groups.txt> [filter] [priority]
```

| | Default | Notes |
|---|---------|-------|
| `filter` | `all` | `lof`, `missense`, `synonymous`, `lof+missense` |
| `priority` | `lof,missense,synonymous` | Resolves duplicate variant rows; use **`keepall`** to keep every annotation |

Conflict rows (same priority) → **`<groups.txt>.conflict`**. Filters missing symbols, bad splice counts, and “all NON_CAN_SPLICE” LoF rows per script rules.

### Step 3 — `step3_merge_and_validate_groups.sh`

```bash
./codes/step3_merge_and_validate_groups.sh <merged_out.txt> [input_dir] [pattern]
```

Merges **`chr*_groups.txt`** in genomic order. Auto-detects **`chr1`** vs **`chr01`** naming, or supply a glob **`pattern`** (quoted), e.g. `"custom_chr*_groups.txt"`. Concatenation only (no re-validation).

### Step 4 — `step4_download_gene_coords.sh`

```bash
./codes/step4_download_gene_coords.sh [ensembl_release] [GRCh38|GRCh37] [out_dir]
```

Defaults: **`115`**, **`GRCh38`**, **`gene_coords`**. Needs network + **`wget`**. Writes `all_genes_coords.txt`, then **splits by chromosome in a single read** of that file, removes downloaded GTF, then **`gene_coords_ensembl<release>.tar.gz`**. **GRCh37:** GTF file is **`Homo_sapiens.GRCh37.<release>.gtf.gz`** matching **`ensembl_release`**.

### Step 5 — `step5_match_genes_to_groups.sh`

```bash
./codes/step5_match_genes_to_groups.sh <gene_coords_dir> <merged_groups.txt> <out_dir> [buffer_kb] [force_regen] [merge_regen]
```

| | Default | Notes |
|---|---------|-------|
| `buffer_kb` | `10` | Added around gene body |
| `force_regen` | `no` | If `yes`, infer intervals from variant positions when symbol missing from coords |
| `merge_regen` | `yes` | Merge recovered genes into main **`chr*_regions.txt`** vs separate **`*_recovered.txt`** |

Output: BED-like **`chr*_regions.txt`**, **`matched_genes.txt`**, **`missing_genes.txt`**, **`matching_summary.txt`**, log under **`out_dir`**. Matching and BED construction share **one pass** over each coordinate file (no per-chr temp match file).

### Step 6 — `step6_extract_genotypes_plink2a.sh`

```bash
./codes/step6_extract_genotypes_plink2a.sh <genotype_dir> <regions_dir> <out_dir> \
  [threads] [output_fmt] [input_fmt] [chr_prefix] [genes_per_chunk]
```

Requires **`plink2a`** (optional **`module`** load). **`output_fmt` / `input_fmt`:** `bfile`, `pgen`, `vcf`, `bgen`; **`auto`** detects inputs. **`chr_prefix`:** e.g. `chr`, `imputed_chr`, or `""` for `1.bed`… **`genes_per_chunk`:** `0` = whole chromosome; `20` → **`chr1_genes_chunk1`** … **`extraction_summary.txt`** + **`extraction.log`**.

### Step 7 — `step7_verify_extraction.sh`

```bash
./codes/step7_verify_extraction.sh <extract_dir> <regions_dir> [report.txt] [format]
```

**`format`** must match step 6 (`bfile`, `pgen`, `vcf`, `bgen`). Report: counts per chr/chunk, sizes, optional MAF histogram if **`plink2a`** runs `--freq`.

### Step 8 — SAIGE-GENE (configuration + run)

| Script | Use |
|--------|-----|
| `codes/step8_pre1_build_saige_config.sh` | Interactive **`key=value`** file |
| `codes/step8_pre2_validate_saige_config.sh <file>` | Path and option checks (`source`s the file) |
| `codes/step8_run_saige_gene_tests.sh <file>` | SPAtests / Step 2 jobs, chunk merge optional |

**Must set:** `GENOTYPE_DIR`, `OUTPUT_DIR`, `GMMAT_MODEL`, `VARIANCE_RATIO`, and either **`GROUP_FILE`** with **`GROUP_FILE_BY_CHR=no`** or **`GROUP_FILE_BY_CHR=yes`** + **`GROUP_DIR`**. Set **`INPUT_FORMAT`**, **`THREADS`**, quoted **`CHROMOSOMES`**, **`CHUNKED_INPUT`** / **`MERGE_CHUNKS`** when using chunks. In-file SAIGE options: use **`r_corr`** (maps to `--r.corr`); quote **`maxMAF_in_groupTest`**; **`annotation_in_groupTest`** uses **colons** between burden masks.

---

### Step 9: Post-Analysis

**Script:** `codes/step9_analyze_results.sh`. Reads SAIGE-GENE chromosome outputs (`chr*_combined_results.txt`, `chr*_chunk*_results.txt`, or `chr*_all_results.txt`). **Default results directory is `.`** — override with **`--results-dir`**, **`-d`**, **`--dir`**, or **`STEP9_RESULTS_DIR`** so the script need not sit next to outputs.

**Wrapper (optional):** `codes/step9_noninteractive_example.sh` — edit the CONFIG block or pass env vars (`RESULTS_DIR`, `STEP9_SCRIPT`, `OPERATIONS`, `PRESET`, etc.); runs `step9_analyze_results.sh` (same repo folder by default). Example:

```bash
chmod +x codes/step9_noninteractive_example.sh
RESULTS_DIR=/path/to/saige_out PRESET=quick ./codes/step9_noninteractive_example.sh
# Equivalent: STEP9_RESULTS_DIR=/path/to/saige_out PRESET=quick ./codes/step9_noninteractive_example.sh
```

**Optional dependencies:** `python3` (Ensembl REST lookups), `Rscript` (PNG QQ/Manhattan; **`ggplot2`** / **`ggrepel`** optional), **`certifi`** on macOS if Ensembl TLS fails.

**Non-interactive** (flags first, then positionals):

```bash
./codes/step9_analyze_results.sh [--results-dir|-d|--dir PATH] OPERATIONS \
  [gene_coords_file] [position_mode] [filter_group] [filter_maf] \
  [coord_source] [ensembl_build] [ensembl_release]
```

| Flag / positional | Meaning |
|-------------------|---------|
| `--results-dir` etc. | Results folder (same as **`STEP9_RESULTS_DIR`**) |
| `OPERATIONS` | One preset (**`standard`**, **`full`**, **`quick`**) or **`+`‑separated** operations (e.g. `mergeall+findsig+top50`) |
| `gene_coords_file` | Optional coordinate TSV for Manhattan (chrom, start, end, gene, …) |
| `position_mode` | Gene position on chromosome: **`start`**, **`end`**, **`midpoint`** |
| `filter_group` | Optional; restrict derived tables/plots to annotation group(s) (e.g. `lof` or `lof;missense`) |
| `filter_maf` | Optional; e.g. **`0.01`** keeps tests with max_MAF ≤ 0.01 |
| `coord_source` | **`sequence`** (default) or **`ensembl`** when no gene coordinate file is given |
| `ensembl_build` | **`37`** or **`38`** when **`coord_source`** is **`ensembl`** |
| `ensembl_release` | Optional Ensembl release (e.g. **`115`**) or custom REST base URL |

**Plot / Cauchy env:** **`STEP9_PLOT_UNFILTERED`** (default **`1`**: one QQ + one Manhattan for all groups × all max_MAF), **`STEP9_MANHATTAN_LABEL_TOP_N`**, **`STEP9_MANHATTAN_LABEL_EXTRA`** (comma-separated symbols always labeled in addition to the top‑N by *P*; case-insensitive; symbols absent from the Manhattan data are skipped), **`STEP9_MANHATTAN_LABEL_EXTRA_FILE`** (tab-separated file merged into the same symbol list; if the first row names a column **`Gene`**, **`Symbol`**, **`gene_symbol`**, or **`HGNC`**, that column is used; otherwise every row including the first is treated as data and column 1 is used). **`STEP9_BONFERRONI_MAF_TESTS`** (Bonferroni divisor for the green line on **`manhattan_plot.png`**, default **`3`**). **`STEP9_CAUCHY_MODE`** = **`off`**, **`plots`**, **`pipeline`**, or **`full`** — progressively omit the Cauchy annotation group from plots, from significance/top/summary tables, or from full-summary totals (**`all_results.txt`** is never edited). Legacy **`STEP9_EXCLUDE_CAUCHY=1`** behaves like **`plots`** if **`STEP9_CAUCHY_MODE`** is unset. **`STEP9_ENSEMBL_SSL_VERIFY`** (`0` = insecure fallback). **Ensembl REST tuning:** **`STEP9_ENSEMBL_BATCH_SIZE`**, **`STEP9_ENSEMBL_PARALLEL`**, **`STEP9_ENSEMBL_POST_TIMEOUT`**.

**Interactive mode** also prompts for Manhattan label count and optional **extra gene labels** (comma list or path to a symbol file). Non-interactive runs use the env vars above only.

**Presets (actual chain in script):** **`standard`** — detect_files, mergeall, listgroups, findsig, top50, groupsum, makeplots, fullsum. **`quick`** — detect_files, mergeall, listgroups, top50, makeplots, fullsum (no findsig or groupsum). **`full`** — detect_files, mergechrom, mergeall, listgroups, findsig, top10, top50, top100, chromsum, groupsum, makeplots, fullsum.

**Common ops:** `mergechrom`, `mergeall`, `listgroups`; `findsig`, `findgws`, `findsug`, `findnom`; `top10`, `top50`, `top100`; `chromsum`, `groupsum`, `fullsum`; `qqdata`, `mandata`; **`plotdata`** and **`makeplots`** are **equivalent** (both call the same code path: refresh `qq_plot_data.txt` / `manhattan_plot_data.txt` and write PNGs and group statistics when R is available).

**Artifacts:** `all_results.txt`, significance-tier files, `top*_genes.txt`, **`qq_plot_data.txt`**, **`manhattan_plot_data.txt`**, **`qq_plot.png`**, **`manhattan_plot.png`** (Bonferroni −log10 *P*), **`group_statistics.txt`** (from makeplots when R runs), `analysis_summary.txt`; helpers `.all_results_no_cauchy.txt` (when Cauchy is omitted from tables), `.gene_coords_lookup.txt`, **`.step9_generate_plots.R`** (written under the results directory when plot generation runs).

```bash
# Interactive (menu + prompts)
./codes/step9_analyze_results.sh

./codes/step9_analyze_results.sh --results-dir /path/to/saige_out standard
STEP9_RESULTS_DIR=/path/to/out ./codes/step9_analyze_results.sh quick
```

---

## Important notes

- Pre-built **gene coordinates** (e.g. `gene_coords_ensembl115.tar.gz`) can replace Step 4 when the Ensembl release matches your VEP build.
- Step **5** optional args: **`buffer_kb`**, **`force_regen`**, **`merge_regen`** — details under Step 5 above.
- Step **6** last argument **`genes_per_chunk`**: **`0`** = one output per chromosome; **`20`**–**`50`** typical when RAM is tight.
- SAIGE config: use **`r_corr`** in the file (not **`r.corr`**); **`annotation_in_groupTest`** uses **colons** between masks.

## Quality Control

**After Step 1:**
```bash
awk -F'\t' 'NR>1 {print $6}' chr1_anno.txt | sort | uniq -c
```

**After Step 2:**
```bash
awk '{
  if(NR%2==1) var_count=NF-2
  else {
    anno_count=NF-2
    if(var_count != anno_count) print "Mismatch at line " NR
  }
}' chr1_groups.txt

awk '{if(NR%2==1) print $1}' chr1_groups.txt | sort | uniq -d
```

**After Step 3:**
```bash
wc -l all_genes_groups.txt

awk '{if(NR%2==1) print $1}' all_genes_groups.txt | sort -u | wc -l
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
  echo "chr${chr}: $(wc -l < gene_pfiles/chr${chr}_genes.pvar) variants"
done
```

---

## Troubleshooting

**Windows line endings:**
```bash
sed -i '' 's/\r$//' codes/*.sh   # GNU sed: sed -i 's/\r$//' codes/*.sh
chmod +x codes/*.sh
```

**LSF cluster submission:**
```bash
for chr in {1..22}; do
  bsub -q normal -n 4 -R "rusage[mem=8GB]" \
    "./codes/step1_vep_ann_clean.sh chr${chr}.vcf.gz chr${chr}_anno.txt 4 vep_lof"
done
```

**Input file not found:**
```bash
# Use absolute paths
./codes/step1_vep_ann_clean.sh /full/path/to/input.vcf.gz output.txt
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
./codes/step6_extract_genotypes_plink2a.sh /data/geno plink_regions out 8 vcf vcf chr 0
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

*Version 2.0.0 | Last updated: 2026-05-01*

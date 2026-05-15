# vep-to-saige: VEP Annotations to SAIGE Gene-Based Association Pipeline

Convert VEP-annotated VCF files into SAIGE-compatible gene group files and run genome-wide set-based association analysis for rare variants.

## Contents

**Introduction**

- [Overview](#overview)

**Getting started**

- [Clone and setup](#clone-and-setup)
- [Quick Start: Full Pipeline](#quick-start-full-pipeline)
- [Quick Start: Alternative Workflows](#quick-start-alternative-workflows)

**Reference**

- [Expected output structure](#expected-output-structure)
- [Requirements](#requirements)
- [Pipeline scripts (steps 1-11)](#pipeline-scripts-steps-1-11) — Step **9** = one results folder; Step **10** = multiple folders (documented after Step 9); Step **11** = SAIGE groups → REGENIE anno/set-list (optional REGENIE burden workflow)

**Notes and troubleshooting**

- [Important notes](#important-notes)
- [Troubleshooting](#troubleshooting)

**Citation and support**

- [Citation](#citation)
- [Support](#support)

## Overview

| Stage | What it does |
|-------|----------------|
| **VEP → groups** | Parse VEP VCFs → SAIGE **`GENE var …` / `GENE anno …`** groups |
| **Coordinates & regions** | Ensembl / matched coords → PLINK-style regions |
| **Genotypes** | Extract per-gene regions (**BGEN**, **PLINK**, **PGEN**, **VCF**) |
| **SAIGE-GENE** | Burden + SKAT-O tests genome-wide |
| **Summarize & plots** | **Step 9** = one output dir · **Step 10** = batch many dirs → **`SummaryResults/`** |
| **Optional: SAIGE → REGENIE** | **Step 11** (optional) — SAIGE **`chr*_group.txt`** → per-chr **`--anno-file`** / **`--set-list`** + **`mask_def.example.txt`** (**`step11_saige_group_to_regenie.sh`**) |

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

```bash
# Steps 1-3: VEP to SAIGE groups
./codes/step1_vep_ann_clean.sh input.vcf.gz chr1_anno.txt 4 vep_lof
./codes/step2_create_gene_groups.sh chr1_anno.txt chr1_groups.txt all keepall
./codes/step3_merge_and_validate_groups.sh all_genes_groups.txt .

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

# Step **9** = one directory of **`chr*_combined_results.txt`**. Step **10** = batch across folders.
./codes/step9_analyze_results.sh --results-dir saige_results standard

# -----------------------------------------------------------------------------
# Step 10 — only if you have MULTIPLE sibling results folders (not one saige_results/).
# One folder of results → Step 9 above is enough; skip Step 10.
# -----------------------------------------------------------------------------
# chmod +x codes/step10_run_step9_all_folders.sh
# ./codes/step10_run_step9_all_folders.sh                     # standard per folder → SummaryResults/
# ./codes/step10_run_step9_all_folders.sh -i AFR EUR AMR       # interactive; Y = same params for all folders
```

*Further detail: [Pipeline scripts (steps 1-11)](#pipeline-scripts-steps-1-11) → [Step 10](#step-10-batch-step-9-multiple-cohort-folders).*

---

## Quick Start: Alternative Workflows

```bash
#--------------------------------------------------------------------------
# Workflow 1: Basic gene-based analysis (no genotype extraction)
#--------------------------------------------------------------------------

# Steps 1-3: Create gene groups from VEP
./codes/step1_vep_ann_clean.sh input.vcf.gz chr1_anno.txt 4 vep_lof
./codes/step2_create_gene_groups.sh chr1_anno.txt chr1_groups.txt all keepall
./codes/step3_merge_and_validate_groups.sh all_genes_groups.txt .

# Use groups directly with your genotypes in SAIGE
# (Skip steps 4-7 if you already have organized genotype files)

#--------------------------------------------------------------------------
# Workflow 2: Large WGS dataset with chunking
#--------------------------------------------------------------------------

# Steps 1-3: Create gene groups
./codes/step1_vep_ann_clean.sh wgs.vcf.gz wgs_anno.txt 16 vep_lof
./codes/step2_create_gene_groups.sh wgs_anno.txt wgs_groups.txt all keepall
./codes/step3_merge_and_validate_groups.sh all_genes_groups.txt .

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

#--------------------------------------------------------------------------
# Optional Step 10 — only if you have MULTIPLE sibling SAIGE results folders
# (e.g. AFR/, EUR/). One folder → Step 9 above is enough; skip Step 10.
#--------------------------------------------------------------------------
# chmod +x codes/step10_run_step9_all_folders.sh
# ./codes/step10_run_step9_all_folders.sh
# ./codes/step10_run_step9_all_folders.sh -i AFR EUR AMR
```

*Further detail: [Pipeline scripts (steps 1-11)](#pipeline-scripts-steps-1-11) → [Step 10](#step-10-batch-step-9-multiple-cohort-folders).*

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
# ├── SummaryResults/            # Step 10 only (multiple cohort folders): prefixed PNGs + tables per cohort
# ├── regenie_masks_from_saige/ # Step 11 optional: chrNN/chrNN.regenie.annotations.txt, chrNN/chrNN.regenie.setlist.txt, mask_def.example.txt
# └── codes/*.sh                 # Pipeline scripts (incl. step9_analyze_results.sh, step9_noninteractive_example.sh, step10_run_step9_all_folders.sh)
```

---

## Requirements

| Tool | Used in |
|------|---------|
| **`bcftools`** | Step 1 |
| **`awk`** (GNU awk recommended) | Many steps |
| **`wget`** | Step 4 |
| **`plink2a`** | Steps 6–7 |
| **Bash 4.0+** | All scripts |

| Optional | Used for |
|----------|----------|
| **`python3`**, **`certifi`** | Step 9 Ensembl lookups |
| **`Rscript`**, **`ggplot2`**, **`ggrepel`** | Step 9 QQ / Manhattan PNGs |

---

## Pipeline scripts (steps 1-11)

**How to run scripts**

| Your layout | Command pattern |
|-------------|-----------------|
| Scripts under **`codes/`** | `./codes/stepN_*.sh …` |
| Scripts next to cohort folders (some clones) | `./step9_*.sh` / `./step10_*.sh` from repo root |

**Which post-SAIGE step?**

| Situation | Section |
|-----------|---------|
| One SAIGE output tree | [Step 9](#step-9-post-analysis) |
| Multiple sibling result folders | [Step 10](#step-10-batch-step-9-multiple-cohort-folders) |

| Step | Script | Role |
|------|--------|------|
| 1 | `step1_vep_ann_clean.sh` | VEP → annotation TSV |
| 2 | `step2_create_gene_groups.sh` | TSV → SAIGE `GENE var …` / `GENE anno …` |
| 3 | `step3_merge_and_validate_groups.sh` | Merge per-chromosome group files |
| 4 | `step4_download_gene_coords.sh` | Ensembl GTF → `gene_coords/chr*_genes.txt` |
| 5 | `step5_match_genes_to_groups.sh` | Coordinates + groups → PLINK regions |
| 6 | `step6_extract_genotypes_plink2a.sh` | Extract cohort genotypes per gene region |
| 7 | `step7_verify_extraction.sh` | QC vs region lists |
| 8 | `step8_pre1_build_saige_config.sh`, `step8_pre2_validate_saige_config.sh`, `step8_run_saige_gene_tests.sh` | SAIGE-GENE config + run test jobs |
| 9 | `step9_analyze_results.sh` (optional `step9_noninteractive_example.sh`) | **Single results folder** — merge & summarize SAIGE-GENE outputs; QQ/Manhattan (**`STEP9_MANHATTAN_P_MODE`**) |
| 10 | `step10_run_step9_all_folders.sh` | **Multiple results folders** — run Step 9 once per folder; **`--help`**; prefixed copies → **`SummaryResults/`** ([Step 10](#step-10-batch-step-9-multiple-cohort-folders)) |
| 11 | `step11_saige_group_to_regenie.sh` | SAIGE **`chr*_group.txt`** → REGENIE **`--anno-file`** / **`--set-list`** per chromosome + example **`--mask-def`** (**Step 11** below) |

**Performance (recent versions)**

| Step | Note |
|------|------|
| **4** | Streams **`all_genes_coords.txt`** once per build |
| **5** | One **`awk`** per chr for filter + BED + gene list |
| **1**, **3**, **6** | Lighter parsing (portable **awk**, no GNU-only **`grep -P`**) |
| **9** `fullsum` | Single pass counts; may reuse **`annotation_group_summary`** from **`groupsum`** |
| **10** | Reuses coord cache **`.step9_batch_coord_cache/<hash>/`** for cohort 2+ when Step 9 args match |

**Reading the step tables:** Positionals are **in order** (`<required>` `[optional]`). **Expected input** = file or value type. **Options** = allowed keywords. **Defaults** apply when you skip optional args.

### Step 1 — `step1_vep_ann_clean.sh`

```bash
./codes/step1_vep_ann_clean.sh <input.vcf.gz> <output.txt> [threads] [lof_mode]
```

| Argument | Expected input | Options / default |
|----------|----------------|-------------------|
| `<input.vcf.gz>` | Path to input VCF (**bgzip** `.vcf.gz`) with VEP annotations in `INFO` | Must be readable by **`bcftools`** |
| `<output.txt>` | Path for the annotation table (created/overwritten) | Tab-separated text |
| `[threads]` | Integer thread count for `bcftools view --threads` | Default **`1`** |
| `[lof_mode]` | Which LoF-like variants to keep | **`loftee_only`** — LOFTEE high-confidence only; **`vep_lof`** (default) — VEP LoF consequences; **`any`** — HC or VEP LoF |

Requires **`bcftools`**. Drops `LoF=LC`, deduplicates identical annotation rows. **Output columns:** `Variant`, `Gene`, `Consequence`, `LoF`, `LoF_flags`, `Group`, `Total_Anno`, `NON_CAN_SPLICE_Count`.

### Step 2 — `step2_create_gene_groups.sh`

```bash
./codes/step2_create_gene_groups.sh <anno.txt> <groups.txt> [filter] [priority]
```

| Argument | Expected input | Options / default |
|----------|----------------|-------------------|
| `<anno.txt>` | Step 1–style annotation TSV | Same columns as Step 1 output |
| `<groups.txt>` | Output path for SAIGE **`GENE var …` / `GENE anno …`** pairs | Variant line then annotation line per variant |
| `[filter]` | Which consequence classes to include | **`all`** (default), **`lof`**, **`missense`**, **`synonymous`**, **`lof+missense`** |
| `[priority]` | Comma list used when one variant maps to multiple annotations | Default **`lof,missense,synonymous`**; use **`keepall`** to retain every annotation row |

Conflict rows (same priority) → **`<groups.txt>.conflict`**. Filters missing symbols, bad splice counts, and “all NON_CAN_SPLICE” LoF rows per script rules.

### Step 3 — `step3_merge_and_validate_groups.sh`

```bash
./codes/step3_merge_and_validate_groups.sh <merged_out.txt> [input_dir] [pattern]
```

| Argument | Expected input | Options / default |
|----------|----------------|-------------------|
| `<merged_out.txt>` | Output path for the merged groups file (e.g. `all_genes_groups.txt`) | Single combined file in genomic order |
| `[input_dir]` | Directory that holds per-chromosome **`chr*_groups.txt`** files | Default **`.`** (current directory) |
| `[pattern]` | Optional **quoted** glob to select input files | Default auto-detects **`chr1`** vs **`chr01`** naming; e.g. `"custom_chr*_groups.txt"` |

Merges inputs in genomic order. Concatenation only (no re-validation).

### Step 4 — `step4_download_gene_coords.sh`

```bash
./codes/step4_download_gene_coords.sh [ensembl_release] [GRCh38|GRCh37] [out_dir]
```

| Argument | Expected input | Options / default |
|----------|----------------|-------------------|
| `[ensembl_release]` | Ensembl **release number** (matches GTF on Ensembl FTP) | Default **`115`** |
| Genome build | Reference assembly for the GTF | **`GRCh38`** or **`GRCh37`** (positional arg 2); default **`GRCh38`** |
| `[out_dir]` | Directory name to create for `chr*_genes.txt` | Default **`gene_coords`** |

Needs network + **`wget`**. Writes `all_genes_coords.txt`, splits by chromosome, removes downloaded GTF, archives **`gene_coords_ensembl<release>.tar.gz`**. **GRCh37:** GTF **`Homo_sapiens.GRCh37.<release>.gtf.gz`** must match **`ensembl_release`**.

### Step 5 — `step5_match_genes_to_groups.sh`

```bash
./codes/step5_match_genes_to_groups.sh <gene_coords_dir> <merged_groups.txt> <out_dir> [buffer_kb] [force_regen] [merge_regen]
```

| Argument | Expected input | Options / default |
|----------|----------------|-------------------|
| `<gene_coords_dir>` | Directory from Step **4** (`gene_coords/`) containing **`chr*_genes.txt`** | Path to folder, not a single file |
| `<merged_groups.txt>` | Step **3** merged groups file (e.g. `all_genes_groups.txt`) | One variant/anno block per gene group |
| `<out_dir>` | Directory to write **`plink_regions/`**-style outputs | Created if needed |
| `[buffer_kb]` | Extra kb padded around each gene interval | Integer ≥ **`0`**; default **`10`** |
| `[force_regen]` | Recover coords from variant positions when symbol missing | **`yes`** / **`no`**; default **`no`** |
| `[merge_regen]` | Merge recovered intervals into main chr files | **`yes`** / **`no`**; default **`yes`** |

Output: BED-like **`chr*_regions.txt`**, **`matched_genes.txt`**, **`missing_genes.txt`**, **`matching_summary.txt`**, log under **`out_dir`**. Matching and BED construction share **one pass** over each coordinate file (no per-chr temp match file).

### Step 6 — `step6_extract_genotypes_plink2a.sh`

```bash
./codes/step6_extract_genotypes_plink2a.sh <genotype_dir> <regions_dir> <out_dir> \
  [threads] [output_fmt] [input_fmt] [chr_prefix] [genes_per_chunk]
```

| Argument | Expected input | Options / default |
|----------|----------------|-------------------|
| `<genotype_dir>` | Directory with cohort genotype files (**`.bed/.bim/.fam`**, **`.pgen`**, **`.bgen`**, **`.vcf.gz`**, …) | Must match **`input_fmt`** (or use **`auto`**) |
| `<regions_dir>` | Step **5** output (folder with **`chr*_regions.txt`**) | Same paths you used in Step 5 |
| `<out_dir>` | Where **`plink2a`** writes per-chr extractions | Created if needed |
| `[threads]` | Parallel threads for extraction | Positive integer; script default if omitted (see script header) |
| `[output_fmt]` | Target format for extracted genotypes | **`bfile`**, **`pgen`**, **`vcf`**, **`bgen`** |
| `[input_fmt]` | Format of inputs under **`genotype_dir`** | **`bfile`**, **`pgen`**, **`vcf`**, **`bgen`**, or **`auto`** to detect |
| `[chr_prefix]` | Filename prefix for chromosomes | e.g. **`chr`**, **`imputed_chr`**, or **`""`** if files are **`1.bed`** … **`22.bed`** |
| `[genes_per_chunk]` | Genes per chunk file | **`0`** = one output set per chromosome; **`20`**–**`50`** typical if RAM-limited → **`chr1_genes_chunk1`** … |

Requires **`plink2a`**. Writes **`extraction_summary.txt`** + **`extraction.log`** under **`out_dir`**.

### Step 7 — `step7_verify_extraction.sh`

```bash
./codes/step7_verify_extraction.sh <extract_dir> <regions_dir> [report.txt] [format]
```

| Argument | Expected input | Options / default |
|----------|----------------|-------------------|
| `<extract_dir>` | Step **6** **`out_dir`** (extracted genotypes) | Must match **`format`** |
| `<regions_dir>` | Same **`regions_dir`** passed to Step **6** | For cross-check vs intended regions |
| `[report.txt]` | Path for text QC report | Optional; default chosen by script if omitted |
| `[format]` | Extracted genotype format | **`bfile`**, **`pgen`**, **`vcf`**, **`bgen`** — **must match Step 6** |

Report: counts per chr/chunk, sizes, optional MAF histogram if **`plink2a`** runs `--freq`.

### Step 8 — SAIGE-GENE (configuration + run)

| Script | Arguments / inputs | Role |
|--------|-------------------|------|
| `step8_pre1_build_saige_config.sh` | Interactive prompts → writes **`key=value`** config path you choose | Build SAIGE config template |
| `step8_pre2_validate_saige_config.sh` | `<config.txt>` — path to **`key=value`** file | Validates paths/options (**`source`s** the file) |
| `step8_run_saige_gene_tests.sh` | `<config.txt>` same validated file | Runs SPAtests / SAIGE-GENE jobs; chunk merge per **`MERGE_CHUNKS`** |

**Config file (`key=value`) — required keys**

| Key | Expected input | Notes |
|-----|----------------|-------|
| **`GENOTYPE_DIR`** | Path to Step **6** genotype extract folder | Same as Step 6 `<out_dir>` |
| **`OUTPUT_DIR`** | SAIGE outputs (will contain **`chr*_combined_results.txt`** when merged) | Writable directory |
| **`GMMAT_MODEL`** | Path to **`.rda`** null model from SAIGE step 1 | From your SAIGE null fit |
| **`VARIANCE_RATIO`** | Path to variance ratio text from SAIGE step 1 | From your SAIGE null fit |
| **`GROUP_FILE`** + **`GROUP_FILE_BY_CHR=no`** | Single merged groups file | **Or** **`GROUP_FILE_BY_CHR=yes`** + **`GROUP_DIR`** pointing at per-chr **`chr*_groups.txt`** |
| **`INPUT_FORMAT`** | **`bfile`**, **`pgen`**, **`bgen`**, **`vcf`** | Must match Step **6** extraction |
| **`THREADS`** | Integer | CPU threads for SAIGE jobs |
| **`CHROMOSOMES`** | Comma list, quoted if needed | e.g. **`1,2,…,22`** |
| **`MERGE_CHUNKS`** / **`CHUNKED_INPUT`** | **`yes`**/**`no`** per your chunking | When **`MERGE_CHUNKS=yes`**, merged **`chr*_combined_results.txt`** |

**Other common keys:** **`r_corr`** (maps to SAIGE **`--r.corr`**); **`maxMAF_in_groupTest`** — comma-separated, **quoted** in file; **`annotation_in_groupTest`** — **colons** between burden masks within the string.

---

### Step 9: Post-Analysis

| Question | Answer |
|----------|--------|
| When to use Step **9**? | All SAIGE chr outputs sit under **one** directory. |
| Results dir default | **`.`** — set **`--results-dir`** / **`STEP9_RESULTS_DIR`** if needed. |
| Input files SAIGE produces | **`chr*_combined_results.txt`**, **`chr*_chunk*_results.txt`**, or **`chr*_all_results.txt`**. |
| Several cohort folders instead? | Use [Step 10](#step-10-batch-step-9-multiple-cohort-folders) — do **not** point Step 9 at multiple cohort trees at once. |

**Optional helper:** `codes/step9_noninteractive_example.sh` — set **`RESULTS_DIR`**, **`PRESET`**, **`STEP9_MANHATTAN_P_MODE`**, etc., then run:

```bash
chmod +x codes/step9_noninteractive_example.sh
RESULTS_DIR=/path/to/saige_out PRESET=quick ./codes/step9_noninteractive_example.sh
# Multi-P QQ/Manhattan: STEP9_MANHATTAN_P_MODE=all PRESET=plots_only ./codes/step9_noninteractive_example.sh
# Equivalent: STEP9_RESULTS_DIR=/path/to/saige_out PRESET=quick ./codes/step9_noninteractive_example.sh
```

**Dependencies for optional features**

| Feature | Needs |
|---------|--------|
| Ensembl coords | **`python3`**, **`certifi`** (macOS TLS) |
| PNG QQ / Manhattan | **`Rscript`**, **`ggplot2`** (**`ggrepel`** optional) |

**CLI** — flags first, then positionals:

```bash
./codes/step9_analyze_results.sh [--results-dir|-d|--dir PATH] OPERATIONS \
  [gene_coords_file] [position_mode] [filter_group] [filter_maf] \
  [coord_source] [ensembl_build] [ensembl_release]
```

| Flag / positional | Expected input | Options / meaning |
|-------------------|----------------|-------------------|
| `--results-dir`, `-d`, `--dir` | Path to **one** SAIGE output directory | Same as env **`STEP9_RESULTS_DIR`**; holds **`chr*_combined_results.txt`** (and related files) |
| `OPERATIONS` | First positional after optional flags | Preset **`standard`**, **`full`**, **`quick`** **or** **`+`**‑joined ops (e.g. `mergeall+findsig+top50`) |
| `gene_coords_file` | Optional tab file with gene coordinates for Manhattan | Typically chr, start, end, gene symbol columns; omit if using sequence-only or Ensembl |
| `position_mode` | Where on the gene interval to plot | **`start`**, **`end`**, **`midpoint`** |
| `filter_group` | Annotation group filter | e.g. **`lof`** or **`lof;missense`**; omit for no filter |
| `filter_maf` | Numeric MAF ceiling | e.g. **`0.01`** → keep tests with **max_MAF ≤ 0.01** |
| `coord_source` | How to place genes without a coord file | **`sequence`** (from SAIGE table) or **`ensembl`** (REST lookup) |
| `ensembl_build` | Genome build for Ensembl | **`37`** or **`38`** when **`coord_source=ensembl`** |
| `ensembl_release` | Ensembl release or REST base URL | e.g. **`115`**; optional |

**Common environment variables (override defaults without editing the script)**

| Variable | Expected input | Typical values |
|----------|----------------|----------------|
| **`STEP9_MANHATTAN_P_MODE`** | Which *P* columns drive QQ/Manhattan / summaries | **`pvalue`** (default), **`burden`**, **`skat`**, **`all`** |
| **`STEP9_CAUCHY_MODE`** | How to treat Cauchy annotation rows | **`off`**, **`plots`**, **`pipeline`**, **`full`** |
| **`STEP9_PLOT_UNFILTERED`** | QQ/Manhattan for all groups × MAF | **`1`** (default) or **`0`** |

**More `STEP9_*` environment variables**

| Topic | Variables |
|-------|-----------|
| **Plots / labels** | **`STEP9_PLOT_UNFILTERED`** (default **`1`**) · **`STEP9_MANHATTAN_LABEL_TOP_N`** · **`STEP9_MANHATTAN_LABEL_EXTRA`** (comma symbols) · **`STEP9_MANHATTAN_LABEL_EXTRA_FILE`** (tab: **`Gene`** / **`Symbol`** / col1) · **`STEP9_BONFERRONI_MAF_TESTS`** (default **`3`**) |
| **Which *P* column** | Tables may have **`Pvalue`**, **`Pvalue_Burden`**, **`Pvalue_SKAT`**. **`STEP9_MANHATTAN_P_MODE`** = **`pvalue`** \| **`burden`** \| **`skat`** \| **`all`**. **`all`** → files suffixed **`_burden`** / **`_skat`**; **`group_statistics.txt`** stays combined **`Pvalue`** only. |
| **Cauchy** | **`STEP9_CAUCHY_MODE`** = **`off`** \| **`plots`** \| **`pipeline`** \| **`full`**. **`all_results.txt`** never edited. Legacy: **`STEP9_EXCLUDE_CAUCHY=1`** ⇒ **`plots`**. |
| **Ensembl tuning** | **`STEP9_ENSEMBL_SSL_VERIFY`** (`0` = insecure). **`STEP9_ENSEMBL_BATCH_SIZE`**, **`STEP9_ENSEMBL_PARALLEL`**, **`STEP9_ENSEMBL_POST_TIMEOUT`** |
| **Reuse lookups** | **`STEP9_REUSE_COORD_LOOKUP`** → saved **`.gene_coords_lookup.txt`**. **`STEP9_REUSE_ENSEMBL_REPORT`** optional. Step **10** sets these on cohort 2+ when args match. |

**Interactive vs non-interactive**

| Mode | Behavior |
|------|----------|
| **Interactive** | Prompts: Manhattan label count, extra symbols, Burden/SKAT **P** mode if columns exist |
| **Non-interactive** | Uses **`STEP9_MANHATTAN_P_MODE`** and env vars above |

**Presets → operations run**

| Preset | Pipeline |
|--------|----------|
| **`standard`** | detect_files → mergeall → listgroups → findsig → top50 → groupsum → makeplots → fullsum |
| **`quick`** | detect_files → mergeall → listgroups → top50 → makeplots → fullsum *(no findsig / groupsum)* |
| **`full`** | detect_files → mergechrom → mergeall → listgroups → findsig → top10 → top50 → top100 → chromsum → groupsum → makeplots → fullsum |

**Single operations (examples)** — chain with **`+`**: `mergechrom`, `mergeall`, `listgroups`, `findsig`, `findgws`, `findsug`, `findnom`, `top10`, `top50`, `top100`, `chromsum`, `groupsum`, `fullsum`, `qqdata`, `mandata`. **`plotdata`** ≡ **`makeplots`**.

**Main output files**

| Category | Files |
|----------|--------|
| **Merged table** | **`all_results.txt`** |
| **Significance tiers** | **`genome_wide_sig*.txt`**, **`suggestive_sig*.txt`**, **`nominal_sig*.txt`**, **`sig_p001*.txt`** (+ **`_*burden`** / **`_*skat`** files when present) |
| **Gene lists** | **`top*_genes.txt`** |
| **Plot data / PNG** | **`qq_plot_data*.txt`**, **`manhattan_plot_data*.txt`**, **`qq_plot*.png`**, **`manhattan_plot*.png`** |
| **Summary** | **`group_statistics.txt`** (makeplots, combined **P**), **`analysis_summary.txt`** (fullsum) |
| **Helpers** | **`.all_results_no_cauchy.txt`**, **`.gene_coords_lookup.txt`**, **`.step9_generate_plots.R`** |

```bash
# Interactive (menu + prompts)
./codes/step9_analyze_results.sh

./codes/step9_analyze_results.sh --results-dir /path/to/saige_out standard
STEP9_RESULTS_DIR=/path/to/out ./codes/step9_analyze_results.sh quick
```

### Step 10 (batch Step 9: multiple cohort folders)

- **Use Step 10:** Several sibling folders under one parent, each with **`chr*_combined_results.txt`**.
- **Skip Step 10:** Only **one** SAIGE results directory → Step **9** alone.

| Situation | Example command |
|-----------|-----------------|
| **One** results dir | `./codes/step9_analyze_results.sh --results-dir saige_results …` |
| **Multiple** dirs | `./codes/step10_run_step9_all_folders.sh …` |

CLI reference: **`./codes/step10_run_step9_all_folders.sh --help`**

**What it does**

| Action | Detail |
|--------|--------|
| Runs | **`step9_analyze_results.sh`** once per cohort subfolder |
| Copies | **`*.png`** + main tables → **`SummaryResults/`** as **`COHORT_filename`** |

| Argument / env | Expected input | Options / meaning |
|----------------|----------------|-------------------|
| **`-i`** / **`--interactive`** | No value | Prompts for Step 9–style options in the terminal; cannot combine with **`-- …`** Step 9 tail |
| **`[COHORT …]`** | Names of child folders under the cohort root | Each folder must contain **`chr*_combined_results.txt`**; omit to auto-discover all qualifying folders |
| **`--`** then Step 9 args | Same order as **`step9_analyze_results.sh`** | Forwarded to every cohort; **do not** pass **`--results-dir`** (the wrapper sets **`STEP9_RESULTS_DIR`** per folder). Omit **`--`** entirely → **`standard`** preset each cohort |
| **`STEP9_COHORT_ROOT`** | Parent directory of **`AFR/`**, **`EUR/`**, … | Default: repo root when the script lives in **`codes/`**, else the script’s directory |
| **`STEP9_SUMMARY_DIR`** | Where to copy prefixed **`COHORT_*`** artifacts | Default **`<cohort_root>/SummaryResults`** |

```bash
chmod +x codes/step10_run_step9_all_folders.sh

# Repository root = parent of cohort folders + SummaryResults (wrapper in codes/ uses parent as cohort root).
./codes/step10_run_step9_all_folders.sh                     # all qualifying cohort folders → step9 standard each
./codes/step10_run_step9_all_folders.sh AFR EUR Pooled      # only these folders

# Non-interactive: step9 args after `--` (omit `--results-dir`; wrapper sets STEP9_RESULTS_DIR per cohort).
./codes/step10_run_step9_all_folders.sh -- quick
./codes/step10_run_step9_all_folders.sh AFR EUR -- standard /path/gene_coords.tsv midpoint lof 0.01

# Interactive (terminal required): shared vs per-cohort parameters.
./codes/step10_run_step9_all_folders.sh --interactive
./codes/step10_run_step9_all_folders.sh -i AFR EUR
```

**Paths**

| Env var | Role | Default |
|---------|------|---------|
| **`STEP9_COHORT_ROOT`** | Parent of **`AFR/`**, **`EUR/`**, … | Repo root if script is under **`codes/`**, else script dir |
| **`STEP9_SUMMARY_DIR`** | Prefixed copies destination | **`<cohort_root>/SummaryResults`** |

**Behavior**

| Topic | Rule |
|-------|------|
| **Default preset** | No **`-- …`** tail → **`standard`** per cohort (same as Step 9 CLI) |
| **After `--`** | Tokens forwarded to **`step9_analyze_results.sh`** every cohort (**omit `--results-dir`**) |
| **Env** | **`STEP9_*`** vars apply to each cohort run |
| **`-i` / `--interactive`** | 2+ cohorts: prompt **same params for all?** — **Y** once, **N** per cohort. **Never** mix **`-i`** with **`-- …`** Step 9 args |

**Performance**

| Mechanism | Effect |
|-----------|--------|
| Step 9 **`fullsum`** per cohort | Single-pass counts; may reuse **`annotation_group_summary`** |
| Coordinate cache | Same Step 9 args → cohort **1** builds lookup; cohort **2+** use **`.step9_batch_coord_cache/<hash>/`** (**`STEP9_COORD_CACHE_DIR`**) |
| Caveat | Sequence-only Manhattan skips cache; missing genes → **NA** on plots |

**Long runs** — log to file:

```bash
nohup ./codes/step10_run_step9_all_folders.sh > step10_all_folders.log 2>&1 &
```

### Step 11 — `step11_saige_group_to_regenie.sh`

Convert SAIGE gene-group files (**`GENE var …`** / **`GENE anno …`** blocks, optionally **`GENE weight …`**) into REGENIE burden-test inputs: tab-separated **annotation** files (**`--anno-file`**), **set-list** files (**`--set-list`**), and an example **mask-definition** file (**`--mask-def`**). **`step11_saige_group_to_regenie.sh`** is bash-only. Align mask tokens with Milind et al.–style LoF masks (e.g. MAF below 1%) when configuring REGENIE.

```bash
./codes/step11_saige_group_to_regenie.sh
```

| Variable / EDIT block | Expected input | Options / default |
|----------|----------------|-------------------|
| `SAIGE_GROUP_DIR` / **`SAIGE_GROUP_DIR_DEFAULT`** | Directory containing per-chromosome **`chr##_group.txt`** or **`chr#_group.txt`** (autosomes **1–22**) | Default **`./saige_groups`** |
| `OUT_DIR` / **`OUT_DIR_DEFAULT`** | Output root for **`chrNN/`** subfolders | Default **`./regenie_masks_from_saige`** |
| **`STRIP_CHR_PREFIX`** | If **`1`**, strip a leading **`chr`** from variant IDs in anno/set-list rows | Default **`0`** |

**Outputs (under `OUT_DIR`)**

| Path | Role |
|------|------|
| **`chrNN/chrNN.regenie.annotations.txt`** | Variant ID ↔ gene/set ↔ annotation token (**`--anno-file`**) |
| **`chrNN/chrNN.regenie.setlist.txt`** | Gene/set ↔ chromosome ↔ min position ↔ comma-separated variant IDs (**`--set-list`**) |
| **`mask_def.example.txt`** | Example **`--mask-def`** template — edit so the second column lists tokens **exactly** as they appear in the anno lines |

Requires **bash** only (no **`awk`**). **Windows CRLF** in the script breaks **`set -o pipefail`** — run **`dos2unix`** on the script if needed.

---

## Important notes

| Topic | Note |
|-------|------|
| Pre-built **`gene_coords_*.tar.gz`** | Can replace Step **4** if Ensembl release matches VEP |
| Step **5** optionals | **`buffer_kb`**, **`force_regen`**, **`merge_regen`** — see **Step 5** section above |
| Step **6** **`genes_per_chunk`** | **`0`** = whole chr per file; **`20`–`50`** if RAM-limited |
| SAIGE config file | Use **`r_corr`** (not **`r.corr`**); **`annotation_in_groupTest`** masks separated by **colons** |

---

## Troubleshooting

#### Windows line endings

```bash
sed -i '' 's/\r$//' codes/*.sh   # GNU sed: sed -i 's/\r$//' codes/*.sh
chmod +x codes/*.sh
```

#### LSF: run Step 1 per chromosome

```bash
for chr in {1..22}; do
  bsub -q normal -n 4 -R "rusage[mem=8GB]" \
    "./codes/step1_vep_ann_clean.sh chr${chr}.vcf.gz chr${chr}_anno.txt 4 vep_lof"
done
```

#### Input file not found

```bash
./codes/step1_vep_ann_clean.sh /full/path/to/input.vcf.gz output.txt
```

#### PLINK2 out of memory

```bash
plink2 --pfile input --memory 16000 --threads 8 ...
```

#### Missing genes (Step 5)

- Match Ensembl version to VEP
- Read **`missing_genes.txt`**
- Try alternate symbols / ENSG IDs

#### Step 6 format auto-detect fails

```bash
./codes/step6_extract_genotypes_plink2a.sh /data/geno plink_regions out 8 vcf vcf chr 0
```

#### PGEN → BGEN for SAIGE

```bash
plink2 --pfile chr1_genes \
       --export bgen-1.2 bits=8 \
       --out chr1_genes_bgen

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

*Version 2.1.2 | Last updated: 2026-05-15*

#!/bin/bash

# Script: step4_download_gene_coords.sh
# Description: Download Ensembl gene coordinates by chromosome
# Usage: step4_download_gene_coords.sh [ensembl_version] [genome_build] [output_dir]
# Internet required

set -euo pipefail

ENSEMBL_VERSION="${1:-115}"
GENOME_BUILD="${2:-GRCh38}"
OUTPUT_DIR="${3:-gene_coords}"

mkdir -p "$OUTPUT_DIR"

echo "[$(date)] Downloading Ensembl gene coordinates..." >&2
echo "  Version: $ENSEMBL_VERSION" >&2
echo "  Build: $GENOME_BUILD" >&2
echo "  Output: $OUTPUT_DIR" >&2

#==========================================
# Detect OS and set commands
#==========================================

if [[ "$OSTYPE" == "darwin"* ]]; then
    ZCAT_CMD="gunzip -c"
    AWK_CMD="awk"
else
    ZCAT_CMD="zcat"
    AWK_CMD="awk"
fi

#==========================================
# Download GTF
#==========================================

if [ "$GENOME_BUILD" = "GRCh38" ]; then
    GTF_URL="ftp://ftp.ensembl.org/pub/release-${ENSEMBL_VERSION}/gtf/homo_sapiens/Homo_sapiens.GRCh38.${ENSEMBL_VERSION}.gtf.gz"
elif [ "$GENOME_BUILD" = "GRCh37" ]; then
    # Filename release tag must match ENSEMBL_VERSION (directory release-N)
    GTF_URL="ftp://ftp.ensembl.org/pub/grch37/release-${ENSEMBL_VERSION}/gtf/homo_sapiens/Homo_sapiens.GRCh37.${ENSEMBL_VERSION}.gtf.gz"
else
    echo "Error: Build must be GRCh38 or GRCh37" >&2
    exit 1
fi

echo "Downloading GTF..." >&2
wget -q --show-progress "$GTF_URL" -O "$OUTPUT_DIR/ensembl.gtf.gz"

#==========================================
# Extract gene coordinates
#==========================================

echo "Extracting gene coordinates..." >&2

$ZCAT_CMD "$OUTPUT_DIR/ensembl.gtf.gz" | \
  $AWK_CMD -F'\t' '
  BEGIN {
    OFS="\t"
    print "Gene_Symbol", "Gene_ID", "Chr", "Start", "End", "Strand", "Biotype"
  }
  $3 == "gene" {
    # Extract gene_id
    gid = ""; gname = ""; gtype = ""
    
    if (match($9, /gene_id "([^"]+)"/)) {
      gid = substr($9, RSTART+9, RLENGTH-10)
    }
    
    # Extract gene_name
    if (match($9, /gene_name "([^"]+)"/)) {
      gname = substr($9, RSTART+11, RLENGTH-12)
    }
    
    # Extract gene_biotype
    if (match($9, /gene_biotype "([^"]+)"/)) {
      gtype = substr($9, RSTART+14, RLENGTH-15)
    }
    
    # Clean chromosome name (remove chr prefix if present)
    chr = $1
    gsub(/^chr/, "", chr)
    
    if (gname != "" && gid != "") {
      print gname, gid, chr, $4, $5, $7, gtype
    }
  }
  ' > "$OUTPUT_DIR/all_genes_coords.txt"

#==========================================
# Split by chromosome
#==========================================

echo "Splitting by chromosome..." >&2

# Create header file
head -1 "$OUTPUT_DIR/all_genes_coords.txt" > "$OUTPUT_DIR/header.txt"

# Split by chromosome (1-22, X, Y, MT)
for chr in {1..22} X Y MT; do
  OUTPUT_FILE="$OUTPUT_DIR/chr${chr}_genes.txt"
  
  # Add header
  cat "$OUTPUT_DIR/header.txt" > "$OUTPUT_FILE"
  
  # Extract chromosome
  $AWK_CMD -F'\t' -v chr="$chr" 'NR>1 && $3==chr' "$OUTPUT_DIR/all_genes_coords.txt" >> "$OUTPUT_FILE"
  
  if [ -s "$OUTPUT_FILE" ]; then
    GENE_COUNT=$(tail -n +2 "$OUTPUT_FILE" | wc -l)
    echo "  chr${chr}: $GENE_COUNT genes" >&2
  else
    rm "$OUTPUT_FILE"
  fi
done

# Cleanup
rm "$OUTPUT_DIR/header.txt"
rm "$OUTPUT_DIR/ensembl.gtf.gz"

#==========================================
# Summary and packaging
#==========================================

TOTAL_GENES=$(tail -n +2 "$OUTPUT_DIR/all_genes_coords.txt" | wc -l)

echo "" >&2
echo "==========================================" >&2
echo "Download Complete" >&2
echo "==========================================" >&2
echo "Total genes: $TOTAL_GENES" >&2
echo "Output directory: $OUTPUT_DIR" >&2
echo "" >&2
echo "Files created:" >&2
ls -1 "$OUTPUT_DIR"/chr*.txt 2>/dev/null | while read f; do
  count=$(tail -n +2 "$f" | wc -l)
  echo "  $(basename $f): $count genes" >&2
done
echo "" >&2
echo "Transfer package:" >&2
echo "  tar -czf gene_coords_ensembl${ENSEMBL_VERSION}.tar.gz $OUTPUT_DIR" >&2
echo "==========================================" >&2

# Create transfer package
tar -czf "gene_coords_ensembl${ENSEMBL_VERSION}.tar.gz" "$OUTPUT_DIR"

echo "" >&2
echo "Ready to transfer:" >&2
echo "  gene_coords_ensembl${ENSEMBL_VERSION}.tar.gz" >&2
echo "  Size: $(du -h gene_coords_ensembl${ENSEMBL_VERSION}.tar.gz | cut -f1)" >&2

exit 0

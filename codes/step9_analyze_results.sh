#!/bin/bash

# step9_analyze_results.sh
# Interactive analysis script for SAIGE-GENE results
# Enhanced version with gene coordinate support and annotation filtering

set -euo pipefail

# Where chromosome / merged SAIGE output files live (--results-dir PATH, STEP9_RESULTS_DIR, or '.').
RESULT_FOLDER="${STEP9_RESULTS_DIR:-.}"

_remain=()
while [ $# -gt 0 ]; do
    case "$1" in
        -d|--results-dir|--dir)
            if [ $# -lt 2 ]; then
                echo "ERROR: $1 requires a directory argument." >&2
                exit 2
            fi
            RESULT_FOLDER="$2"
            shift 2
            ;;
        *)
            _remain+=("$1")
            shift
            ;;
    esac
done

# Avoid unbound/expansion issues with set -u and an empty ${_remain[@]}.
if [ "${#_remain[@]}" -gt 0 ]; then
    set -- "${_remain[@]}"
else
    set --
fi

OUTPUT_DIR=""

finalize_results_directory() {
    local rf="${1:-.}"
    if [ -z "$rf" ]; then
        rf="."
    fi
    if [ ! -d "$rf" ]; then
        echo "ERROR: Results folder is not a directory: $rf" >&2
        exit 1
    fi
    OUTPUT_DIR="$(cd "$rf" && pwd)"
}

print_startup_banner() {
    echo "=========================================="
    echo "SAIGE-GENE Results Analysis Tool"
    echo "Enhanced Version with Annotations"
    echo "=========================================="
    echo "Started: $(date)"
    echo "Results directory: $OUTPUT_DIR"
    echo ""
}

GENE_COORD_FILE=""
POSITION_MODE="midpoint"  # Options: start, end, midpoint
COORD_SOURCE="sequence"   # Options: sequence, ensembl
ENSEMBL_SERVER="https://rest.ensembl.org"
ENSEMBL_BUILD="38"        # Options: 37, 38
ENSEMBL_RELEASE=""        # Optional release number (e.g., 115)
FILTER_GROUP=""  # Filter by specific annotation group
FILTER_MAF=""    # Filter by max_MAF threshold

# Plot outputs (makeplots / plotdata): set via STEP9_PLOT_UNFILTERED or interactive prompts (see Main Execution).
PLOT_DO_UNFILTERED=""

#==========================================
# Plot filename helpers
#==========================================

plot_slug_group() {
    local s
    s=$(echo "$1" | sed 's/[\/;]/_/g' | tr ' ' '_' | tr -cd 'A-Za-z0-9_.-')
    if [ -z "$s" ]; then
        echo "group"
    else
        echo "$s"
    fi
}

plot_slug_maf() {
    echo "$1" | sed 's/^0\./0p/' | tr -d ' '
}

#==========================================
# Display Banner
#==========================================

#==========================================
# Parse Command Line Arguments
#==========================================

if [ $# -eq 0 ]; then
    # Interactive mode
    INTERACTIVE=true
else
    # Command mode
    INTERACTIVE=false
    OPERATIONS="$1"
    
    # Check for gene coordinate file
    if [ $# -ge 2 ]; then
        GENE_COORD_FILE="$2"
    fi
    
    # Check for position mode
    if [ $# -ge 3 ]; then
        POSITION_MODE="$3"
    fi
    
    # Check for group filter
    if [ $# -ge 4 ]; then
        FILTER_GROUP="$4"
    fi
    
    # Check for MAF filter
    if [ $# -ge 5 ]; then
        FILTER_MAF="$5"
    fi
    
    # Optional coordinate source override: sequence|ensembl
    if [ $# -ge 6 ]; then
        COORD_SOURCE="$6"
    fi
    
    # Optional Ensembl build: 37|38
    if [ $# -ge 7 ]; then
        ENSEMBL_BUILD="$7"
    fi
    
    # Optional Ensembl release (e.g. 115) or custom URL
    if [ $# -ge 8 ]; then
        ENSEMBL_RELEASE="$8"
    fi
fi

# Non-interactive / flag-based runs resolve the results folder now.
if [ "$INTERACTIVE" != true ]; then
    finalize_results_directory "$RESULT_FOLDER"
    print_startup_banner
fi

#==========================================
# Detect Available Annotation Groups
#==========================================

detect_groups() {
    local file="$1"
    
    if [ ! -f "$file" ]; then
        echo "  Warning: File $file not found"
        return 1
    fi
    
    echo "Available annotation groups in data:"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    
    tail -n +2 "$file" | cut -f2 | sort -u | while read -r group; do
        count=$(tail -n +2 "$file" | awk -F'\t' -v g="$group" '$2 == g' | wc -l)
        printf "  %-30s %8d genes\n" "$group" "$count"
    done
    echo ""
}

#==========================================
# Gene Coordinate Functions
#==========================================

setup_gene_coords() {
    if [ -n "$GENE_COORD_FILE" ] && [ -f "$GENE_COORD_FILE" ]; then
        echo "Gene coordinate file: $GENE_COORD_FILE"
        echo "Position mode: $POSITION_MODE"
        
        # Create lookup file with chr and position
        awk -v mode="$POSITION_MODE" 'NR > 1 {
            chr = $1
            start = $2
            end = $3
            gene = $4
            
            if (mode == "start") {
                pos = start
            } else if (mode == "end") {
                pos = end
            } else {
                pos = int((start + end) / 2)
            }
            
            print gene "\t" chr "\t" pos
        }' "$GENE_COORD_FILE" > "$OUTPUT_DIR/.gene_coords_lookup.txt"
        
        GENE_COUNT=$(wc -l < "$OUTPUT_DIR/.gene_coords_lookup.txt")
        echo "  Loaded coordinates for $GENE_COUNT genes"
        echo ""
        return 0
    else
        if [ -n "$GENE_COORD_FILE" ]; then
            echo "Warning: Gene coordinate file not found: $GENE_COORD_FILE"
        fi
        echo "  Will extract chromosome from filenames"
        echo ""
        return 1
    fi
}

setup_gene_coords_ensembl() {
    local all_results_file="$OUTPUT_DIR/all_results.txt"
    local lookup_file="$OUTPUT_DIR/.gene_coords_lookup.txt"
    local tmp_report="$OUTPUT_DIR/.ensembl_lookup_report.txt"
    
    # Configure server from selected build/release.
    if [[ "$ENSEMBL_RELEASE" =~ ^https?:// ]]; then
        ENSEMBL_SERVER="$ENSEMBL_RELEASE"
    else
        if [ "$ENSEMBL_BUILD" = "37" ]; then
            ENSEMBL_SERVER="https://grch37.rest.ensembl.org"
        else
            ENSEMBL_SERVER="https://rest.ensembl.org"
        fi
    fi
    
    echo "Using Ensembl coordinates"
    echo "  If SSL fails on macOS/Python: pip3 install certifi, "
    echo "    or STEP9_ENSEMBL_SSL_VERIFY=0 (not recommended for production)."
    echo "  Genome build: GRCh$ENSEMBL_BUILD"
    if [ -n "$ENSEMBL_RELEASE" ] && [[ ! "$ENSEMBL_RELEASE" =~ ^https?:// ]]; then
        echo "  Requested release: $ENSEMBL_RELEASE"
        echo "  Note: release number is recorded for tracking; exact archive routing requires a custom Ensembl URL."
    fi
    echo "  Server: $ENSEMBL_SERVER"
    echo "  Position mode: $POSITION_MODE"
    echo "  Speed: STEP9_ENSEMBL_BATCH_SIZE (default 500), STEP9_ENSEMBL_PARALLEL workers (default 4)"
    
    # Check what releases are actually available on selected server.
    python3 - "$ENSEMBL_SERVER" "$ENSEMBL_RELEASE" << 'PY'
import json
import os
import ssl
import sys
import urllib.error
import urllib.request


def ssl_context():
    insecure = os.environ.get("STEP9_ENSEMBL_SSL_VERIFY", "1").strip().lower() in ("0", "false", "no")
    if insecure:
        print(
            "  Warning: STEP9_ENSEMBL_SSL_VERIFY=0 — HTTPS certificate verification disabled.",
            flush=True,
        )
        ctx = ssl.create_default_context()
        ctx.check_hostname = False
        ctx.verify_mode = ssl.CERT_NONE
        return ctx
    try:
        import certifi

        return ssl.create_default_context(cafile=certifi.where())
    except Exception:
        pass
    for path in (
        "/etc/ssl/cert.pem",
        "/opt/homebrew/etc/openssl@3/cert.pem",
        "/usr/local/etc/openssl@3/cert.pem",
        "/etc/pki/tls/certs/ca-bundle.crt",
        "/etc/ssl/certs/ca-certificates.crt",
    ):
        try:
            if os.path.isfile(path):
                return ssl.create_default_context(cafile=path)
        except Exception:
            pass
    return ssl.create_default_context()


_SSL = ssl_context()


def _ssl_recoverable(err):
    if isinstance(err, ssl.SSLCertVerificationError):
        return True
    msg = " ".join(
        [
            repr(err),
            str(err),
            repr(getattr(err, "reason", "")),
            str(getattr(err, "reason", "")),
        ]
    ).upper()
    return "CERTIFICATE_VERIFY_FAILED" in msg or "CERTIFICATE VERIFY FAILED" in msg or "VERIFY FAILED: UNABLE TO GET LOCAL ISSUER" in msg


def urlopen_https(req, timeout):
    global _SSL
    try:
        return urllib.request.urlopen(req, timeout=timeout, context=_SSL)
    except Exception as e:
        if not _ssl_recoverable(e):
            raise
    try:
        import certifi

        _SSL = ssl.create_default_context(cafile=certifi.where())
        return urllib.request.urlopen(req, timeout=timeout, context=_SSL)
    except Exception:
        pass
    for path in (
        "/etc/ssl/cert.pem",
        "/opt/homebrew/etc/openssl@3/cert.pem",
        "/usr/local/etc/openssl@3/cert.pem",
    ):
        if os.path.isfile(path):
            try:
                _SSL = ssl.create_default_context(cafile=path)
                return urllib.request.urlopen(req, timeout=timeout, context=_SSL)
            except Exception:
                continue
    insecure = (
        "\n"
        + "  Hint: SSL verification failed. Try: pip3 install certifi\n"
        + "       Or run macOS installer: Applications/Python → Install Certificates.command\n"
        + "       Or temporarily: export STEP9_ENSEMBL_SSL_VERIFY=0 (not recommended)"
    )
    print(f"  Could not query available releases from server.{insecure} Continuing anyway.")
    sys.exit(0)


server = sys.argv[1].rstrip("/")
requested = sys.argv[2].strip()
url = f"{server}/info/data?content-type=application/json"
req = urllib.request.Request(url, headers={"Accept": "application/json", "User-Agent": "step9-analyze-results"})

try:
    with urlopen_https(req, timeout=20) as resp:
        data = json.loads(resp.read().decode("utf-8"))
except Exception as e:
    print(f"  Could not query available releases from server ({e}). Continuing anyway.")
    sys.exit(0)

releases = data.get("releases", [])
if isinstance(releases, list) and releases:
    releases_txt = ", ".join(str(r) for r in releases)
    print(f"  Available release(s) on this server: {releases_txt}")
    if requested and requested.isdigit() and int(requested) not in [int(r) for r in releases if str(r).isdigit()]:
        print(f"  Warning: Requested release {requested} is not listed on this server.")
else:
    print("  Server did not report release list via info/data.")
PY
    
    if [ ! -f "$all_results_file" ]; then
        echo "  all_results.txt not found; generating it first..."
        detect_files
        op_mergeall
    fi
    
    if [ ! -f "$all_results_file" ]; then
        echo "  ERROR: Could not prepare all_results.txt for Ensembl lookup"
        return 1
    fi
    
    get_columns "$all_results_file" > /dev/null 2>&1 || {
        echo "  ERROR: Could not detect required columns in all_results.txt"
        return 1
    }
    
    if ! command -v python3 >/dev/null 2>&1; then
        echo "  ERROR: python3 is required for Ensembl coordinate lookup"
        return 1
    fi
    
    python3 - "$all_results_file" "$lookup_file" "$tmp_report" "$POSITION_MODE" "$ENSEMBL_SERVER" "$REGIONCOL" << 'PY'
import concurrent.futures
import csv
import json
import os
import ssl
import sys
import time
import urllib.error
import urllib.parse
import urllib.request


def ssl_context():
    insecure = os.environ.get("STEP9_ENSEMBL_SSL_VERIFY", "1").strip().lower() in ("0", "false", "no")
    if insecure:
        print(
            "  Warning: STEP9_ENSEMBL_SSL_VERIFY=0 — HTTPS certificate verification disabled.",
            flush=True,
        )
        ctx = ssl.create_default_context()
        ctx.check_hostname = False
        ctx.verify_mode = ssl.CERT_NONE
        return ctx
    try:
        import certifi

        return ssl.create_default_context(cafile=certifi.where())
    except Exception:
        pass
    for path in (
        "/etc/ssl/cert.pem",
        "/opt/homebrew/etc/openssl@3/cert.pem",
        "/usr/local/etc/openssl@3/cert.pem",
        "/etc/pki/tls/certs/ca-bundle.crt",
        "/etc/ssl/certs/ca-certificates.crt",
    ):
        try:
            if os.path.isfile(path):
                return ssl.create_default_context(cafile=path)
        except Exception:
            pass
    return ssl.create_default_context()


_SSL = ssl_context()


def _ssl_recoverable(err):
    if isinstance(err, ssl.SSLCertVerificationError):
        return True
    msg = " ".join(
        [
            repr(err),
            str(err),
            repr(getattr(err, "reason", "")),
            str(getattr(err, "reason", "")),
        ]
    ).upper()
    return "CERTIFICATE_VERIFY_FAILED" in msg or "CERTIFICATE VERIFY FAILED" in msg or "VERIFY FAILED: UNABLE TO GET LOCAL ISSUER" in msg


def urlopen_https(req, timeout):
    global _SSL
    try:
        return urllib.request.urlopen(req, timeout=timeout, context=_SSL)
    except Exception as e:
        if not _ssl_recoverable(e):
            raise
    try:
        import certifi

        _SSL = ssl.create_default_context(cafile=certifi.where())
        return urllib.request.urlopen(req, timeout=timeout, context=_SSL)
    except Exception:
        pass
    for path in (
        "/etc/ssl/cert.pem",
        "/opt/homebrew/etc/openssl@3/cert.pem",
        "/usr/local/etc/openssl@3/cert.pem",
    ):
        if os.path.isfile(path):
            try:
                _SSL = ssl.create_default_context(cafile=path)
                return urllib.request.urlopen(req, timeout=timeout, context=_SSL)
            except Exception:
                continue
    raise RuntimeError(
        "HTTPS certificate verification failed. Try: pip3 install certifi, "
        "or macOS: Install Certificates.command, "
        "or STEP9_ENSEMBL_SSL_VERIFY=0 (not recommended)."
    )

all_results_file, lookup_file, report_file, position_mode, ensembl_server, region_col = sys.argv[1:]
region_idx = int(region_col) - 1

def extract_gene(region: str) -> str:
    if ":" in region:
        token = region.split(":")[-1]
        if token.isdigit():
            return region
        return token
    return region

ordered_genes = []
seen = set()
with open(all_results_file, "r", newline="") as f:
    reader = csv.reader(f, delimiter="\t")
    header = next(reader, None)
    for row in reader:
        if region_idx >= len(row):
            continue
        gene = extract_gene(row[region_idx]).strip()
        if not gene:
            continue
        if gene not in seen:
            seen.add(gene)
            ordered_genes.append(gene)

def _extract_coord(data: dict):
    if not isinstance(data, dict):
        return None
    if any(k not in data for k in ("seq_region_name", "start", "end")):
        return None
    chr_raw = str(data["seq_region_name"]).replace("chr", "")
    start = int(data["start"])
    end = int(data["end"])
    if position_mode == "start":
        pos = start
    elif position_mode == "end":
        pos = end
    else:
        pos = (start + end) // 2
    return (chr_raw, pos)

def fetch_gene_single(gene: str):
    safe_gene = urllib.parse.quote(gene, safe="")
    url = f"{ensembl_server}/lookup/symbol/homo_sapiens/{safe_gene}?content-type=application/json"
    req = urllib.request.Request(
        url,
        headers={
            "Accept": "application/json",
            "Content-Type": "application/json",
            "User-Agent": "step9-analyze-results",
        },
    )
    try:
        with urlopen_https(req, timeout=15) as resp:
            if resp.status != 200:
                return None
            data = json.loads(resp.read().decode("utf-8"))
    except Exception:
        return None
    return _extract_coord(data)

def fetch_batch(genes):
    url = f"{ensembl_server}/lookup/symbol/homo_sapiens?content-type=application/json"
    payload = json.dumps({"symbols": genes}).encode("utf-8")
    req = urllib.request.Request(
        url,
        data=payload,
        method="POST",
        headers={
            "Accept": "application/json",
            "Content-Type": "application/json",
            "User-Agent": "step9-analyze-results",
        },
    )
    out = {}
    try:
        post_timeout = int(os.environ.get("STEP9_ENSEMBL_POST_TIMEOUT", "60").strip())
    except Exception:
        post_timeout = 60
    post_timeout = max(15, min(300, post_timeout))
    try:
        with urlopen_https(req, timeout=post_timeout) as resp:
            if resp.status != 200:
                return out
            data = json.loads(resp.read().decode("utf-8"))
    except Exception:
        data = {}
    if isinstance(data, dict):
        for g in genes:
            c = _extract_coord(data.get(g))
            if c is not None:
                out[g] = c
    return out

def fetch_chunk_complete(chunk: list) -> dict:
    """One POST batch plus per-gene GET only for symbols the batch JSON skipped."""
    batch_coords = fetch_batch(chunk)
    out = {}
    for gene, val in batch_coords.items():
        out[gene] = {"chr": val[0], "pos": val[1], "source": "ensembl"}
    for gene in chunk:
        if gene in out:
            continue
        val = fetch_gene_single(gene)
        if val is not None:
            out[gene] = {"chr": val[0], "pos": val[1], "source": "ensembl"}
    return out

def _env_int(name: str, default: int) -> int:
    try:
        v = int(os.environ.get(name, str(default)).strip())
        return v if v > 0 else default
    except Exception:
        return default

BATCH_SIZE = max(10, _env_int("STEP9_ENSEMBL_BATCH_SIZE", 500))
MAX_WORKERS = max(1, _env_int("STEP9_ENSEMBL_PARALLEL", 4))
total_genes = len(ordered_genes)
coords = {}
chunks = [ordered_genes[i : i + BATCH_SIZE] for i in range(0, len(ordered_genes), BATCH_SIZE)]
n_chunks = len(chunks)
print(
    f"[Ensembl] Starting batched lookup for {total_genes} genes "
    f"({n_chunks} requests, batch size {BATCH_SIZE}, up to {min(MAX_WORKERS, n_chunks or 1)} parallel)",
    flush=True,
)
start_ts = time.time()
if n_chunks == 0:
    pass
elif MAX_WORKERS <= 1:
    processed = 0
    for ch in chunks:
        coords.update(fetch_chunk_complete(ch))
        processed += len(ch)
        elapsed = max(time.time() - start_ts, 1e-6)
        rate = processed / elapsed
        remaining = max(total_genes - processed, 0)
        eta_seconds = int(remaining / rate) if rate > 0 else -1
        eta_txt = f"{eta_seconds}s" if eta_seconds >= 0 else "NA"
        print(
            f"[Ensembl] Processed {processed}/{total_genes}; mapped so far: {len(coords)}; ETA: {eta_txt}",
            flush=True,
        )
else:
    workers = min(MAX_WORKERS, n_chunks)
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as ex:
        fut_to_meta = {ex.submit(fetch_chunk_complete, ch): len(ch) for ch in chunks}
        processed = 0
        for fut in concurrent.futures.as_completed(fut_to_meta):
            processed += fut_to_meta[fut]
            coords.update(fut.result())
            elapsed = max(time.time() - start_ts, 1e-6)
            rate = processed / elapsed
            remaining = max(total_genes - processed, 0)
            eta_seconds = int(remaining / rate) if rate > 0 else -1
            eta_txt = f"{eta_seconds}s" if eta_seconds >= 0 else "NA"
            print(
                f"[Ensembl] Processed {processed}/{total_genes}; mapped so far: {len(coords)}; ETA: {eta_txt}",
                flush=True,
            )

missing = [g for g in ordered_genes if g not in coords]
for i, gene in enumerate(ordered_genes):
    if gene in coords:
        continue
    prev = None
    nxt = None
    j = i - 1
    while j >= 0:
        g = ordered_genes[j]
        if g in coords:
            prev = coords[g]
            break
        j -= 1
    j = i + 1
    while j < len(ordered_genes):
        g = ordered_genes[j]
        if g in coords:
            nxt = coords[g]
            break
        j += 1
    if prev and nxt and prev["chr"] == nxt["chr"]:
        coords[gene] = {"chr": prev["chr"], "pos": (prev["pos"] + nxt["pos"]) // 2, "source": "neighbor_interp"}
    elif prev:
        coords[gene] = {"chr": prev["chr"], "pos": prev["pos"], "source": "neighbor_prev"}
    elif nxt:
        coords[gene] = {"chr": nxt["chr"], "pos": nxt["pos"], "source": "neighbor_next"}

with open(lookup_file, "w", newline="") as out:
    for gene in ordered_genes:
        c = coords.get(gene)
        if c:
            out.write(f"{gene}\t{c['chr']}\t{c['pos']}\n")

resolved = sum(1 for g in ordered_genes if g in coords)
interpolated = sum(1 for g in ordered_genes if g in coords and coords[g]["source"] != "ensembl")
unresolved = len(ordered_genes) - resolved

with open(report_file, "w") as rep:
    rep.write(f"TOTAL_GENES={len(ordered_genes)}\n")
    rep.write(f"ENSEMBL_FOUND={len(ordered_genes) - len(missing)}\n")
    rep.write(f"INTERPOLATED={interpolated}\n")
    rep.write(f"UNRESOLVED={unresolved}\n")
PY
    
    if [ ! -f "$lookup_file" ]; then
        echo "  ERROR: Ensembl coordinate lookup failed"
        return 1
    fi
    
    # shellcheck disable=SC1090
    source "$tmp_report"
    echo "  Loaded coordinates for $TOTAL_GENES genes"
    echo "  Found in Ensembl: $ENSEMBL_FOUND"
    echo "  Filled from neighboring genes: $INTERPOLATED"
    echo "  Still unresolved: $UNRESOLVED"
    echo ""
    return 0
}

extract_gene_name() {
    local region="$1"
    
    # Region format examples:
    # chr1:12345-67890
    # GENE_NAME
    # chr1:12345-67890:GENE_NAME
    
    # If contains colon, extract gene name after last colon or before first colon
    if [[ "$region" == *":"* ]]; then
        # Try to get gene name after last colon
        gene=$(echo "$region" | awk -F: '{print $NF}')
        # If that looks like a number, use the whole region
        if [[ "$gene" =~ ^[0-9]+$ ]]; then
            gene="$region"
        fi
    else
        gene="$region"
    fi
    
    echo "$gene"
}

extract_chr_from_filename() {
    local filename="$1"
    
    # Extract chromosome from filename like chr1_combined.txt or chr22_chunk1_results.txt
    if [[ "$filename" =~ chr([0-9]+|X|Y) ]]; then
        echo "${BASH_REMATCH[1]}"
    else
        echo "NA"
    fi
}

annotate_with_coords() {
    local infile="$1"
    local outfile="$2"
    local lookup="$OUTPUT_DIR/.gene_coords_lookup.txt"
    
    if [ ! -f "$lookup" ]; then
        echo "  Warning: Gene coordinates not available"
        cp "$infile" "$outfile"
        return 1
    fi
    
    # Ensure Region column index (same logic as extract_gene_name / Region field).
    if [ -z "${REGIONCOL:-}" ] && [ -f "$infile" ]; then
        get_columns "$infile" >/dev/null 2>&1 || true
    fi
    local rcol="${REGIONCOL:-1}"
    
    echo "  Annotating with gene coordinates (single-pass join)..."
    
    local tmp_out="${outfile}.tmp.$$"
    awk -v lookup="$lookup" -v rcol="$rcol" '
        function gene_from_region(r,    n, t, a) {
            if (index(r, ":") == 0)
                return r
            n = split(r, a, /:/)
            t = a[n]
            if (t ~ /^[0-9]+$/)
                return r
            return t
        }
        BEGIN {
            FS = "\t"
            while ((getline line < lookup) > 0) {
                split(line, a, "\t")
                if (a[1] != "") {
                    _chr[a[1]] = a[2]
                    _pos[a[1]] = a[3]
                }
            }
            close(lookup)
        }
        NR == 1 {
            print "CHR\tPOS\t" $0
            next
        }
        {
            region = $rcol
            g = gene_from_region(region)
            if (g in _chr) {
                print _chr[g] "\t" _pos[g] "\t" $0
            } else {
                print "NA\tNA\t" $0
            }
        }
    ' "$infile" > "$tmp_out"
    awk_ec=$?
    if [ "$awk_ec" -ne 0 ]; then
        rm -f "$tmp_out"
        echo "  ERROR: Failed while annotating coordinates (awk exit $awk_ec)." >&2
        return 1
    fi
    mv -f "$tmp_out" "$outfile"
    
    ANNOTATED=$(tail -n +2 "$outfile" | awk -F'\t' '$1 != "NA"' | wc -l)
    TOTAL=$(tail -n +2 "$outfile" | wc -l)
    echo "  Annotated $ANNOTATED of $TOTAL genes with coordinates"
}

#==========================================
# Available Operations
#==========================================

show_menu() {
    echo "Available Operations:"
    echo "===================="
    echo ""
    echo "Data Combination:"
    echo "  mergechrom     - Merge chunked results by chromosome"
    echo "  mergeall       - Combine all chromosomes into single file"
    echo "  listgroups     - List all available annotation groups"
    echo ""
    echo "Significance Filtering:"
    echo "  findsig        - Extract significant results (multiple thresholds)"
    echo "  findgws        - Extract genome-wide significant only (p < 5e-8)"
    echo "  findsug        - Extract suggestive results (p < 1e-5)"
    echo "  findnom        - Extract nominal results (p < 0.05)"
    echo ""
    echo "Gene Ranking:"
    echo "  top10          - Extract top 10 genes"
    echo "  top50          - Extract top 50 genes"
    echo "  top100         - Extract top 100 genes"
    echo ""
    echo "Summary Statistics:"
    echo "  chromsum       - Per-chromosome summary table"
    echo "  groupsum       - Per-annotation-group summary table"
    echo "  fullsum        - Complete summary report"
    echo ""
    echo "Visualization Data:"
    echo "  qqdata         - Generate QQ plot data"
    echo "  mandata        - Generate Manhattan plot data"
    echo "  plotdata       - Generate both QQ and Manhattan data"
    echo "  makeplots      - Generate PNG plots + group statistics table"
    echo ""
    echo "Comprehensive Analysis:"
    echo "  standard       - Standard analysis (mergeall+findsig+top50+fullsum)"
    echo "  full           - Full analysis (all operations)"
    echo "  quick          - Quick analysis (mergeall+top50+fullsum)"
    echo ""
    echo "Usage Examples:"
    echo "  ./step9_analyze_results.sh"
    echo "  ./step9_analyze_results.sh standard"
    echo "  ./step9_analyze_results.sh --results-dir /path/to/results standard"
    echo "  STEP9_RESULTS_DIR=/path/to/results ./step9_analyze_results.sh standard"
    echo "  ./step9_analyze_results.sh standard gene_coords.txt midpoint"
    echo "  ./step9_analyze_results.sh standard gene_coords.txt midpoint 'LOF' 0.01"
    echo "  ./step9_analyze_results.sh standard '' midpoint '' '' ensembl 37 115"
    echo "  ./step9_analyze_results.sh mergeall+listgroups+findsig"
    echo ""
    echo "Results folder (chromosome *.txt outputs):"
    echo "  Default is the current directory. Use:"
    echo "    --results-dir PATH  or  -d PATH"
    echo "    env STEP9_RESULTS_DIR=PATH"
    echo ""
    echo "Ensembl REST speed (coordinate source=ensembl):"
    echo "  STEP9_ENSEMBL_BATCH_SIZE=N   symbols per POST (default 500)"
    echo "  STEP9_ENSEMBL_PARALLEL=N     concurrent batch workers (default 4; set 1 to serialize)"
    echo "  STEP9_ENSEMBL_POST_TIMEOUT=N POST timeout seconds (default 60, max 300)"
    echo ""
    echo "Plot PNG options (makeplots / plotdata; non-interactive):"
    echo "  STEP9_PLOT_UNFILTERED=1|0   include qq + manhattan (all groups, all max_MAF; default 1)"
    echo "  STEP9_BONFERRONI_MAF_TESTS=N   Manhattan Bonferroni divisor for 0.05/(genes×N), default 3"
    echo "  STEP9_MANHATTAN_LABEL_TOP_N=N   label top N genes on Manhattan (default 10; interactive can raise)"
    echo "  STEP9_MANHATTAN_FDR_ALPHA=N      BH-FDR alpha for manhattan_*_fdr.png (default 0.1)"
    echo ""
    echo "Arguments:"
    echo "  1. Operations (required in non-interactive mode)"
    echo "  2. Gene coordinates file (optional)"
    echo "  3. Position mode: start/end/midpoint (default: midpoint)"
    echo "  4. Filter by annotation group (optional)"
    echo "  5. Filter by max_MAF threshold (optional, e.g., 0.01)"
    echo "  6. Coordinate source: sequence or ensembl (optional)"
    echo "  7. Ensembl genome build: 37 or 38 (optional, used if arg6=ensembl)"
    echo "  8. Ensembl release (e.g. 115) or custom URL (optional)"
    echo ""
    echo "Position Modes:"
    echo "  start          - Use gene start position"
    echo "  end            - Use gene end position"
    echo "  midpoint       - Use gene midpoint (default)"
    echo ""
}

#==========================================
# Detect Result Files
#==========================================

detect_files() {
    echo "Detecting result files..."
    
    if ls "$OUTPUT_DIR"/chr*_combined_results.txt 1> /dev/null 2>&1; then
        RESULT_TYPE="combined"
        RESULT_PATTERN="chr*_combined_results.txt"
        echo "  Found combined results (chunks were merged)"
    elif ls "$OUTPUT_DIR"/chr*_chunk*_results.txt 1> /dev/null 2>&1; then
        RESULT_TYPE="chunked"
        RESULT_PATTERN="chr*_chunk*_results.txt"
        echo "  Found chunked results"
    elif ls "$OUTPUT_DIR"/chr*_all_results.txt 1> /dev/null 2>&1; then
        RESULT_TYPE="all"
        RESULT_PATTERN="chr*_all_results.txt"
        echo "  Found chromosome-level results"
    else
        echo "ERROR: No result files found in $OUTPUT_DIR" >&2
        exit 1
    fi
    
    echo ""
}

#==========================================
# Determine Column Numbers
#==========================================

get_columns() {
    local file="$1"
    
    if [ ! -f "$file" ]; then
        echo "ERROR: File $file not found" >&2
        return 1
    fi
    
    HEADER=$(head -1 "$file")
    
    # Find Region column (gene name)
    REGIONCOL=$(echo "$HEADER" | tr '\t' '\n' | grep -n -E '^Region$' | cut -d: -f1 | head -1)
    
    # Find Group column (annotation)
    GROUPCOL=$(echo "$HEADER" | tr '\t' '\n' | grep -n -E '^Group$' | cut -d: -f1 | head -1)
    
    # Find Pvalue column (combined p-value)
    PCOL=$(echo "$HEADER" | tr '\t' '\n' | grep -n -E '^Pvalue$' | cut -d: -f1 | head -1)
    
    # Find max_MAF column
    MAFCOL=$(echo "$HEADER" | tr '\t' '\n' | grep -n -E '^max_MAF$' | cut -d: -f1 | head -1)
    
    # Find MAC column
    MACCOL=$(echo "$HEADER" | tr '\t' '\n' | grep -n -E '^MAC$' | cut -d: -f1 | head -1)
    
    if [ -z "$REGIONCOL" ]; then
        echo "  ERROR: Could not find 'Region' column" >&2
        return 1
    fi
    
    if [ -z "$GROUPCOL" ]; then
        echo "  ERROR: Could not find 'Group' column" >&2
        return 1
    fi
    
    if [ -z "$PCOL" ]; then
        echo "  ERROR: Could not find 'Pvalue' column" >&2
        return 1
    fi
    
    echo "  Region column: $REGIONCOL (gene name)"
    echo "  Group column: $GROUPCOL (annotation)"
    echo "  P-value column: $PCOL (Pvalue)"
    [ -n "$MAFCOL" ] && echo "  max_MAF column: $MAFCOL"
    [ -n "$MACCOL" ] && echo "  MAC column: $MACCOL"
}

#==========================================
# Filter by Group and MAF
#==========================================

apply_filters() {
    local infile="$1"
    local outfile="$2"
    local group="$3"
    local maf="$4"
    
    # If no filters, just copy
    if [ -z "$group" ] && [ -z "$maf" ]; then
        cp "$infile" "$outfile"
        return 0
    fi
    
    echo "  Applying filters:"
    [ -n "$group" ] && echo "    Annotation group: $group"
    [ -n "$maf" ] && echo "    max_MAF <= $maf"
    
    head -1 "$infile" > "$outfile"
    
    # Build awk filter
    local awk_filter='NR > 1'
    
    if [ -n "$group" ] && [ -n "$GROUPCOL" ]; then
        awk_filter="$awk_filter && \$$GROUPCOL == \"$group\""
    fi
    
    if [ -n "$maf" ] && [ -n "$MAFCOL" ]; then
        awk_filter="$awk_filter && \$$MAFCOL != \"NA\" && \$$MAFCOL <= $maf"
    fi
    
    tail -n +2 "$infile" | awk -F'\t' "$awk_filter" >> "$outfile"
    
    FILTERED_COUNT=$(tail -n +2 "$outfile" | wc -l)
    echo "  Filtered to $FILTERED_COUNT genes"
}

#==========================================
# Operation: List Groups
#==========================================

op_listgroups() {
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "Listing annotation groups..."
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    
    if [ ! -f "$OUTPUT_DIR/all_results.txt" ]; then
        echo "  ERROR: all_results.txt not found. Run mergeall first."
        echo ""
        return 1
    fi
    
    get_columns "$OUTPUT_DIR/all_results.txt"
    echo ""
    
    detect_groups "$OUTPUT_DIR/all_results.txt"
    
    # Save to file
    {
        echo "Annotation Group Summary"
        echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
        tail -n +2 "$OUTPUT_DIR/all_results.txt" | cut -f"$GROUPCOL" | sort | uniq -c | \
            awk '{printf "%-30s %8d genes\n", $2, $1}'
    } > "$OUTPUT_DIR/annotation_groups.txt"
    
    echo "  ✓ Saved to: annotation_groups.txt"
    echo ""
}

#==========================================
# Operation: Merge Chromosome Chunks
#==========================================

op_mergechrom() {
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "Merging chunked results by chromosome..."
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    
    if [ "$RESULT_TYPE" != "chunked" ]; then
        echo "  Skipping: Results are not chunked"
        echo ""
        return
    fi
    
    for chr in {1..22} X Y; do
        if ls "$OUTPUT_DIR"/chr${chr}_chunk*_results.txt 1> /dev/null 2>&1; then
            echo "  Processing chr${chr}..."
            
            FIRST_CHUNK=$(ls "$OUTPUT_DIR"/chr${chr}_chunk*_results.txt | sort -V | head -1)
            
            head -1 "$FIRST_CHUNK" > "$OUTPUT_DIR/chr${chr}_combined.txt"
            tail -n +2 -q "$OUTPUT_DIR"/chr${chr}_chunk*_results.txt >> "$OUTPUT_DIR/chr${chr}_combined.txt"
            
            COUNT=$(tail -n +2 "$OUTPUT_DIR/chr${chr}_combined.txt" | wc -l)
            echo "    Combined $COUNT genes"
        fi
    done
    
    RESULT_PATTERN="chr*_combined.txt"
    RESULT_TYPE="combined"
    echo "  ✓ Chromosome merging complete"
    echo ""
}

#==========================================
# Operation: Merge All Chromosomes
#==========================================

op_mergeall() {
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "Combining all chromosomes..."
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    
    if [ -f "$OUTPUT_DIR/all_results.txt" ]; then
        echo "  all_results.txt already exists"
        TOTAL_GENES=$(tail -n +2 "$OUTPUT_DIR/all_results.txt" | wc -l)
        echo "  Total genes: $TOTAL_GENES"
        
        # Show annotation groups
        get_columns "$OUTPUT_DIR/all_results.txt"
        echo ""
        detect_groups "$OUTPUT_DIR/all_results.txt"
        
        if [ -f "$OUTPUT_DIR/.gene_coords_lookup.txt" ]; then
            echo ""
            annotate_with_coords "$OUTPUT_DIR/all_results.txt" "$OUTPUT_DIR/all_results_annotated.txt"
            echo "  ✓ Saved annotated version: all_results_annotated.txt"
        fi
        echo ""
        return
    fi
    
    FIRST_FILE=$(ls "$OUTPUT_DIR"/$RESULT_PATTERN 2>/dev/null | sort -V | head -1)
    
    if [ -z "$FIRST_FILE" ]; then
        echo "ERROR: No result files found" >&2
        exit 1
    fi
    
    # Add CHR column if no coordinate file
    if [ ! -f "$OUTPUT_DIR/.gene_coords_lookup.txt" ]; then
        echo "  Adding chromosome information from filenames..."
        
        # Get header and add CHR column
        head -1 "$FIRST_FILE" | awk -F'\t' '{print "CHR\t" $0}' > "$OUTPUT_DIR/all_results.txt"
        
        # Process each chromosome file
        for file in "$OUTPUT_DIR"/$RESULT_PATTERN; do
            CHR=$(extract_chr_from_filename "$(basename "$file")")
            tail -n +2 "$file" | awk -v chr="$CHR" -F'\t' '{print chr "\t" $0}' >> "$OUTPUT_DIR/all_results.txt"
        done
    else
        # Regular merge without CHR column
        head -1 "$FIRST_FILE" > "$OUTPUT_DIR/all_results.txt"
        tail -n +2 -q "$OUTPUT_DIR"/$RESULT_PATTERN >> "$OUTPUT_DIR/all_results.txt"
    fi
    
    TOTAL_GENES=$(tail -n +2 "$OUTPUT_DIR/all_results.txt" | wc -l)
    echo "  Total genes tested: $TOTAL_GENES"
    echo "  ✓ Saved to: all_results.txt"
    echo ""
    
    # Show available annotation groups
    get_columns "$OUTPUT_DIR/all_results.txt"
    echo ""
    detect_groups "$OUTPUT_DIR/all_results.txt"
    
    # Annotate with coordinates if available
    if [ -f "$OUTPUT_DIR/.gene_coords_lookup.txt" ]; then
        echo ""
        annotate_with_coords "$OUTPUT_DIR/all_results.txt" "$OUTPUT_DIR/all_results_annotated.txt"
        echo "  ✓ Saved annotated version: all_results_annotated.txt"
    fi
    
    echo ""
}

#==========================================
# Operation: Extract Significant Results
#==========================================

op_findsig() {
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "Extracting significant results..."
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    
    if [ ! -f "$OUTPUT_DIR/all_results.txt" ]; then
        echo "  ERROR: all_results.txt not found. Run mergeall first."
        echo ""
        return 1
    fi
    
    get_columns "$OUTPUT_DIR/all_results.txt"
    echo ""
    
    # Apply filters if specified
    WORK_FILE="$OUTPUT_DIR/all_results.txt"
    if [ -n "$FILTER_GROUP" ] || [ -n "$FILTER_MAF" ]; then
        WORK_FILE="$OUTPUT_DIR/all_results_filtered.txt"
        apply_filters "$OUTPUT_DIR/all_results.txt" "$WORK_FILE" "$FILTER_GROUP" "$FILTER_MAF"
        echo ""
    fi
    
    # Genome-wide significant (p < 5e-8)
    head -1 "$WORK_FILE" > "$OUTPUT_DIR/genome_wide_sig.txt"
    tail -n +2 "$WORK_FILE" | awk -v pcol="$PCOL" -F'\t' '$pcol != "NA" && $pcol != "" && $pcol < 5e-8' >> "$OUTPUT_DIR/genome_wide_sig.txt"
    GWS_COUNT=$(tail -n +2 "$OUTPUT_DIR/genome_wide_sig.txt" | wc -l)
    echo "  Genome-wide (p < 5e-8):  $GWS_COUNT genes"
    
    # Suggestive (p < 1e-5)
    head -1 "$WORK_FILE" > "$OUTPUT_DIR/suggestive_sig.txt"
    tail -n +2 "$WORK_FILE" | awk -v pcol="$PCOL" -F'\t' '$pcol != "NA" && $pcol != "" && $pcol < 1e-5' >> "$OUTPUT_DIR/suggestive_sig.txt"
    SUG_COUNT=$(tail -n +2 "$OUTPUT_DIR/suggestive_sig.txt" | wc -l)
    echo "  Suggestive (p < 1e-5):   $SUG_COUNT genes"
    
    # Nominal (p < 0.05)
    head -1 "$WORK_FILE" > "$OUTPUT_DIR/nominal_sig.txt"
    tail -n +2 "$WORK_FILE" | awk -v pcol="$PCOL" -F'\t' '$pcol != "NA" && $pcol != "" && $pcol < 0.05' >> "$OUTPUT_DIR/nominal_sig.txt"
    NOM_COUNT=$(tail -n +2 "$OUTPUT_DIR/nominal_sig.txt" | wc -l)
    echo "  Nominal (p < 0.05):      $NOM_COUNT genes"
    
    # P < 0.01
    head -1 "$WORK_FILE" > "$OUTPUT_DIR/sig_p001.txt"
    tail -n +2 "$WORK_FILE" | awk -v pcol="$PCOL" -F'\t' '$pcol != "NA" && $pcol != "" && $pcol < 0.01' >> "$OUTPUT_DIR/sig_p001.txt"
    P001_COUNT=$(tail -n +2 "$OUTPUT_DIR/sig_p001.txt" | wc -l)
    echo "  P < 0.01:                $P001_COUNT genes"
    
    echo "  ✓ Significance filtering complete"
    echo ""
}

#==========================================
# Operation: Extract Genome-wide Significant
#==========================================

op_findgws() {
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "Extracting genome-wide significant results..."
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    
    if [ ! -f "$OUTPUT_DIR/all_results.txt" ]; then
        echo "  ERROR: all_results.txt not found. Run mergeall first."
        echo ""
        return 1
    fi
    
    get_columns "$OUTPUT_DIR/all_results.txt"
    echo ""
    
    # Apply filters if specified
    WORK_FILE="$OUTPUT_DIR/all_results.txt"
    if [ -n "$FILTER_GROUP" ] || [ -n "$FILTER_MAF" ]; then
        WORK_FILE="$OUTPUT_DIR/all_results_filtered.txt"
        apply_filters "$OUTPUT_DIR/all_results.txt" "$WORK_FILE" "$FILTER_GROUP" "$FILTER_MAF"
        echo ""
    fi
    
    head -1 "$WORK_FILE" > "$OUTPUT_DIR/genome_wide_sig.txt"
    tail -n +2 "$WORK_FILE" | awk -v pcol="$PCOL" -F'\t' '$pcol != "NA" && $pcol != "" && $pcol < 5e-8' >> "$OUTPUT_DIR/genome_wide_sig.txt"
    GWS_COUNT=$(tail -n +2 "$OUTPUT_DIR/genome_wide_sig.txt" | wc -l)
    echo "  Genome-wide significant: $GWS_COUNT genes"
    echo "  ✓ Saved to: genome_wide_sig.txt"
    echo ""
}

#==========================================
# Operation: Extract Suggestive Results
#==========================================

op_findsug() {
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "Extracting suggestive results..."
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    
    if [ ! -f "$OUTPUT_DIR/all_results.txt" ]; then
        echo "  ERROR: all_results.txt not found. Run mergeall first."
        echo ""
        return 1
    fi
    
    get_columns "$OUTPUT_DIR/all_results.txt"
    echo ""
    
    # Apply filters if specified
    WORK_FILE="$OUTPUT_DIR/all_results.txt"
    if [ -n "$FILTER_GROUP" ] || [ -n "$FILTER_MAF" ]; then
        WORK_FILE="$OUTPUT_DIR/all_results_filtered.txt"
        apply_filters "$OUTPUT_DIR/all_results.txt" "$WORK_FILE" "$FILTER_GROUP" "$FILTER_MAF"
        echo ""
    fi
    
    head -1 "$WORK_FILE" > "$OUTPUT_DIR/suggestive_sig.txt"
    tail -n +2 "$WORK_FILE" | awk -v pcol="$PCOL" -F'\t' '$pcol != "NA" && $pcol != "" && $pcol < 1e-5' >> "$OUTPUT_DIR/suggestive_sig.txt"
    SUG_COUNT=$(tail -n +2 "$OUTPUT_DIR/suggestive_sig.txt" | wc -l)
    echo "  Suggestive results: $SUG_COUNT genes"
    echo "  ✓ Saved to: suggestive_sig.txt"
    echo ""
}

#==========================================
# Operation: Extract Nominal Significant
#==========================================

op_findnom() {
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "Extracting nominal significant results..."
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    
    if [ ! -f "$OUTPUT_DIR/all_results.txt" ]; then
        echo "  ERROR: all_results.txt not found. Run mergeall first."
        echo ""
        return 1
    fi
    
    get_columns "$OUTPUT_DIR/all_results.txt"
    echo ""
    
    # Apply filters if specified
    WORK_FILE="$OUTPUT_DIR/all_results.txt"
    if [ -n "$FILTER_GROUP" ] || [ -n "$FILTER_MAF" ]; then
        WORK_FILE="$OUTPUT_DIR/all_results_filtered.txt"
        apply_filters "$OUTPUT_DIR/all_results.txt" "$WORK_FILE" "$FILTER_GROUP" "$FILTER_MAF"
        echo ""
    fi
    
    head -1 "$WORK_FILE" > "$OUTPUT_DIR/nominal_sig.txt"
    tail -n +2 "$WORK_FILE" | awk -v pcol="$PCOL" -F'\t' '$pcol != "NA" && $pcol != "" && $pcol < 0.05' >> "$OUTPUT_DIR/nominal_sig.txt"
    NOM_COUNT=$(tail -n +2 "$OUTPUT_DIR/nominal_sig.txt" | wc -l)
    echo "  Nominal significant: $NOM_COUNT genes"
    echo "  ✓ Saved to: nominal_sig.txt"
    echo ""
}

#==========================================
# Operation: Extract Top Genes
#==========================================

op_topgenes() {
    local n=$1
    local outfile="$OUTPUT_DIR/top${n}_genes.txt"
    
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "Extracting top $n genes..."
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    
    if [ ! -f "$OUTPUT_DIR/all_results.txt" ]; then
        echo "  ERROR: all_results.txt not found. Run mergeall first."
        echo ""
        return 1
    fi
    
    get_columns "$OUTPUT_DIR/all_results.txt"
    echo ""
    
    # Apply filters if specified
    WORK_FILE="$OUTPUT_DIR/all_results.txt"
    if [ -n "$FILTER_GROUP" ] || [ -n "$FILTER_MAF" ]; then
        WORK_FILE="$OUTPUT_DIR/all_results_filtered.txt"
        apply_filters "$OUTPUT_DIR/all_results.txt" "$WORK_FILE" "$FILTER_GROUP" "$FILTER_MAF"
        echo ""
    fi
    
    head -1 "$WORK_FILE" > "$outfile"
    tail -n +2 "$WORK_FILE" | \
        awk -v pcol="$PCOL" -F'\t' '$pcol != "NA" && $pcol != ""' | \
        sort -t$'\t' -k"$PCOL","$PCOL"g | \
        sed -n "1,${n}p" >> "$outfile"
    
    echo "  ✓ Saved to: top${n}_genes.txt"
    
    # Display top 20 or fewer
    DISPLAY_N=$n
    [ $n -gt 20 ] && DISPLAY_N=20
    
    echo ""
    echo "  Top $DISPLAY_N genes:"
    echo "  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    
    if command -v column &> /dev/null; then
        # Avoid SIGPIPE under "set -o pipefail" by not piping formatted output to a second head.
        head -$((DISPLAY_N + 1)) "$outfile" | column -t -s$'\t'
    else
        head -$((DISPLAY_N + 1)) "$outfile"
    fi
    
    echo ""
}

#==========================================
# Operation: Chromosome Summary
#==========================================

op_chromsum() {
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "Generating per-chromosome summary..."
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    
    {
        printf "%-5s %-12s %-12s %-12s %-12s %-12s\n" \
            "Chr" "Total_Genes" "GWS_p<5e-8" "Sug_p<1e-5" "Nom_p<0.05" "P<0.01"
        echo "------------------------------------------------------------------------"
        
        for chr in {1..22} X Y; do
            CHR_FILE=""
            if [ -f "$OUTPUT_DIR/chr${chr}_combined_results.txt" ]; then
                CHR_FILE="$OUTPUT_DIR/chr${chr}_combined_results.txt"
            elif [ -f "$OUTPUT_DIR/chr${chr}_combined.txt" ]; then
                CHR_FILE="$OUTPUT_DIR/chr${chr}_combined.txt"
            elif [ -f "$OUTPUT_DIR/chr${chr}_all_results.txt" ]; then
                CHR_FILE="$OUTPUT_DIR/chr${chr}_all_results.txt"
            fi
            
            if [ -n "$CHR_FILE" ] && [ -f "$CHR_FILE" ]; then
                get_columns "$CHR_FILE" > /dev/null 2>&1
                
                TOTAL=$(tail -n +2 "$CHR_FILE" | wc -l)
                GWS=$(tail -n +2 "$CHR_FILE" | awk -v pcol="$PCOL" -F'\t' '$pcol != "NA" && $pcol != "" && $pcol < 5e-8' | wc -l)
                SUG=$(tail -n +2 "$CHR_FILE" | awk -v pcol="$PCOL" -F'\t' '$pcol != "NA" && $pcol != "" && $pcol < 1e-5' | wc -l)
                NOM=$(tail -n +2 "$CHR_FILE" | awk -v pcol="$PCOL" -F'\t' '$pcol != "NA" && $pcol != "" && $pcol < 0.05' | wc -l)
                P001=$(tail -n +2 "$CHR_FILE" | awk -v pcol="$PCOL" -F'\t' '$pcol != "NA" && $pcol != "" && $pcol < 0.01' | wc -l)
                
                printf "%-5s %-12s %-12s %-12s %-12s %-12s\n" \
                    "$chr" "$TOTAL" "$GWS" "$SUG" "$NOM" "$P001"
            fi
        done
    } > "$OUTPUT_DIR/chromosome_summary.txt"
    
    cat "$OUTPUT_DIR/chromosome_summary.txt"
    echo ""
    echo "  ✓ Saved to: chromosome_summary.txt"
    echo ""
}

#==========================================
# Operation: Annotation Group Summary
#==========================================

op_groupsum() {
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "Generating per-annotation-group summary..."
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    
    if [ ! -f "$OUTPUT_DIR/all_results.txt" ]; then
        echo "  ERROR: all_results.txt not found. Run mergeall first."
        echo ""
        return 1
    fi
    
    get_columns "$OUTPUT_DIR/all_results.txt"
    echo ""
    
    {
        printf "%-30s %-12s %-12s %-12s %-12s %-12s\n" \
            "Annotation_Group" "Total_Genes" "GWS_p<5e-8" "Sug_p<1e-5" "Nom_p<0.05" "P<0.01"
        echo "----------------------------------------------------------------------------------------"
        
        # Get unique groups
        tail -n +2 "$OUTPUT_DIR/all_results.txt" | cut -f"$GROUPCOL" | sort -u | while read -r group; do
            TOTAL=$(tail -n +2 "$OUTPUT_DIR/all_results.txt" | awk -v gcol="$GROUPCOL" -v grp="$group" -F'\t' '$gcol == grp' | wc -l)
            GWS=$(tail -n +2 "$OUTPUT_DIR/all_results.txt" | awk -v gcol="$GROUPCOL" -v pcol="$PCOL" -v grp="$group" -F'\t' '$gcol == grp && $pcol != "NA" && $pcol != "" && $pcol < 5e-8' | wc -l)
            SUG=$(tail -n +2 "$OUTPUT_DIR/all_results.txt" | awk -v gcol="$GROUPCOL" -v pcol="$PCOL" -v grp="$group" -F'\t' '$gcol == grp && $pcol != "NA" && $pcol != "" && $pcol < 1e-5' | wc -l)
            NOM=$(tail -n +2 "$OUTPUT_DIR/all_results.txt" | awk -v gcol="$GROUPCOL" -v pcol="$PCOL" -v grp="$group" -F'\t' '$gcol == grp && $pcol != "NA" && $pcol != "" && $pcol < 0.05' | wc -l)
            P001=$(tail -n +2 "$OUTPUT_DIR/all_results.txt" | awk -v gcol="$GROUPCOL" -v pcol="$PCOL" -v grp="$group" -F'\t' '$gcol == grp && $pcol != "NA" && $pcol != "" && $pcol < 0.01' | wc -l)
            
            printf "%-30s %-12s %-12s %-12s %-12s %-12s\n" \
                "$group" "$TOTAL" "$GWS" "$SUG" "$NOM" "$P001"
        done
    } > "$OUTPUT_DIR/annotation_group_summary.txt"
    
    cat "$OUTPUT_DIR/annotation_group_summary.txt"
    echo ""
    echo "  ✓ Saved to: annotation_group_summary.txt"
    echo ""
}

#==========================================
# Operation: Full Summary Report
#==========================================

op_fullsum() {
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "Generating full summary report..."
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    
    if [ ! -f "$OUTPUT_DIR/all_results.txt" ]; then
        echo "  ERROR: all_results.txt not found. Run mergeall first."
        echo ""
        return 1
    fi
    
    get_columns "$OUTPUT_DIR/all_results.txt"
    echo ""
    
    TOTAL_GENES=$(tail -n +2 "$OUTPUT_DIR/all_results.txt" | wc -l)
    
    # Count significant results
    GWS_COUNT=$(tail -n +2 "$OUTPUT_DIR/all_results.txt" | awk -v pcol="$PCOL" -F'\t' '$pcol != "NA" && $pcol != "" && $pcol < 5e-8' | wc -l)
    SUG_COUNT=$(tail -n +2 "$OUTPUT_DIR/all_results.txt" | awk -v pcol="$PCOL" -F'\t' '$pcol != "NA" && $pcol != "" && $pcol < 1e-5' | wc -l)
    NOM_COUNT=$(tail -n +2 "$OUTPUT_DIR/all_results.txt" | awk -v pcol="$PCOL" -F'\t' '$pcol != "NA" && $pcol != "" && $pcol < 0.05' | wc -l)
    P001_COUNT=$(tail -n +2 "$OUTPUT_DIR/all_results.txt" | awk -v pcol="$PCOL" -F'\t' '$pcol != "NA" && $pcol != "" && $pcol < 0.01' | wc -l)
    
    {
        echo "=========================================="
        echo "SAIGE-GENE Results Summary"
        echo "=========================================="
        echo "Analysis Date: $(date)"
        echo "Results Directory: $OUTPUT_DIR"
        echo ""
        
        if [ -n "$FILTER_GROUP" ]; then
            echo "Filter Applied: Annotation Group = '$FILTER_GROUP'"
        fi
        if [ -n "$FILTER_MAF" ]; then
            echo "Filter Applied: max_MAF <= $FILTER_MAF"
        fi
        if [ -n "$FILTER_GROUP" ] || [ -n "$FILTER_MAF" ]; then
            echo ""
        fi
        
        echo "Total Genes Tested: $TOTAL_GENES"
        echo ""
        echo "Significance Thresholds:"
        echo "----------------------------------------"
        printf "  %-25s %8s %8s\n" "Threshold" "Count" "Percent"
        echo "----------------------------------------"
        printf "  %-25s %8s %7.2f%%\n" "Genome-wide (p < 5e-8)" "$GWS_COUNT" "$(awk -v g="$GWS_COUNT" -v t="$TOTAL_GENES" 'BEGIN {printf "%.2f", (g/t)*100}')"
        printf "  %-25s %8s %7.2f%%\n" "Suggestive (p < 1e-5)" "$SUG_COUNT" "$(awk -v g="$SUG_COUNT" -v t="$TOTAL_GENES" 'BEGIN {printf "%.2f", (g/t)*100}')"
        printf "  %-25s %8s %7.2f%%\n" "P < 0.01" "$P001_COUNT" "$(awk -v g="$P001_COUNT" -v t="$TOTAL_GENES" 'BEGIN {printf "%.2f", (g/t)*100}')"
        printf "  %-25s %8s %7.2f%%\n" "Nominal (p < 0.05)" "$NOM_COUNT" "$(awk -v g="$NOM_COUNT" -v t="$TOTAL_GENES" 'BEGIN {printf "%.2f", (g/t)*100}')"
        echo ""
        
        echo "Annotation Groups:"
        echo "----------------------------------------"
        tail -n +2 "$OUTPUT_DIR/all_results.txt" | cut -f"$GROUPCOL" | sort | uniq -c | \
            awk '{printf "  %-30s %8d genes\n", $2, $1}'
        echo ""
        
        echo "Output Files Generated:"
        echo "----------------------------------------"
        
        if [ -f "$OUTPUT_DIR/all_results.txt" ]; then
            echo "  ✓ all_results.txt"
        fi
        if [ -f "$OUTPUT_DIR/all_results_annotated.txt" ]; then
            echo "  ✓ all_results_annotated.txt (with coordinates)"
        fi
        if [ -f "$OUTPUT_DIR/genome_wide_sig.txt" ]; then
            echo "  ✓ genome_wide_sig.txt"
        fi
        if [ -f "$OUTPUT_DIR/suggestive_sig.txt" ]; then
            echo "  ✓ suggestive_sig.txt"
        fi
        if [ -f "$OUTPUT_DIR/nominal_sig.txt" ]; then
            echo "  ✓ nominal_sig.txt"
        fi
        if [ -f "$OUTPUT_DIR/sig_p001.txt" ]; then
            echo "  ✓ sig_p001.txt"
        fi
        if [ -f "$OUTPUT_DIR/top50_genes.txt" ]; then
            echo "  ✓ top50_genes.txt"
        fi
        if [ -f "$OUTPUT_DIR/chromosome_summary.txt" ]; then
            echo "  ✓ chromosome_summary.txt"
        fi
        if [ -f "$OUTPUT_DIR/annotation_group_summary.txt" ]; then
            echo "  ✓ annotation_group_summary.txt"
        fi
        if [ -f "$OUTPUT_DIR/qq_plot_data.txt" ]; then
            echo "  ✓ qq_plot_data.txt"
        fi
        if [ -f "$OUTPUT_DIR/manhattan_plot_data.txt" ]; then
            echo "  ✓ manhattan_plot_data.txt"
        fi
        echo ""
        echo "=========================================="
    } | tee "$OUTPUT_DIR/analysis_summary.txt"
    
    echo ""
    echo "  ✓ Summary saved to: analysis_summary.txt"
    echo ""
}

#==========================================
# Operation: QQ Plot Data
#==========================================

op_qqdata() {
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "Generating QQ plot data..."
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    
    if [ ! -f "$OUTPUT_DIR/all_results.txt" ]; then
        echo "  ERROR: all_results.txt not found. Run mergeall first."
        echo ""
        return 1
    fi
    
    get_columns "$OUTPUT_DIR/all_results.txt"
    echo ""
    
    # Apply filters if specified
    WORK_FILE="$OUTPUT_DIR/all_results.txt"
    if [ -n "$FILTER_GROUP" ] || [ -n "$FILTER_MAF" ]; then
        WORK_FILE="$OUTPUT_DIR/all_results_filtered.txt"
        apply_filters "$OUTPUT_DIR/all_results.txt" "$WORK_FILE" "$FILTER_GROUP" "$FILTER_MAF"
        echo ""
    fi
    
    TOTAL_GENES=$(tail -n +2 "$WORK_FILE" | awk -v pcol="$PCOL" -F'\t' '$pcol != "NA" && $pcol != "" && $pcol > 0' | wc -l)
    
    tail -n +2 "$WORK_FILE" | \
        awk -v pcol="$PCOL" -F'\t' '$pcol != "NA" && $pcol != "" && $pcol > 0 {print $pcol}' | \
        sort -g | \
        awk -v total="$TOTAL_GENES" 'BEGIN {
            print "Expected_log10P\tObserved_log10P"
        } {
            obs = -log($1)/log(10)
            expected = -log((NR-0.5)/total)/log(10)
            print expected "\t" obs
        }' > "$OUTPUT_DIR/qq_plot_data.txt"
    
    POINTS=$(tail -n +2 "$OUTPUT_DIR/qq_plot_data.txt" | wc -l)
    echo "  Generated $POINTS data points"
    echo "  ✓ Saved to: qq_plot_data.txt"
    echo ""
    
    # Calculate genomic inflation factor (lambda)
    LAMBDA=$(tail -n +2 "$WORK_FILE" | \
        awk -v pcol="$PCOL" -F'\t' '$pcol != "NA" && $pcol != "" && $pcol > 0 {print $pcol}' | \
        sort -g | \
        awk 'BEGIN{n=0} {a[n++]=$1} END{
            if (n > 0) {
                median_p = a[int(n/2)]
                chi2 = -2 * log(median_p)
                lambda = chi2 / 0.4549364
                printf "%.4f", lambda
            }
        }')
    
    if [ -n "$LAMBDA" ]; then
        echo "  Genomic Inflation Factor (λ): $LAMBDA"
        echo ""
    fi
}

#==========================================
# Operation: Manhattan Plot Data
#==========================================

op_mandata() {
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "Generating Manhattan plot data..."
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    
    if [ ! -f "$OUTPUT_DIR/all_results.txt" ]; then
        echo "  ERROR: all_results.txt not found. Run mergeall first."
        echo ""
        return 1
    fi
    
    get_columns "$OUTPUT_DIR/all_results.txt"
    echo ""
    
    # Apply filters if specified
    WORK_FILE="$OUTPUT_DIR/all_results.txt"
    if [ -n "$FILTER_GROUP" ] || [ -n "$FILTER_MAF" ]; then
        WORK_FILE="$OUTPUT_DIR/all_results_filtered.txt"
        apply_filters "$OUTPUT_DIR/all_results.txt" "$WORK_FILE" "$FILTER_GROUP" "$FILTER_MAF"
        echo ""
    fi
    
    # Check if we have annotated file with coordinates
    if [ -f "$OUTPUT_DIR/all_results_annotated.txt" ] && [ -f "$OUTPUT_DIR/.gene_coords_lookup.txt" ]; then
        echo "  Using gene coordinates from annotation file"
        
        # Use annotated file
        ANNO_WORK_FILE="$OUTPUT_DIR/all_results_annotated.txt"
        
        # Apply filters to annotated file if needed
        if [ -n "$FILTER_GROUP" ] || [ -n "$FILTER_MAF" ]; then
            get_columns "$ANNO_WORK_FILE" > /dev/null 2>&1
            
            head -1 "$ANNO_WORK_FILE" > "$OUTPUT_DIR/all_results_annotated_filtered.txt"
            
            # Adjust column numbers for annotated file (CHR and POS are prepended)
            local adj_groupcol=$((GROUPCOL+2))
            local adj_pcol=$((PCOL+2))
            local adj_mafcol=""
            [ -n "$MAFCOL" ] && adj_mafcol=$((MAFCOL+2))
            
            # Build filter
            local filter_cmd='NR > 1'
            if [ -n "$FILTER_GROUP" ]; then
                filter_cmd="$filter_cmd && \$$adj_groupcol == \"$FILTER_GROUP\""
            fi
            if [ -n "$FILTER_MAF" ] && [ -n "$adj_mafcol" ]; then
                filter_cmd="$filter_cmd && \$$adj_mafcol != \"NA\" && \$$adj_mafcol <= $FILTER_MAF"
            fi
            
            tail -n +2 "$ANNO_WORK_FILE" | awk -F'\t' "$filter_cmd" >> "$OUTPUT_DIR/all_results_annotated_filtered.txt"
            ANNO_WORK_FILE="$OUTPUT_DIR/all_results_annotated_filtered.txt"
        fi
        
        # Create header (Max_MAF for shape coding on Manhattan plots)
        echo -e "Gene\tCHR\tPOS\tPvalue\tNegLog10P\tGroup\tMax_MAF" > "$OUTPUT_DIR/manhattan_plot_data.txt"

        _maf_a=0
        [ -n "${MAFCOL:-}" ] && _maf_a=$((MAFCOL + 2))

        tail -n +2 "$ANNO_WORK_FILE" | \
            awk -v rcol="$((REGIONCOL+2))" -v gcol="$((GROUPCOL+2))" -v pcol="$((PCOL+2))" -v mafc="$_maf_a" -F'\t' '
            $pcol != "NA" && $pcol != "" && $pcol > 0 {
                chr = $1
                pos = $2
                gene = $rcol
                group = $gcol
                pval = $pcol
                maf = (mafc > 0 && mafc <= NF) ? $mafc : "NA"
                
                # Extract just gene name from region
                if (gene ~ /:/) {
                    split(gene, arr, /:/)
                    gene = arr[length(arr)]
                    if (gene ~ /^[0-9]+$/) gene = $rcol
                }
                
                # Calculate -log10(p)
                log10p = -log(pval)/log(10)
                
                print gene "\t" chr "\t" pos "\t" pval "\t" log10p "\t" group "\t" maf
            }' >> "$OUTPUT_DIR/manhattan_plot_data.txt"
        
    else
        echo "  No gene coordinates available, using chromosome from merged data"
        
        # Create header (Max_MAF encodes SAIGE mask for point shapes)
        echo -e "Gene\tCHR\tPOS\tPvalue\tNegLog10P\tGroup\tMax_MAF" > "$OUTPUT_DIR/manhattan_plot_data.txt"
        
        # Check if CHR column exists in all_results.txt
        HAS_CHR_COL=false
        if head -1 "$WORK_FILE" | grep -q "^CHR"; then
            HAS_CHR_COL=true
            echo "  Using CHR column from merged results"
            # Adjust column indices
            CHRCOL=1
            adj_regioncol=$((REGIONCOL+1))
            adj_groupcol=$((GROUPCOL+1))
            adj_pcol=$((PCOL+1))
            if [ -n "${MAFCOL:-}" ]; then adj_mafcol=$((MAFCOL+1)); else adj_mafcol=0; fi
        else
            echo "  Warning: No CHR information available, using gene index for position"
            adj_regioncol=$REGIONCOL
            adj_groupcol=$GROUPCOL
            adj_pcol=$PCOL
            if [ -n "${MAFCOL:-}" ]; then adj_mafcol=$MAFCOL; else adj_mafcol=0; fi
        fi
        
        if [ "$HAS_CHR_COL" = true ]; then
            # Use CHR column from data
            tail -n +2 "$WORK_FILE" | \
                awk -v chrcol="$CHRCOL" -v rcol="$adj_regioncol" -v gcol="$adj_groupcol" -v pcol="$adj_pcol" -v mafc="$adj_mafcol" -F'\t' '
                $pcol != "NA" && $pcol != "" && $pcol > 0 {
                    chr = $chrcol
                    gene = $rcol
                    group = $gcol
                    pval = $pcol
                    maf = (mafc > 0 && mafc <= NF) ? $mafc : "NA"
                    
                    # Extract gene name from region
                    if (gene ~ /:/) {
                        split(gene, arr, /:/)
                        gene = arr[length(arr)]
                        if (gene ~ /^[0-9]+$/) gene = $rcol
                    }
                    
                    # Calculate -log10(p)
                    log10p = -log(pval)/log(10)
                    
                    # Use gene index as position
                    pos = NR
                    
                    print gene "\t" chr "\t" pos "\t" pval "\t" log10p "\t" group "\t" maf
                }' >> "$OUTPUT_DIR/manhattan_plot_data.txt"
        else
            # No CHR info, extract from region or use NA
            tail -n +2 "$WORK_FILE" | \
                awk -v rcol="$adj_regioncol" -v gcol="$adj_groupcol" -v pcol="$adj_pcol" -v mafc="$adj_mafcol" -F'\t' '
                $pcol != "NA" && $pcol != "" && $pcol > 0 {
                    gene = $rcol
                    group = $gcol
                    pval = $pcol
                    maf = (mafc > 0 && mafc <= NF) ? $mafc : "NA"
                    
                    # Try to extract chromosome from region
                    chr = "NA"
                    if ($rcol ~ /^chr/) {
                        split($rcol, arr, /:/)
                        chr = arr[1]
                        gsub("chr", "", chr)
                    }
                    
                    # Extract gene name from region
                    if (gene ~ /:/) {
                        split(gene, arr, /:/)
                        gene = arr[length(arr)]
                        if (gene ~ /^[0-9]+$/) gene = $rcol
                    }
                    
                    # Calculate -log10(p)
                    log10p = -log(pval)/log(10)
                    
                    # Use gene index as position
                    pos = NR
                    
                    print gene "\t" chr "\t" pos "\t" pval "\t" log10p "\t" group "\t" maf
                }' >> "$OUTPUT_DIR/manhattan_plot_data.txt"
        fi
    fi
    
    POINTS=$(tail -n +2 "$OUTPUT_DIR/manhattan_plot_data.txt" | wc -l)
    echo "  Generated $POINTS data points"
    echo "  ✓ Saved to: manhattan_plot_data.txt"
    echo "  Columns: Gene, CHR, POS, Pvalue, NegLog10P, Group, Max_MAF"
    echo ""
}



#==========================================
# Operation: Generate Plots and Tables
#==========================================

op_makeplots() {
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "Generating plots and summary tables..."
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

    if ! command -v Rscript >/dev/null 2>&1; then
        echo "  WARNING: Rscript not found. Skipping PNG/table generation."
        echo "  Install R and rerun with operation: makeplots"
        echo ""
        return 0
    fi

    if [ ! -f "$OUTPUT_DIR/all_results.txt" ]; then
        echo "  ERROR: all_results.txt not found. Run mergeall before makeplots."
        echo ""
        return 1
    fi

    get_columns "$OUTPUT_DIR/all_results.txt" >/dev/null 2>&1 || {
        echo "  ERROR: Could not read column headers in all_results.txt"
        return 1
    }

    # Unique genes in merged results (Region column), for Manhattan Bonferroni thresholds
    N_GENES_MERGED=$(awk -F'\t' -v rc="$REGIONCOL" '
      NR > 1 && rc >= 1 && rc <= NF {
        r = $rc
        if (r == "" || r == "NA") next
        if (r ~ /:/) {
          split(r, arr, /:/)
          r = arr[length(arr)]
          if (r ~ /^[0-9]+$/) r = $rc
        }
        seen[r] = 1
      }
      END { n = 0; for (k in seen) n++; print n + 0 }
    ' "$OUTPUT_DIR/all_results.txt")
    [ -z "$N_GENES_MERGED" ] && N_GENES_MERGED=0
    echo "  Unique genes (merged all_results.txt, by Region): $N_GENES_MERGED"

    MANHATTAN_LABEL_TOP_N="${STEP9_MANHATTAN_LABEL_TOP_N:-10}"
    if [ "$INTERACTIVE" = true ]; then
        echo "  Manhattan plot: by default the top ${MANHATTAN_LABEL_TOP_N} lowest-P genes are labeled."
        read -r -p "  Label more genes than that? [y/N] (Enter=no): " ML_MORE
        case "${ML_MORE:-N}" in
            y|Y|yes|YES)
                read -r -p "  Total number of genes to label by name [25]: " MLN
                if [ -z "${MLN:-}" ]; then
                    MANHATTAN_LABEL_TOP_N=25
                elif [ "$MLN" -eq "$MLN" ] 2>/dev/null && [ "$MLN" -gt 0 ]; then
                    MANHATTAN_LABEL_TOP_N="$MLN"
                else
                    echo "  (invalid number; using 25)"
                    MANHATTAN_LABEL_TOP_N=25
                fi
                ;;
            *) ;;
        esac
    fi
    export STEP9_MANHATTAN_LABEL_TOP_N="$MANHATTAN_LABEL_TOP_N"
    echo "  Manhattan gene labels (top by P): $MANHATTAN_LABEL_TOP_N"

    MANHATTAN_FDR_ALPHA="${STEP9_MANHATTAN_FDR_ALPHA:-0.1}"
    if [ "$INTERACTIVE" = true ]; then
        echo "  Manhattan: two PNGs (Bonferroni thresholds + BH-FDR threshold)."
        read -r -p "  BH-FDR significance alpha [${MANHATTAN_FDR_ALPHA}]: " FA_IN
        if [ -n "${FA_IN:-}" ]; then
            if awk -v a="$FA_IN" 'BEGIN { exit !(a > 0 && a <= 1) }'; then
                MANHATTAN_FDR_ALPHA="$FA_IN"
            else
                echo "  (alpha must be in (0,1]; keeping ${MANHATTAN_FDR_ALPHA})"
            fi
        fi
    fi
    export STEP9_MANHATTAN_FDR_ALPHA="$MANHATTAN_FDR_ALPHA"
    echo "  Manhattan BH-FDR alpha (second PNG): $MANHATTAN_FDR_ALPHA"

    if [ "$PLOT_DO_UNFILTERED" != "1" ]; then
        echo "  Note: Unfiltered plot mode was off; enabling single Manhattan/QQ (all groups, all max_MAF)."
        PLOT_DO_UNFILTERED=1
    fi

    local r_script="$OUTPUT_DIR/.step9_generate_plots.R"
    cat > "$r_script" << 'RSCRIPT'
args <- commandArgs(trailingOnly = TRUE)
legacy <- length(args) < 4
n_genes_merged <- NA_integer_
if (!legacy && length(args) >= 5) {
  n_genes_merged <- suppressWarnings(as.integer(args[5]))
}
if (!legacy) {
  out_dir <- args[1]
  tag <- args[2]
  qq_fn <- args[3]
  man_fn <- args[4]
} else {
  out_dir <- args[1]
  tag <- "default"
  qq_fn <- "qq_plot_data.txt"
  man_fn <- "manhattan_plot_data.txt"
}
setwd(out_dir)

safe_dev_off <- function() {
  while (!is.null(dev.list())) dev.off()
}

plot_qq <- function(inpath, pngpath, mtitle) {
  if (!file.exists(inpath)) {
    return(invisible(FALSE))
  }
  qq_data <- tryCatch(
    read.table(inpath, header = TRUE, sep = "\t", stringsAsFactors = FALSE),
    error = function(e) NULL
  )
  if (is.null(qq_data) || nrow(qq_data) == 0) {
    return(invisible(FALSE))
  }
  png(pngpath, width = 1000, height = 1000, res = 150)
  graphics::par(bg = "white", col.axis = "gray25", col.lab = "gray20")
  plot(
    qq_data$Expected_log10P, qq_data$Observed_log10P,
    xlab = expression(-log[10](Expected~P)),
    ylab = expression(-log[10](Observed~P)),
    main = mtitle,
    pch = 20, cex = 0.8, col = rgb(0, 0, 0, 0.5)
  )
  abline(0, 1, col = "red", lwd = 2)
  safe_dev_off()
  invisible(TRUE)
}

plot_man <- function(
  inpath,
  pngpath,
  mtitle,
  tag,
  n_genes_merged = NA_integer_,
  mt_mode = NULL,
  write_group_stats = TRUE
) {
  if (!file.exists(inpath)) {
    return(invisible(FALSE))
  }
  man <- tryCatch(
    read.table(inpath, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) NULL
  )
  if (is.null(man)) {
    return(invisible(FALSE))
  }
  man <- man[!is.na(man$NegLog10P) & !is.na(man$Pvalue) & man$Pvalue > 0, ]
  if (nrow(man) == 0) {
    return(invisible(FALSE))
  }

  if (!("Max_MAF" %in% names(man))) {
    man$Max_MAF <- NA_character_
  }
  man$Max_MAF_lab <- ifelse(
    is.na(man$Max_MAF) | man$Max_MAF == "",
    "NA",
    as.character(man$Max_MAF)
  )

  if (!is.na(n_genes_merged) && !is.na(suppressWarnings(as.integer(n_genes_merged))) &&
      as.integer(n_genes_merged) > 0L) {
    n_genes <- as.integer(n_genes_merged)
  } else {
    n_genes <- length(unique(man$Gene))
  }
  maf_tests <- suppressWarnings(as.numeric(Sys.getenv("STEP9_BONFERRONI_MAF_TESTS", "3")))
  if (!is.finite(maf_tests) || maf_tests < 1) {
    maf_tests <- 3
  }
  p_bonf_gene <- 0.05 / max(n_genes, 1L)
  p_bonf_gene_maf <- 0.05 / max(as.numeric(n_genes) * maf_tests, 1)
  y_bonf_gene <- -log10(p_bonf_gene)
  y_bonf_maf <- -log10(p_bonf_gene_maf)

  label_n <- suppressWarnings(as.integer(Sys.getenv("STEP9_MANHATTAN_LABEL_TOP_N", "10")))
  if (is.na(label_n) || label_n < 1L) {
    label_n <- 10L
  }
  bonf_line1 <- sprintf("Orange: 0.05 / %d genes  ~  -log10 = %.3f", n_genes, y_bonf_gene)
  bonf_line2 <- sprintf(
    "Green: 0.05 / (%d x %.0f)  ~  -log10 = %.3f",
    n_genes, maf_tests, y_bonf_maf
  )

  chr_txt <- as.character(man$CHR)
  man$CHR_num <- suppressWarnings(
    as.numeric(ifelse(chr_txt == "X", "23", ifelse(chr_txt == "Y", "24", chr_txt)))
  )
  ok_chr <- !is.na(man$CHR_num) & !(toupper(chr_txt) %in% c("NA", "")) & chr_txt != "NA"

  chr_vlines <- numeric(0)
  chr_axis_labels <- character(0)

  if (sum(ok_chr) == 0) {
    man$CHR_num <- 1L
    man$BP_cum <- seq_len(nrow(man))
    chr_levels <- 1
    midpoints <- median(man$BP_cum)
    chr_axis_labels <- "1"
  } else {
    man <- man[ok_chr, ]
    if (nrow(man) == 0) {
      return(invisible(FALSE))
    }
    man <- man[order(man$CHR_num, man$POS), ]
    chr_levels <- sort(unique(man$CHR_num))
    nchr <- length(chr_levels)
    chr_span <- vapply(chr_levels, function(ch) {
      mp <- max(man$POS[man$CHR_num == ch], na.rm = TRUE)
      if (!is.finite(mp) || mp <= 0) {
        mp <- max(sum(man$CHR_num == ch), 1)
      }
      mp
    }, numeric(1))
    names(chr_span) <- as.character(chr_levels)
    total_span <- sum(chr_span)
    gap <- if (is.finite(total_span) && total_span > 0) {
      max(total_span * 0.02 / max(nchr, 1L), stats::median(chr_span) * 0.015, 1)
    } else {
      1
    }
    starts <- numeric(nchr)
    names(starts) <- as.character(chr_levels)
    offset <- 0
    for (i in seq_along(chr_levels)) {
      ch <- chr_levels[i]
      starts[as.character(ch)] <- offset
      offset <- offset + chr_span[as.character(ch)] + gap
    }
    man$BP_cum <- man$POS + as.numeric(starts[as.character(man$CHR_num)])
    midpoints <- vapply(chr_levels, function(ch) {
      idx <- man$CHR_num == ch
      (min(man$BP_cum[idx]) + max(man$BP_cum[idx])) / 2
    }, numeric(1))
    if (nchr > 1L) {
      chr_vlines <- starts[as.character(chr_levels[-1L])]
    }
    chr_axis_labels <- vapply(chr_levels, function(ch) {
      if (ch == 23L) {
        "X"
      } else if (ch == 24L) {
        "Y"
      } else {
        as.character(ch)
      }
    }, character(1))
  }

  if (!is.null(mt_mode) && nzchar(as.character(mt_mode))) {
    use_fdr <- tolower(trimws(as.character(mt_mode))) %in% c("fdr", "bh")
  } else {
    mt_raw <- tolower(trimws(Sys.getenv("STEP9_MANHATTAN_MT_METHOD", "bonferroni")))
    use_fdr <- mt_raw %in% c("fdr", "bh")
  }
  fdr_alpha <- suppressWarnings(as.numeric(Sys.getenv("STEP9_MANHATTAN_FDR_ALPHA", "0.1")))
  if (!is.finite(fdr_alpha) || fdr_alpha <= 0 || fdr_alpha > 1) {
    fdr_alpha <- 0.1
  }

  m_tests <- nrow(man)
  bh_crit_p <- NA_real_
  y_fdr <- NA_real_
  n_fdr_sig <- 0L
  if (use_fdr && m_tests >= 1L) {
    pv <- man$Pvalue
    q_bh <- stats::p.adjust(pv, "BH")
    n_fdr_sig <- as.integer(sum(q_bh <= fdr_alpha, na.rm = TRUE))
    ps <- sort(pv)
    mlen <- length(ps)
    thr <- (seq_len(mlen) / mlen) * fdr_alpha
    ok <- ps <= thr
    if (any(ok)) {
      k <- max(which(ok))
      bh_crit_p <- ps[k]
      y_fdr <- -log10(bh_crit_p)
    }
  }
  fdr_line1 <- sprintf("BH-FDR (Benjamini-Hochberg): %d plotted tests, alpha=%g", m_tests, fdr_alpha)
  fdr_line2 <- if (!use_fdr) {
    ""
  } else if (is.na(bh_crit_p)) {
    sprintf("No threshold line (no BH rejections at alpha=%g)", fdr_alpha)
  } else {
    sprintf(
      "Purple: BH gate raw P=%.3g  ~  -log10=%.3f (%d points with BH q <= %g)",
      bh_crit_p, y_fdr, n_fdr_sig, fdr_alpha
    )
  }

  corner_txt <- if (!use_fdr) {
    paste(bonf_line1, bonf_line2, sep = "\n")
  } else {
    paste(fdr_line1, fdr_line2, sep = "\n")
  }

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    m <- man
    m$Group <- as.factor(m$Group)
    m$Max_MAF_lab <- factor(m$Max_MAF_lab)
    gg <- ggplot2::ggplot(m, ggplot2::aes(
      x = BP_cum, y = NegLog10P,
      color = Group, shape = Max_MAF_lab
    ))
    if (length(chr_vlines) > 0) {
      gg <- gg +
        ggplot2::geom_vline(xintercept = chr_vlines, colour = "grey88", linewidth = 0.35)
    }
    gg <- gg +
      ggplot2::geom_point(alpha = 0.55, size = 0.85, stroke = 0.25)

    if (!use_fdr) {
      gg <- gg +
        ggplot2::geom_hline(yintercept = y_bonf_gene, color = "#E69F00", size = 0.55, linetype = "dashed") +
        ggplot2::geom_hline(yintercept = y_bonf_maf, color = "#009E73", size = 0.55, linetype = "dashed")
    } else if (is.finite(y_fdr)) {
      gg <- gg +
        ggplot2::geom_hline(yintercept = y_fdr, color = "#7570B3", size = 0.55, linetype = "dashed")
    }

    gg <- gg +
      ggplot2::scale_x_continuous(
        breaks = midpoints,
        labels = chr_axis_labels,
        expand = ggplot2::expansion(mult = c(0.01, 0.018))
      ) +
      ggplot2::labs(
        title = mtitle,
        x = "Chromosome",
        y = expression(-log[10](P)),
        color = "Annotation group",
        shape = "max_MAF"
      ) +
      ggplot2::guides(x = ggplot2::guide_axis(check.overlap = FALSE)) +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::theme(
        panel.background = ggplot2::element_rect(fill = "white", colour = NA),
        plot.background = ggplot2::element_rect(fill = "white", colour = NA),
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        panel.border = ggplot2::element_rect(fill = NA, colour = "grey45", linewidth = 0.45),
        legend.position = "right",
        legend.justification = "center",
        legend.box = "vertical",
        legend.background = ggplot2::element_rect(fill = "white", colour = NA),
        axis.text.x = ggplot2::element_text(angle = 55, hjust = 1, vjust = 1, colour = "grey15", size = 8),
        axis.text.y = ggplot2::element_text(colour = "grey15"),
        plot.title = ggplot2::element_text(face = "bold"),
        plot.margin = ggplot2::margin(10, 14, 22, 10)
      )

    m_ord <- m[order(m$Pvalue), , drop = FALSE]
    top_df <- head(m_ord, min(as.integer(label_n), nrow(m_ord)))
    if (nrow(top_df) > 0) {
      if (requireNamespace("ggrepel", quietly = TRUE)) {
        gg <- gg + ggrepel::geom_text_repel(
          data = top_df,
          ggplot2::aes(x = BP_cum, y = NegLog10P, label = Gene),
          inherit.aes = FALSE,
          size = 2.6,
          color = "gray15",
          segment.size = 0.25,
          segment.color = "gray45",
          min.segment.length = 0,
          max.overlaps = 50,
          box.padding = 0.35,
          point.padding = 0.15,
          show.legend = FALSE
        )
      } else {
        gg <- gg + ggplot2::geom_text(
          data = top_df,
          ggplot2::aes(x = BP_cum, y = NegLog10P, label = Gene),
          inherit.aes = FALSE,
          size = 2.6,
          vjust = -0.35,
          color = "gray15",
          show.legend = FALSE
        )
      }
    }

    gg <- gg +
      ggplot2::annotate(
        "text",
        x = Inf,
        y = -Inf,
        label = corner_txt,
        hjust = 1,
        vjust = 0,
        size = 2.35,
        lineheight = 1.05,
        color = "gray38"
      )

    png(pngpath, width = 1950, height = 950, res = 150)
    print(gg)
    safe_dev_off()
  } else {
    png(pngpath, width = 1950, height = 950, res = 150)
    ug <- as.factor(man$Group)
    cols <- grDevices::rainbow(length(levels(ug)))[as.integer(ug)]
    nm <- length(unique(as.factor(man$Max_MAF_lab)))
    pchs <- c(16, 17, 15, 18, 8, 3, 4, 10)[as.integer(as.factor(man$Max_MAF_lab))]
    xr <- range(man$BP_cum)
    yr <- range(man$NegLog10P)
    padx <- if (diff(xr) > 0) diff(xr) * 0.012 else 1
    graphics::par(bg = "white", col.axis = "gray30", col.lab = "gray25")
    plot(
      NA_real_,
      NA_real_,
      xlim = c(xr[1] - padx, xr[2] + padx),
      ylim = yr,
      xlab = "Chromosome",
      ylab = expression(-log[10](P)),
      main = mtitle,
      xaxt = "n",
      frame.plot = TRUE
    )
    if (length(chr_vlines) > 0) {
      graphics::abline(v = chr_vlines, col = "grey90", lwd = 1)
    }
    graphics::points(man$BP_cum, man$NegLog10P, pch = pchs, cex = 0.42, col = cols)
    graphics::axis(1, at = midpoints, labels = chr_axis_labels, las = 2, tick = TRUE, cex.axis = 0.62)
    if (!use_fdr) {
      abline(h = y_bonf_gene, col = "#E69F00", lwd = 2, lty = 2)
      abline(h = y_bonf_maf, col = "#009E73", lwd = 2, lty = 2)
    } else if (is.finite(y_fdr)) {
      abline(h = y_fdr, col = "#7570B3", lwd = 2, lty = 2)
    }
    u <- graphics::par("usr")
    graphics::text(
      u[2], u[3],
      labels = corner_txt,
      adj = c(1, 0),
      cex = 0.45,
      col = "gray38",
      xpd = FALSE
    )
    top_n <- min(as.integer(label_n), nrow(man))
    if (top_n > 0) {
      top_genes <- man[order(man$Pvalue), ][seq_len(top_n), , drop = FALSE]
      graphics::text(top_genes$BP_cum, top_genes$NegLog10P, labels = top_genes$Gene,
           pos = 3, cex = 0.48, col = "gray20")
    }
    safe_dev_off()
  }

  if (write_group_stats && (tag == "allGroups_allMaf" || tag == "default")) {
    groups <- sort(unique(man$Group))
    stats_list <- lapply(groups, function(g) {
      d <- man[man$Group == g, ]
      data.frame(
        Group = g,
        N = nrow(d),
        Min_P = min(d$Pvalue, na.rm = TRUE),
        Median_P = median(d$Pvalue, na.rm = TRUE),
        N_GWS = sum(d$Pvalue < 5e-8, na.rm = TRUE),
        N_Sug = sum(d$Pvalue < 1e-5, na.rm = TRUE),
        N_Nom = sum(d$Pvalue < 0.05, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    })
    group_stats <- do.call(rbind, stats_list)
    write.table(group_stats, "group_statistics.txt", quote = FALSE, row.names = FALSE, sep = "\t")
  }

  invisible(TRUE)
}

qq_png <- paste0("qq_plot_", tag, ".png")
man_png_b <- paste0("manhattan_plot_", tag, "_bonferroni.png")
man_png_f <- paste0("manhattan_plot_", tag, "_fdr.png")
plot_qq(qq_fn, qq_png, paste("QQ —", tag))
plot_man(man_fn, man_png_b, paste("Manhattan (Bonferroni) —", tag), tag, n_genes_merged, "bonferroni", TRUE)
plot_man(man_fn, man_png_f, paste("Manhattan (BH-FDR) —", tag), tag, n_genes_merged, "fdr", FALSE)

if (legacy && file.exists(qq_png)) {
  file.copy(qq_png, "qq_plot.png", overwrite = TRUE)
}
if (legacy && file.exists(man_png_b)) {
  file.copy(man_png_b, "manhattan_plot.png", overwrite = TRUE)
}
if (legacy && file.exists(man_png_f)) {
  file.copy(man_png_f, "manhattan_plot_fdr.png", overwrite = TRUE)
}
if (tag == "allGroups_allMaf") {
  if (file.exists(qq_png)) {
    file.copy(qq_png, "qq_plot.png", overwrite = TRUE)
  }
  if (file.exists(man_png_b)) {
    file.copy(man_png_b, "manhattan_plot.png", overwrite = TRUE)
  }
  if (file.exists(man_png_f)) {
    file.copy(man_png_f, "manhattan_plot_fdr.png", overwrite = TRUE)
  }
}
RSCRIPT

    local _fg_save="${FILTER_GROUP:-}"
    local _fm_save="${FILTER_MAF:-}"

    regen_and_plot() {
        local tag="$1"
        local label="$2"
        echo "  $label"
        op_qqdata
        op_mandata
        if ! Rscript "$r_script" "$OUTPUT_DIR" "$tag" "qq_plot_data.txt" "manhattan_plot_data.txt" "${N_GENES_MERGED:-0}"; then
            echo "    WARNING: R plotting failed for tag=$tag"
            return 1
        fi
        [ -f "$OUTPUT_DIR/qq_plot_${tag}.png" ] && echo "    ✓ qq_plot_${tag}.png"
        [ -f "$OUTPUT_DIR/manhattan_plot_${tag}_bonferroni.png" ] && echo "    ✓ manhattan_plot_${tag}_bonferroni.png"
        [ -f "$OUTPUT_DIR/manhattan_plot_${tag}_fdr.png" ] && echo "    ✓ manhattan_plot_${tag}_fdr.png"
        return 0
    }

    if [ "$PLOT_DO_UNFILTERED" = "1" ]; then
        FILTER_GROUP=""
        FILTER_MAF=""
        regen_and_plot "allGroups_allMaf" "[all groups, all max_MAF] → qq + manhattan (_bonferroni + _fdr PNGs)"
    fi

    FILTER_GROUP="$_fg_save"
    FILTER_MAF="$_fm_save"

    if [ -f "$OUTPUT_DIR/group_statistics.txt" ]; then
        echo "  ✓ group_statistics.txt (from all-groups plot when enabled)"
    fi
    echo ""
}

#==========================================
# Operation Dispatcher
#==========================================

run_operation() {
    local op=$1
    
    case $op in
        mergechrom)
            op_mergechrom
            ;;
        mergeall)
            op_mergeall
            ;;
        listgroups)
            op_listgroups
            ;;
        findsig)
            op_findsig
            ;;
        findgws)
            op_findgws
            ;;
        findsug)
            op_findsug
            ;;
        findnom)
            op_findnom
            ;;
        top10)
            op_topgenes 10
            ;;
        top50)
            op_topgenes 50
            ;;
        top100)
            op_topgenes 100
            ;;
        chromsum)
            op_chromsum
            ;;
        groupsum)
            op_groupsum
            ;;
        fullsum)
            op_fullsum
            ;;
        qqdata)
            op_qqdata
            ;;
        mandata)
            op_mandata
            ;;
        plotdata)
            op_makeplots
            ;;
        makeplots)
            op_makeplots
            ;;
        standard)
            detect_files
            op_mergeall
            op_listgroups
            op_findsig
            op_topgenes 50
            op_groupsum
            op_makeplots
            op_fullsum
            ;;
        full)
            detect_files
            op_mergechrom
            op_mergeall
            op_listgroups
            op_findsig
            op_topgenes 10
            op_topgenes 50
            op_topgenes 100
            op_chromsum
            op_groupsum
            op_makeplots
            op_fullsum
            ;;
        quick)
            detect_files
            op_mergeall
            op_listgroups
            op_topgenes 50
            op_makeplots
            op_fullsum
            ;;
        *)
            echo "ERROR: Unknown operation '$op'" >&2
            echo ""
            show_menu
            exit 1
            ;;
    esac
}

#==========================================
# Main Execution
#==========================================

if [ "$INTERACTIVE" = true ]; then
    echo "Example: /data/run123/saige_gene_out   or   .   for this folder"
    read -r -p "Path to SAIGE-GENE results folder [${RESULT_FOLDER}]: " RF_PROMPT
    RESULT_FOLDER="${RF_PROMPT:-$RESULT_FOLDER}"
    finalize_results_directory "$RESULT_FOLDER"
    print_startup_banner
    
    show_menu
    echo ""
    
    echo "Examples: midpoint (usual) | start | end — used when gene coordinates are available."
    read -r -p "Position mode for gene plotting (start/end/midpoint) [midpoint]: " POS_MODE_INPUT
    POSITION_MODE="${POS_MODE_INPUT:-midpoint}"
    
    # Ask for gene coordinate file
    echo "Example: /refs/gencode_genes.tsv  (tab: chr, start, end, gene) — or Enter to skip."
    read -r -p "Gene coordinate file (optional): " GENE_COORD_INPUT
    if [ -n "$GENE_COORD_INPUT" ]; then
        GENE_COORD_FILE="$GENE_COORD_INPUT"
    else
        echo "No coordinate file provided."
        echo "Example: type 'ensembl' to fetch coordinates from Ensembl REST, or Enter for sequence order in results."
        read -r -p "Coordinate source [sequence/ensembl] (default: sequence): " COORD_INPUT
        COORD_SOURCE="${COORD_INPUT:-sequence}"
        
        if [ "$COORD_SOURCE" = "ensembl" ]; then
            echo "Example: 38 for GRCh38 / current REST; 37 for GRCh37 server."
            read -r -p "Ensembl genome build [37/38] (default: 38): " ENS_BUILD_INPUT
            ENSEMBL_BUILD="${ENS_BUILD_INPUT:-38}"
            
            echo "Example: 115 — optional; server lists available releases via /info/data."
            read -r -p "Ensembl release (optional, e.g. 115; Enter to skip): " ENS_REL_INPUT
            ENSEMBL_RELEASE="${ENS_REL_INPUT:-}"
        else
            COORD_SOURCE="sequence"
        fi
    fi
    
    # Ask for annotation group filter
    echo "Example: lof  or  lof;missense  — filters analysis tables only."
    read -r -p "Filter by annotation group (optional; Enter = no filter): " FILTER_INPUT
    if [ -n "$FILTER_INPUT" ]; then
        FILTER_GROUP="$FILTER_INPUT"
    fi
    
    # Ask for MAF filter
    echo "Example: 0.01  (keep tests with max_MAF <= 0.01). Enter = no filter."
    read -r -p "Filter by max_MAF threshold (optional): " MAF_INPUT
    if [ -n "$MAF_INPUT" ]; then
        FILTER_MAF="$MAF_INPUT"
    fi
    
    echo ""
    echo "--- PNG plots (makeplots / plotdata): one QQ + one Manhattan (all groups, all max_MAF) ---"
    echo "  qq_plot_allGroups_allMaf.png"
    echo "  manhattan_plot_allGroups_allMaf_bonferroni.png   manhattan_plot_allGroups_allMaf_fdr.png"
    read -r -p "Generate those plots when using plotdata/makeplots? [Y/n] (Enter=yes): " PU
    case "${PU:-Y}" in
        n|N|no|NO)
            PLOT_DO_UNFILTERED=0
            ;;
        *)
            PLOT_DO_UNFILTERED=1
            ;;
    esac
    
    echo ""
    echo "Examples: standard | quick | plotdata | mergeall+findsig+top50"
    read -r -p "Enter operations: " OPERATIONS
    
    if [ -z "$OPERATIONS" ]; then
        echo "No operations specified. Exiting."
        exit 0
    fi
fi

# Plot defaults from environment when not set interactively
if [ -z "${PLOT_DO_UNFILTERED:-}" ]; then
    PLOT_DO_UNFILTERED="${STEP9_PLOT_UNFILTERED:-1}"
fi

# Normalize user inputs
COORD_SOURCE=$(echo "$COORD_SOURCE" | tr '[:upper:]' '[:lower:]')
POSITION_MODE=$(echo "$POSITION_MODE" | tr '[:upper:]' '[:lower:]')
if [[ "$POSITION_MODE" != "start" && "$POSITION_MODE" != "end" && "$POSITION_MODE" != "midpoint" ]]; then
    echo "Warning: Invalid position mode '$POSITION_MODE', using midpoint."
    POSITION_MODE="midpoint"
fi
if [[ "$ENSEMBL_BUILD" != "37" && "$ENSEMBL_BUILD" != "38" ]]; then
    echo "Warning: Invalid Ensembl build '$ENSEMBL_BUILD', using 38."
    ENSEMBL_BUILD="38"
fi

# Detect files first (needed before optional Ensembl coordinate generation)
detect_files

# Setup gene coordinates
if [ -n "$GENE_COORD_FILE" ]; then
    setup_gene_coords
elif [ "$COORD_SOURCE" = "ensembl" ]; then
    setup_gene_coords_ensembl || echo "  Falling back to sequence-based coordinates."
fi

# Parse and execute operations
IFS='+' read -ra OPS <<< "$OPERATIONS"

echo "Executing operations: ${OPS[*]}"
if [ -n "$FILTER_GROUP" ]; then
    echo "Filtering by annotation group: $FILTER_GROUP"
fi
if [ -n "$FILTER_MAF" ]; then
    echo "Filtering by max_MAF <= $FILTER_MAF"
fi
echo "PNG plots (makeplots): single QQ + Manhattan (all groups × all max_MAF), unfiltered=$PLOT_DO_UNFILTERED"
echo ""

for op in "${OPS[@]}"; do
    run_operation "$op"
done

#==========================================
# Final Summary
#==========================================

echo "=========================================="
echo "Analysis Complete"
echo "=========================================="
echo "Completed: $(date)"
echo ""

if [ -f "$OUTPUT_DIR/analysis_summary.txt" ]; then
    echo "Summary report available at:"
    echo "  $OUTPUT_DIR/analysis_summary.txt"
    echo ""
fi

echo "Generated files in $OUTPUT_DIR:"
ls -lh "$OUTPUT_DIR"/*.txt 2>/dev/null | grep -v "chunk" | awk '{print "  " $9 " (" $5 ")"}'
echo ""

#==========================================
# Next Steps Suggestions
#==========================================

echo "=========================================="
echo "Next Steps"
echo "=========================================="
echo ""
echo "Primary outputs:"
echo "  - analysis_summary.txt"
echo "  - top50_genes.txt"
echo "  - qq_plot_data.txt and manhattan_plot_data.txt"
echo "  - qq_plot.png; manhattan_plot.png = Bonferroni; manhattan_plot_fdr.png = BH-FDR (when Rscript is available)"
echo ""

if [ -f "$OUTPUT_DIR/annotation_groups.txt" ]; then
    echo "Available annotation groups:"
    awk 'NR<=10 {print}' "$OUTPUT_DIR/annotation_groups.txt"
    echo ""
fi

echo "Automatic plot generation:
  - Included in standard, full, quick, and plotdata
  - Or run directly: ./step9_analyze_results.sh makeplots"

echo ""
echo "=========================================="
echo "Analysis script completed successfully"
echo "=========================================="
echo ""

exit 0

"""
02_download_opentargets.py
==========================
Purpose : Download the Open Targets Platform datasets needed for the project
          and write a schema report (column names, row counts, sample rows) so
          the next scripts can be written against the REAL columns.

Datasets (release set in RELEASE):
  disease                         disease IDs / names / ontology
  target                          Ensembl gene ID -> gene symbol (for mapping)
  drug_molecule                   drug IDs / names
  drug_mechanism_of_action        drug -> target genes (drug targets)
  clinical_indication             drug -> disease (known indications)
  association_by_datatype_direct  gene-disease scores per EVIDENCE TYPE. We will
                                  later keep only GENETIC evidence types, so drug
                                  information never leaks into disease gene lists.

Output  : data/raw/opentargets/<dataset>/*.parquet     (not committed to GitHub)
          results/opentargets_download_manifest.csv
          results/opentargets_schema_report.txt
Safety  : a dataset larger than its MAX_GB cap is skipped (not downloaded) and
          reported, so you never silently pull tens of GB.
Re-run  : safe. Files already downloaded with the right size are skipped.
"""

import csv
import re
import time
from pathlib import Path
from urllib.parse import unquote, urljoin
from urllib.request import Request, urlopen

# ------------------------------------------------------------------
# 1. SETTINGS
# ------------------------------------------------------------------
BASE_DIR = Path(
    r"C:\Users\Vyomini\Documents\MU Engineering_Btech [CB]"
    r"\Projects\Course_Projects\Sem7_Projects"
    r"\Sys Bio & Network Model Project"
    r"\Network_Drug_Repurposing"
)

RELEASE  = "26.09"
BASE_URL = f"https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/{RELEASE}/output/"

# dataset name -> maximum size (GB) we are willing to download
DATASETS = {
    "disease": 1,
    "target": 3,
    "drug_molecule": 1,
    "drug_mechanism_of_action": 1,
    "clinical_indication": 2,
    "association_by_datatype_direct": 8,
}

RAW_DIR     = BASE_DIR / "data" / "raw" / "opentargets"
RESULTS_DIR = BASE_DIR / "results"
RAW_DIR.mkdir(parents=True, exist_ok=True)
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

HEADERS = {"User-Agent": "Mozilla/5.0 (course project downloader)"}
TIMEOUT = 120
RETRIES = 3
MAX_DEPTH = 3          # how many sub-folder levels to crawl inside a dataset

report_lines = []      # everything printed is also saved to the schema report


def log(text=""):
    print(text)
    report_lines.append(str(text))


# ------------------------------------------------------------------
# 2. HELPERS: list a web folder, get a file size, download a file
# ------------------------------------------------------------------
def fetch_text(url):
    req = Request(url, headers=HEADERS)
    with urlopen(req, timeout=TIMEOUT) as r:
        return r.read().decode("utf-8", errors="replace")


def list_parquet_files(url, depth=0):
    """Return URLs of all .parquet files under a folder URL (recursive)."""
    html = fetch_text(url)
    found = []
    for href in re.findall(r'href="([^"]+)"', html):
        # skip sort links, parent folder, absolute paths, anchors
        if href.startswith(("?", "/", "#", "http", "..")):
            continue
        full = urljoin(url, href)
        if href.endswith("/"):
            if depth < MAX_DEPTH:
                found += list_parquet_files(full, depth + 1)
        elif href.lower().endswith(".parquet"):
            found.append(full)
    return found


def remote_size(url):
    """File size in bytes from a HEAD request (None if the server does not say)."""
    try:
        req = Request(url, headers=HEADERS, method="HEAD")
        with urlopen(req, timeout=TIMEOUT) as r:
            length = r.headers.get("Content-Length")
            return int(length) if length else None
    except Exception:
        return None


def download_file(url, dest, expected_size):
    """Download url -> dest with retries; skip if already complete."""
    if dest.exists() and expected_size and dest.stat().st_size == expected_size:
        return "skipped (already complete)"
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_name(dest.name + ".part")
    for attempt in range(1, RETRIES + 1):
        try:
            req = Request(url, headers=HEADERS)
            with urlopen(req, timeout=TIMEOUT) as r, open(tmp, "wb") as f:
                while True:
                    chunk = r.read(1024 * 1024)
                    if not chunk:
                        break
                    f.write(chunk)
            tmp.replace(dest)
            return "downloaded"
        except Exception as err:
            log(f"    attempt {attempt}/{RETRIES} failed: {type(err).__name__}: {err}")
            time.sleep(3)
    return "FAILED"


# ------------------------------------------------------------------
# 3. DOWNLOAD EACH DATASET
# ------------------------------------------------------------------
manifest = []
log("=" * 70)
log(f"OPEN TARGETS RELEASE {RELEASE} - DOWNLOAD")
log("=" * 70)

for dataset, cap_gb in DATASETS.items():
    log(f"\n--- {dataset} ---")
    dataset_url = urljoin(BASE_URL, dataset + "/")
    try:
        files = list_parquet_files(dataset_url)
    except Exception as err:
        log(f"  Could not list the folder: {type(err).__name__}: {err}")
        manifest.append([dataset, 0, "", "listing failed"])
        continue

    sizes = [remote_size(u) for u in files]
    known = [s for s in sizes if s]
    total_gb = sum(known) / 1024 ** 3
    log(f"  parquet files: {len(files)}   total size: {total_gb:.2f} GB "
        f"(cap {cap_gb} GB)" + ("" if len(known) == len(files) else "   [some sizes unknown]"))

    if not files:
        manifest.append([dataset, 0, "", "no parquet files found"])
        continue
    if total_gb > cap_gb:
        log(f"  SKIPPED: larger than the {cap_gb} GB cap. Tell me and we will filter it instead.")
        manifest.append([dataset, len(files), round(total_gb, 2), "skipped: over cap"])
        continue

    statuses = {"downloaded": 0, "skipped (already complete)": 0, "FAILED": 0}
    for i, (url, size) in enumerate(zip(files, sizes), 1):
        rel = unquote(url[len(dataset_url):])
        status = download_file(url, RAW_DIR / dataset / rel, size)
        statuses[status] += 1
        if i % 10 == 0 or i == len(files):
            log(f"  {i}/{len(files)} files processed")
    log(f"  result: {statuses}")
    manifest.append([dataset, len(files), round(total_gb, 2),
                    "OK" if statuses["FAILED"] == 0 else "some files FAILED"])

with open(RESULTS_DIR / "opentargets_download_manifest.csv", "w",
          newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow(["dataset", "n_files", "size_gb", "status"])
    w.writerows(manifest)

# ------------------------------------------------------------------
# 4. SCHEMA REPORT (needs pyarrow)
# ------------------------------------------------------------------
log("\n" + "=" * 70)
log("SCHEMA REPORT  (send this to me)")
log("=" * 70)

try:
    import pyarrow.compute as pc
    import pyarrow.dataset as pads

    for dataset in DATASETS:
        folder = RAW_DIR / dataset
        log(f"\n### {dataset}")
        if not folder.exists() or not any(folder.rglob("*.parquet")):
            log("  (not downloaded)")
            continue
        try:
            d = pads.dataset(str(folder), format="parquet", partitioning="hive")
            names = d.schema.names
            log(f"  rows: {d.count_rows():,}")
            log("  columns:")
            for field in d.schema:
                log(f"    - {field.name}: {str(field.type)[:110]}")
            # sample rows
            for k, row in enumerate(d.head(2).to_pylist(), 1):
                log(f"  sample row {k}: {str(row)[:700]}")
            # distinct values of a few key label columns, if they exist
            for col in ("datatypeId", "datasourceId", "clinicalStage", "maxClinicalStage",
                        "maximumClinicalStage", "phase", "status", "drugType"):
                if col in names:
                    vals = pc.unique(d.to_table(columns=[col]).column(col)).to_pylist()
                    log(f"  distinct {col} ({len(vals)}): {str(sorted(map(str, vals)))[:600]}")
        except Exception as err:
            log(f"  could not read: {type(err).__name__}: {err}")
except ImportError:
    log("pyarrow is not installed: run  conda install -c conda-forge pyarrow")

(RESULTS_DIR / "opentargets_schema_report.txt").write_text(
    "\n".join(report_lines), encoding="utf-8")
log("\nSaved: results/opentargets_schema_report.txt")
log("Script 02 finished.")
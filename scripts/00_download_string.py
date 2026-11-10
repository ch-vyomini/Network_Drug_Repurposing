"""
00_download_string.py
=====================
Purpose : Download the STRING v12.5 human files needed by script 01 into data/raw/string/.
          (Open Targets data is downloaded by script 02.)
Files   : 9606.protein.physical.links.v12.5.txt.gz   physical protein-protein interactions with scores
          9606.protein.info.v12.5.txt.gz             STRING protein id -> preferred gene symbol
Notes   : * If a matching file is already in data/raw/string/ the download is skipped.
          * The URLs follow STRING's usual naming pattern. If a download fails, get the same files by hand
            from https://string-db.org/cgi/download  (organism: Homo sapiens, version 12.5) and put them
            in data/raw/string/ - script 01 finds them by name.
          * Script 01 keeps only interactions with combined_score >= 400, so the unfiltered file gives the same networks.
"""
import os
import shutil
from pathlib import Path
from urllib.request import Request, urlopen

BASE_DIR = Path(os.environ.get("NDR_BASE", Path(__file__).resolve().parents[1]))
OUT = BASE_DIR / "data" / "raw" / "string"
OUT.mkdir(parents=True, exist_ok=True)

VERSION = "12.5"
FILES = {   # glob pattern that means 'already have it'  ->  download URL
    "9606.protein.physical.links*": f"https://stringdb-downloads.org/download/protein.physical.links.v{VERSION}/9606.protein.physical.links.v{VERSION}.txt.gz",
    "9606.protein.info*": f"https://stringdb-downloads.org/download/protein.info.v{VERSION}/9606.protein.info.v{VERSION}.txt.gz",
}

for pattern, url in FILES.items():
    have = [p for p in OUT.glob(pattern) if p.is_file()]
    if have:
        print(f"already present, skipped: {have[0].name}")
        continue
    dest = OUT / url.split("/")[-1]
    print("downloading", url)
    try:
        req = Request(url, headers={"User-Agent": "Mozilla/5.0 (course project downloader)"})
        with urlopen(req, timeout=300) as r, open(str(dest) + ".part", "wb") as f:
            shutil.copyfileobj(r, f, 1024 * 1024)
        os.replace(str(dest) + ".part", dest)
        print("  saved", dest.name, f"({dest.stat().st_size / 1e6:.0f} MB)")
    except Exception as err:
        print(f"  FAILED: {type(err).__name__}: {err}")
        print("  -> download the file by hand from https://string-db.org/cgi/download "
              f"(Homo sapiens, v{VERSION}) and place it in {OUT}")
print("Script 00 finished.")

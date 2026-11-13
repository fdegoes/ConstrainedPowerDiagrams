#!/usr/bin/env python3
"""Downloads the parts of the public 3D grain datasets that the reconstruction
tests need (per-grain volumes and centroids, or a segmented grain label map
to compute them from). The record ids come from an unverified literature
search: run with --list first and check each record's title, license and
files before downloading.

  zenodo 14094899  Nugget sandstone, 3DXRD + X-ray CT: the per-grain 3DXRD
                   table (3DXRD.mat, 0.2 MB) and one CT segmentation
                   (Step3_XRCT_Segmentation.tif, 83 MB)
  zenodo 20323210  pure iron, lab DCT, 5 time steps: one zipped grain map,
                   T0_AGG100_crop_gm3d.h5.zip (808 MB)
  zenodo 3785085   Zircaloy-4 + 7Nb, serial-section 3D EBSD: only the
                   metadata by default; the segmented DREAM.3D file is
                   dream3d.zip (3.4 GB), fetched with
                   --record 3785085 --files dream3d.zip

Each record goes into its own new directory <out>/zenodo_<id>/ with a
record.json of its metadata (title, license, DOI, file list). Files are
checked against their Zenodo MD5 checksums. Only the standard library is
used.

Run it outside the Claude Code sandbox, which blocks zenodo.org:
  python3 external/download.py --list
  python3 external/download.py [--record ID ...] [--files NAME ...]
      [--max-gb 5] [--out data/grains]
"""

import argparse
import fnmatch
import hashlib
import json
import os
import sys
import urllib.request

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
API = "https://zenodo.org/api/records/%d"
DOCS = ["*readme*", "*README*", "*.md", "*.txt", "*.pdf"]

# Per record: the file patterns to download, and at most how many of the
# data files (in name order) to take; the documentation is always taken.
RECORDS = {
    14094899: dict(name="Nugget sandstone, 3DXRD + CT", data=["3DXRD.mat", "Step3_XRCT_Segmentation.tif"],
                   maxData=None),
    20323210: dict(name="pure iron, lab DCT", data=["*_gm3d.h5.zip"], maxData=1),
    3785085: dict(name="Zircaloy-4 + 7Nb, 3D EBSD", data=[], maxData=None),
}


def get_json(url):
    with urllib.request.urlopen(url, timeout=60) as r:
        return json.load(r)


def matches(name, patterns):
    return any(fnmatch.fnmatch(name.lower(), p.lower()) for p in patterns)


def select(files, spec, explicit):
    """The files to download: the explicit names if given, otherwise the
    documentation and the first maxData data files."""
    if explicit:
        return [f for f in files if f["key"] in explicit]
    docs = [f for f in files if matches(f["key"], DOCS)]
    data = sorted((f for f in files if matches(f["key"], spec["data"])), key=lambda f: f["key"])
    if spec["maxData"] is not None:
        data = data[:spec["maxData"]]
    return docs + [f for f in data if f not in docs]


def download(f, path):
    """Streams f to path (through path.part), checking its MD5."""
    url = f["links"].get("content") or f["links"]["self"]
    algo, _, want = f.get("checksum", ":").partition(":")
    h = hashlib.md5()
    tmp = path + ".part"
    with urllib.request.urlopen(url, timeout=120) as r, open(tmp, "wb") as out:
        done = 0
        while True:
            chunk = r.read(1 << 20)
            if not chunk:
                break
            out.write(chunk)
            h.update(chunk)
            done += len(chunk)
            print("\r     %.1f / %.1f MB" % (done / 1e6, f["size"] / 1e6), end="", flush=True)
    print()
    if algo == "md5" and h.hexdigest() != want:
        os.remove(tmp)
        raise RuntimeError("checksum mismatch for %s" % f["key"])
    os.replace(tmp, path)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--record", nargs="+", type=int, default=list(RECORDS), choices=list(RECORDS))
    ap.add_argument("--files", nargs="+", default=None, help="download these file names instead of the default selection")
    ap.add_argument("--list", action="store_true", help="only print each record's metadata and files")
    ap.add_argument("--max-gb", type=float, default=5.0, help="skip files larger than this (GB)")
    ap.add_argument("--out", default=os.path.join(ROOT, "data", "grains"))
    args = ap.parse_args()

    for rid in args.record:
        spec = RECORDS[rid]
        print("\n=== zenodo %d: %s" % (rid, spec["name"]))
        try:
            rec = get_json(API % rid)
        except Exception as e:  # noqa: BLE001
            print("  cannot read the record: %s" % e)
            continue
        md = rec.get("metadata", {})
        lic = md.get("license") or md.get("rights") or {}
        print("  title:   %s" % md.get("title"))
        print("  license: %s" % (lic.get("id") if isinstance(lic, dict) else lic))
        print("  doi:     %s" % rec.get("doi"))
        files = rec.get("files", [])
        chosen = select(files, spec, args.files)
        for f in files:
            mark = "*" if f in chosen else " "
            print("  %s %10.1f MB  %s" % (mark, f["size"] / 1e6, f["key"]))
        if args.list:
            continue

        d = os.path.join(args.out, "zenodo_%d" % rid)
        os.makedirs(d, exist_ok=True)
        with open(os.path.join(d, "record.json"), "w") as out:
            json.dump(dict(id=rid, doi=rec.get("doi"), title=md.get("title"), license=lic,
                           files=[dict(key=f["key"], size=f["size"], checksum=f.get("checksum")) for f in files]),
                      out, indent=1)
        for f in chosen:
            path = os.path.join(d, os.path.basename(f["key"]))
            if os.path.exists(path):
                print("  -> %s exists, skipped" % f["key"])
            elif f["size"] > args.max_gb * 1e9:
                print("  -> %s is %.1f GB, above --max-gb, skipped" % (f["key"], f["size"] / 1e9))
            else:
                print("  -> %s" % f["key"])
                try:
                    download(f, path)
                except Exception as e:  # noqa: BLE001
                    print("     failed: %s" % e, file=sys.stderr)
    print("\nFiles are in %s" % args.out)


if __name__ == "__main__":
    main()

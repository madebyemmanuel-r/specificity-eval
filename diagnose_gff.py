"""
diagnose_gff.py
===============
Prints the raw GFF response for known genic coordinates so we can
inspect the exact attribute format NCBI returns. Run this first.

Usage:
    python diagnose_gff.py
    python diagnose_gff.py --api-key YOUR_KEY
"""

import argparse
import time
import requests


def fetch_gff(accession, start, end, api_key=None):
    time.sleep(0.15 if api_key else 0.4)
    params = {
        "db": "nuccore",
        "id": accession,
        "seq_start": start,
        "seq_stop": end,
        "rettype": "gff",
        "retmode": "text",
    }
    if api_key:
        params["api_key"] = api_key
    r = requests.get(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
        params=params,
        timeout=15.0,
    )
    return r.status_code, r.text


def main(api_key=None):
    cases = [
        ("fliC core window",    "NC_000913.3", 2002500, 2002800),
        ("UL27 core window",    "NC_001806.2", 53500,   55000),
        ("Invalid accession",   "XX_GARBAGE",  1,       100),
    ]

    for label, acc, start, end in cases:
        print(f"\n{'='*65}")
        print(f"  {label}  ({acc}:{start}-{end})")
        print(f"{'='*65}")
        status, text = fetch_gff(acc, start, end, api_key)
        print(f"  HTTP {status}")
        # Print full response (strip ##FASTA for brevity)
        gff_body = text.split("##FASTA")[0]
        print(gff_body[:3000])
        if len(gff_body) > 3000:
            print(f"  ... [truncated, total {len(gff_body)} chars]")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--api-key", dest="api_key", default=None)
    args = parser.parse_args()
    main(args.api_key)
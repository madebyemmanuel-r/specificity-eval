"""
lookup_intersecting_gene.py
============================
Queries the NCBI E-utilities efetch API (GenBank flatfile format) to identify 
the gene or genomic feature overlapping a given coordinate window on a 
RefSeq nucleotide accession.

Design & Formatting Notes
-------------------------
Although rettype=gff is documented by NCBI, the API silently falls back to 
GenBank flatfile format for coordinate-bounded nuccore requests. This module 
directly parses GenBank flatfile syntax (/qualifier=) using optimised regexes.

Priority for name extraction: /gene= > /product= > /locus_tag=
"""

import re
import time
import requests
from typing import Optional, Final

# Compiled regexes for performance across repeat calls
_RE_GENE: Final = re.compile(r'/gene="([^"]+)"')
_RE_PRODUCT: Final = re.compile(r'/product="([^"]+)"')
_RE_LOCUS_TAG: Final = re.compile(r'/locus_tag="([^"]+)"')

# Qualifier values that are structural labels, not informative gene symbols
_UNINFORMATIVE_GENES: Final = frozenset({
    "gene", "cds", "mrna", "exon", "region", "regulatory", "recombination_region"
})

def _parse_genbank_features(text: str) -> Optional[str]:
    """
    Extract the best gene name from a GenBank flatfile response using 
    strict feature qualifier priority hierarchy.

    Parameters
    ----------
    text : str
        The raw text response from NCBI efetch.

    Returns
    -------
    Optional[str]
        The extracted feature identifier, or None if the region is intergenic
        or the response contains an API error message.
    """
    # NCBI returns HTTP 200 with space-padded ASCII errors for invalid accessions.
    # Collapse whitespace to flag errors early and safely.
    collapsed = text.replace(" ", "").replace("\n", "").replace("\r", "")
    if collapsed.startswith("Error:") or "Failedtounderstand" in collapsed:
        return None

    # Priority 1: /gene="SYMBOL" — Canonical, concise gene symbols
    for match in _RE_GENE.finditer(text):
        val = match.group(1).strip()
        if val.lower() not in _UNINFORMATIVE_GENES:
            return val

    # Priority 2: /product="DESCRIPTION" — Informative protein/RNA product details
    for match in _RE_PRODUCT.finditer(text):
        val = match.group(1).strip()
        if val.lower() not in _UNINFORMATIVE_GENES:
            return val

    # Priority 3: /locus_tag="TAG" — Systematic genomic identifier fallback
    for match in _RE_LOCUS_TAG.finditer(text):
        val = match.group(1).strip()
        if val.lower() not in _UNINFORMATIVE_GENES:
            return val

    return None


def lookup_intersecting_gene(
    accession: str,
    hit_start: int,
    hit_end: int,
    api_key: Optional[str] = None,
) -> str:
    """
    Identify the gene or genomic feature at a coordinate window via NCBI.

    Normalises coordinate orientation, applies adaptive API rate throttling,
    and safely catches network, timeout, or schema errors.

    Parameters
    ----------
    accession : str
        RefSeq nucleotide accession (e.g., "NC_000913.3").
    hit_start : int
        Start coordinate of the window (1-based index).
    hit_end : int
        End coordinate of the window (1-based index).
    api_key : Optional[str], default=None
        NCBI API key. Authenticated requests use a 0.15s throttle (~7 req/s);
        unauthenticated requests throttle at 0.40s (~2.5 req/s).

    Returns
    -------
    str
        Gene symbol, product, or locus tag if an annotated feature overlaps 
        the window; "Intergenic Region" if empty. On failure, returns an explicit 
        "Lookup Error (...)" string.
    """
    # Input sanity validation
    if not accession or not isinstance(accession, str):
        return "Lookup Error (Invalid Accession Argument)"
    
    # Adaptive rate-limiting throttle to respect NCBI guidelines
    time.sleep(0.15 if api_key else 0.40)

    # Normalise orientation dynamically (BLAST hits can be flipped)
    hit_min = min(hit_start, hit_end)
    hit_max = max(hit_start, hit_end)

    # Protect against invalid or zero/negative coordinates
    if hit_min <= 0:
        hit_min = 1

    params = {
        "db": "nuccore",
        "id": accession.strip(),
        "seq_start": hit_min,
        "seq_stop": hit_max,
        "rettype": "gb",      # Explicitly targeting GenBank flatfile
        "retmode": "text",
    }
    if api_key:
        params["api_key"] = api_key.strip()

    efetch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

    try:
        res = requests.get(efetch_url, params=params, timeout=15.0)

        if res.status_code == 429:
            return "Lookup Error (Rate Limited)"
        if res.status_code != 200:
            return f"Lookup Error (HTTP {res.status_code})"

        # Handle empty payloads safely
        if not res.text or res.text.isspace():
            return "Lookup Error (Empty NCBI Response)"

        result = _parse_genbank_features(res.text)

        if result is None:
            # Differentiate a real biological gap from a quiet API crash layout
            collapsed = res.text.replace(" ", "").replace("\n", "").replace("\r", "")
            if "Failedtounderstand" in collapsed or collapsed.startswith("Error:"):
                return "Lookup Error (Bad Accession)"
            return "Intergenic Region"

        return result

    except requests.exceptions.Timeout:
        return "Lookup Error (Timeout)"
    except requests.exceptions.ConnectionError:
        return "Lookup Error (Connection Failed)"
    except Exception as e:
        # Catch-all fallback for unhandled exceptions or unexpected runtime drops
        return f"Lookup Error (Unexpected: {type(e).__name__})"

"""
test_lookup_gene.py
===================
End-to-end integration tests for lookup_intersecting_gene() against the live 
NCBI E-utilities API.

Ground-truth coordinates sourced directly from NCBI Gene records:
  fliC  → Gene ID 949101   → NC_000913.3 : 2002110..2003606 (complement)
  UL27  → Gene ID 24271469 → NC_001806.2 : 53059..55795    (complement)
  UL28  → Gene ID 2703457  → NC_001806.2 : 53059..58160    (complement)
  UL30  → Gene ID 2703459  → NC_001806.2 : 63799..67680    (complement)
"""

import re
import time
import pytest
import requests

# ---------------------------------------------------------------------------
# Fixtures & Helpers
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def api_key(request):
    """Optionally reads API key from command line option or env fallback."""
    return request.config.getoption("--api-key", default=None)


def _fetch_gb_raw(accession: str, start: int, end: int, api_key: str = None) -> str:
    """Fetch raw flatfile directly to discover features independent of production parsing logic."""
    time.sleep(0.15 if api_key else 0.40)
    params = {
        "db": "nuccore",
        "id": accession,
        "seq_start": start,
        "seq_stop": end,
        "rettype": "gb",
        "retmode": "text",
    }
    if api_key:
        params["api_key"] = api_key
    r = requests.get(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
        params=params,
        timeout=15.0,
    )
    return r.text


# ---------------------------------------------------------------------------
# Group 1 & 2: Targeted Core Detection & Coordinate Robustness
# ---------------------------------------------------------------------------

@pytest.mark.parametrize(
    "label, accession, start, end, expected_options, is_startswith",
    [
        # 1a: Standard crisp window inside fliC
        ("fliC_core", "NC_000913.3", 2002500, 2002800, ["fliC"], False),
        # 1b: Full coordinate window boundary
        ("fliC_full", "NC_000913.3", 2002110, 2003606, ["fliC"], False),
        # 1c: Overlapping architecture window (UL28 overlaps with UL27)
        ("UL27_window", "NC_001806.2", 53500, 55000, ["ul27", "ul28"], False),
        # 1d: Isolated unambiguous viral gene
        ("UL30_core", "NC_001806.2", 64500, 66000, ["ul30"], False),
        # 2a: Flipped coordinates orientation sanity check
        ("reversed_coords", "NC_000913.3", 2003606, 2002110, ["fliC"], False),
        # 2b: Single-base window snapshot inside feature boundary
        ("single_nucleotide", "NC_000913.3", 2002500, 2002500, ["fliC"], False),
        # 3: Eukaryotic / Intronic Edge Cases - BRCA1 Intron 10 window (validates that it parses the encompassing gene envelope)
        ("BRCA1_intronic", "NC_000017.11", 43095000, 43095500, ["BRCA1"], False),
        # 4a: Garbage input layout error handling
        ("invalid_accession", "XX_GARBAGE999", 1, 100, ["lookup error"], True),
    ],
)

def test_lookup_cases(label, accession, start, end, expected_options, is_startswith, api_key):
    res = lookup_intersecting_gene(accession, start, end, api_key=api_key)
    res_lower = res.lower()
    
    if is_startswith:
        assert any(res_lower.startswith(opt.lower()) for opt in expected_options), f"Failed on {label}: got {res}"
    else:
        assert any(opt.lower() in res_lower for opt in expected_options), f"Failed on {label}: got {res}"

# ---------------------------------------------------------------------------
# Group 3: Dynamic Intergenic Discovery
# ---------------------------------------------------------------------------

def test_intergenic_region_dynamic(api_key):
    """Queries absolute start coordinates to extract intervals, sorting out a gap dynamically."""
    accession = "NC_000913.3"
    text = _fetch_gb_raw(accession, 1, 4000, api_key=api_key)
    
    raw_locs = re.findall(
        r'^\s{5}(?:gene|CDS)\s+(?:complement\()?(\d+)\.\.(\d+)',
        text,
        re.MULTILINE,
    )
    intervals = sorted([(int(s), int(e)) for s, e in raw_locs])
    
    gap_start, gap_end = 270, 280  # Hardcoded static fallback (thrL / thrA gap)

    for i in range(len(intervals) - 1):
        current_end = intervals[i][1]
        next_start = intervals[i+1][0]
        if (next_start - current_end) >= 20:
            mid = (current_end + next_start) // 2
            gap_start, gap_end = mid - 5, mid + 5
            break

    res = lookup_intersecting_gene(accession, gap_start, gap_end, api_key=api_key)
    assert "intergenic region" in res.lower(), f"Expected intergenic window at {gap_start}-{gap_end}, got: {res}"


# ---------------------------------------------------------------------------
# Group 4: Edge Case Robustness Layouts
# ---------------------------------------------------------------------------

def test_zero_width_coordinate(api_key):
    res = lookup_intersecting_gene("NC_001806.2", 1, 1, api_key=api_key)
    assert isinstance(res, str) and len(res) > 0


def test_massive_window_parsing(api_key):
    res = lookup_intersecting_gene("NC_001806.2", 1, 50000, api_key=api_key)
    assert isinstance(res, str) and len(res) > 0


# ---------------------------------------------------------------------------
# Group 5: Rate Limiter Consistency Checks
# ---------------------------------------------------------------------------

def test_api_call_consistency(api_key):
    """Executes consecutive validation calls to ensure parsing maps accurately each pass."""
    repeats = [
        lookup_intersecting_gene("NC_000913.3", 2002500, 2002800, api_key=api_key)
        for _ in range(10)
    ]
    assert len(set(repeats)) == 1, f"Inconsistent returns generated across requests: {repeats}"
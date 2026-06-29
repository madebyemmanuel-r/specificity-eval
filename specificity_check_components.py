from __future__ import annotations
from collections import defaultdict
from dataclasses import dataclass, field, asdict
from functools import lru_cache
from pathlib import Path
from typing import Optional, Dict, List, Any
from NCBI import ncbi_session
import logging
import math
import time
import requests
import xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import re
import matplotlib

matplotlib.use('Agg')  # Headless backend safety for server/headless execution

# Import our native thermodynamic matrix data
from thermo_data import NEAREST_NEIGHBOUR_DATA, TERMINAL_CORRECTIONS

logger = logging.getLogger(__name__)

@dataclass(slots=True)
class AssayConditions:
    """Standard RPA reaction assay conditions for calculating nearest-neighbour thermodynamics."""
    monovalent_salt_conc: float = 100.0  # mM (KAc)
    divalent_salt_conc: float = 14.0     # mM (MgAc)
    dntp_conc: float = 1.8               # mM
    oligo_conc: float = 0.5              # uM
    template_conc: float = 1e-6          # uM (1 pM)

@dataclass(slots=True)
class BlastHit:
    """Represents a single parsed BLAST hit, unified with profile metrics and alignment parameters."""
    seq_id: str
    primer_direction: str
    pct_identity: float
    q_start: int
    q_end: int
    back_seq_start: int
    back_seq_end: int
    evalue: float
    num_gaps: int
    strand: str
    primer_pair: str = "Unknown"   # specific forward/reverse pair identifier 
    assay_name: str = "Unknown"    # gene/locus-level grouping across multiple pairs  
    alignment: Optional[str] = field(default=None)
    query_sequence: Optional[str] = field(default=None)
    subject_sequence: Optional[str] = field(default=None)
    theta: Optional[float] = field(default=None)
    dG: Optional[float] = field(default=None)
    three_prime_window: Optional[str] = field(default=None)
    three_prime_stable: Optional[bool] = field(default=None)
    three_prime_mismatches: Optional[int] = field(default=None)

@dataclass(slots=True)
class PairRisk:
    """
    Structured evaluation format for verified convergent structural pairings.
    
    This class aggregates paired forward and reverse off-target alignment data
    to assess the probability of physical, spurious amplification. A risk is 
    categorized based on geometric constraints (distance, orientation) and 3' stability, 
    and further prioritized by thermodynamic fractional occupancy (theta).
    """
    primer_pair: str
    assay_name: str
    accession: str
    fwd_strand: str
    rev_strand: str
    fwd_back_start: int
    fwd_back_end: int
    rev_back_start: int
    rev_back_end: int
    amplicon_size_bp: int
    fwd_pct_identity: float
    rev_pct_identity: float
    fwd_theta: Optional[float]
    rev_theta: Optional[float]
    fwd_dG: Optional[float]
    rev_dG: Optional[float]
    fwd_3prime_window: Optional[str]
    rev_3prime_window: Optional[str]
    fwd_3prime_mismatches: Optional[int]
    rev_3prime_mismatches: Optional[int]
    convergent_pair: bool
    amplification_risk: bool
    stable_amplification_risk: Optional[bool]
    risk_classification: str
    risk_basis: str
    intersected_gene: str = "Not Evaluated"

# Convert thermodynamics lookup array into operational dataframe in memory
THERMO_DF = pd.DataFrame(NEAREST_NEIGHBOUR_DATA, columns=['Pair One', 'Pair Two', 'Delta H', 'Delta H Uncertainty', 'Delta S', 'Delta S Uncertainty'])
THERMO_DF['Gibbs free energy'] = THERMO_DF['Delta H'] - (310.15 * (THERMO_DF['Delta S'] / 1000.0))

DNA_COMP = str.maketrans({"A": "T", "T": "A", "C": "G", "G": "C"})

def reverse_complement(seq: str) -> str:
    return seq.upper().translate(DNA_COMP)[::-1]

# Structural GFF terms to ignore so they are not mistaken for real gene names
_GENERIC_ATTRIBUTE_VALUES = frozenset({
        "gene", "cds", "mrna", "exon", "region",
        "sequence_feature", "biological_region", "primary_transcript",
        "lnc_rna", "transcript",
 })

def lookup_intersecting_gene(
    accession: str,
    hit_start: int,
    hit_end: int,
    api_key: str = None,
) -> str:
    """Identify all genes or genomic features at a coordinate window via NCBI.

    Queries the E-utilities efetch endpoint for the given accession and
    coordinate range, then parses the returned GenBank flatfile for all
    informative feature names overlapping that window.

    Where multiple distinct genes are present (e.g. overlapping genes, operon
    junctions, or multi-gene amplicons), all are returned joined by ' | ' in
    chromosomal order.
    """
    time.sleep(0.15 if api_key else 0.5)

    params = {
        "db": "nuccore",
        "id": accession,
        "seq_start": min(hit_start, hit_end),
        "seq_stop": max(hit_start, hit_end),
        "rettype": "gb",
        "retmode": "text",
    }
    if api_key:
        params["api_key"] = api_key

    try:
        res = ncbi_session.get(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
            params=params,
            timeout=20.0,
        )

        if res.status_code == 429:
            return "Lookup Error (Rate Limited)"
        if res.status_code != 200:
            return f"Lookup Error (HTTP {res.status_code})"

        text = res.text

        # NCBI returns HTTP 200 with a space-padded ASCII error string for bad
        # accessions. Collapse whitespace to detect this reliably.
        collapsed = text.replace(" ", "").replace("\n", "")
        if collapsed.startswith("Error:") or "Failedtounderstand" in collapsed:
            return "Lookup Error (Bad Accession)"

        # Priority 1: collect all distinct /gene= symbols in flatfile order.
        seen: set[str] = set()
        genes: list[str] = []
        for match in re.finditer(r'/gene="([^"]+)"', text):
            val = match.group(1).strip()
            if val.lower() not in _GENERIC_ATTRIBUTE_VALUES and val not in seen:
                seen.add(val)
                genes.append(val)
        if genes:
            return " | ".join(genes)

        # Priority 2: collect all distinct /product= descriptions in flatfile order.
        seen.clear()
        products: list[str] = []
        for match in re.finditer(r'/product="([^"]+)"', text):
            val = match.group(1).strip()
            if val not in seen:
                seen.add(val)
                products.append(val)
        if products:
            return " | ".join(products)

        # Priority 3: first /locus_tag= as last resort (always present in RefSeq).
        locus = re.search(r'/locus_tag="([^"]+)"', text)
        if locus:
            return locus.group(1).strip()

        return "Intergenic Region"

    except requests.exceptions.Timeout:
        return "Lookup Error (Timeout)"
    except Exception as e:
        import traceback
        traceback.print_exc()
        return "Lookup Error (E-utils)"
        
def calculate_thermodynamics(plus_sense: str, minus_sense: str, conditions: AssayConditions) -> Optional[dict]:
    """Computes nearest-neighbour thermodynamic parameters for a primer-template duplex."""
    plus_sense = plus_sense.upper()
    minus_sense = minus_sense.upper()
    output = {"Enthalpy": 0.2, "Entropy": -5.7, "Gibbs": 1.96}

    if len(plus_sense) < 2 or len(minus_sense) < 2:
        return None

    nn_pairs = zip(
        [plus_sense[i : i + 2] for i in range(len(plus_sense) - 1)],
        [minus_sense[i : i + 2] for i in range(len(minus_sense) - 1)],
    )

    for p1, p2 in nn_pairs:
        match = THERMO_DF[(THERMO_DF["Pair One"] == p1) & (THERMO_DF["Pair Two"] == p2)]
        if match.empty:
            match = THERMO_DF[(THERMO_DF["Pair One"] == p2[::-1]) & (THERMO_DF["Pair Two"] == p1[::-1])]
        if match.empty:
            return None

        row = match.iloc[0]
        output["Enthalpy"] += float(row["Delta H"])
        output["Entropy"] += float(row["Delta S"])
        output["Gibbs"] += float(row["Gibbs free energy"])

    for terminal in [plus_sense[0] + minus_sense[0], plus_sense[-1] + minus_sense[-1]]:
        if terminal in {"AT", "TA"}:
            output["Enthalpy"] += TERMINAL_CORRECTIONS["terminal_at_enthalpy"]
            output["Entropy"] += TERMINAL_CORRECTIONS["terminal_at_entropy"]
        elif terminal in {"CG", "GC"}:
            output["Enthalpy"] += TERMINAL_CORRECTIONS["terminal_cg_enthalpy"]
            output["Entropy"] += TERMINAL_CORRECTIONS["terminal_cg_entropy"]

    if plus_sense == reverse_complement(plus_sense) and minus_sense == plus_sense:
        output["Entropy"] += TERMINAL_CORRECTIONS["symmetry_correction_entropy"]
 
    if conditions.divalent_salt_conc > conditions.dntp_conc:
        salt_conc = (conditions.monovalent_salt_conc + 120.0 * math.sqrt(conditions.divalent_salt_conc - conditions.dntp_conc)) * 1e-3
    else:
        salt_conc = conditions.monovalent_salt_conc * 1e-3

    phosphates_present = (len(plus_sense) * 2) - 2
    output["Adjusted Entropy"] = output["Entropy"] + (0.368 * (phosphates_present / 2.0) * math.log(salt_conc))
    output["Adjusted Gibbs"] = output["Enthalpy"] - (310.15 * output["Adjusted Entropy"] / 1000.0)
    return output

def calculate_equilibrium_constant(dH_cal: float, dS_cal: float, temperature_c: float) -> float:
    R = 1.9872  # cal/(mol*K)
    T_kelvin = temperature_c + 273.15
    dG = dH_cal - (T_kelvin * dS_cal)
    return float(np.exp(-dG / (R * T_kelvin)))

def solve_fraction_bound(K, C_A, C_B):
    if K < 1e-10 or C_A <= 0 or C_B <= 0: 
        return 0.0 # Heuristic override for extremely unstable pairs 
    
    a = K
    b = -(K * C_A + K * C_B + 1)
    c = K * C_A * C_B
    
    # Prevents math domain errors from floating point precision limits
    discriminant = max(0.0, b**2 - 4*a*c)
    x = (-b - np.sqrt(discriminant)) / (2*a)
    
    # Calculate fraction and clip strictly to [0.0, 1.0] to protect against float jitter
    theta = x / C_A
    return float(max(0.0, min(1.0, theta)))

def parse_alignment_window(hit: BlastHit, n_bases: int, max_3prime_mismatches: int) -> None:
    """
    Isolates the 3' terminal sequence window from the BLAST alignment block.
    Maps lines correctly: Line 0 is Subject (Template), Line 2 is Query (Primer).

    The 3' terminal tip of a primer must bind stably to the template sequence for 
    the polymerase enzyme to initiate extension. If the 3' end is frayed due to too 
    many mismatches, amplification will fail, even if the rest of the primer matches.
    """
    if not hit.alignment:
        return
        
    lines = [line.strip() for line in hit.alignment.split("\n") if line.strip()]
    if len(lines) < 3:
        return

    # CORRECT BIOLOGICAL MAPPING (Verified by static query sequence on Line 2)
    subject_sequence = lines[0]  # Background Template
    match_string = lines[1]      # Midline pipes/dashes
    query_sequence = lines[2]    # Primer Query

    # Truncate match midline to align perfectly with the Primer Query length
    clean_match = match_string[:len(query_sequence)]
    
    # Isolate the final 3' bases of the primer matching region
    window = clean_match[-n_bases:] if len(clean_match) >= n_bases else clean_match

    # Calculate mismatches within the 3' window
    mismatches = len(window) - window.count("|")
    
    # Assign properties back to the BlastHit container object cleanly
    hit.query_sequence = query_sequence
    hit.subject_sequence = subject_sequence
    hit.three_prime_window = window
    hit.three_prime_mismatches = mismatches
    
    # A hit is an amplification risk if its 3' end is stable (mismatches under threshold)
    hit.three_prime_stable = mismatches <= max_3prime_mismatches
    
def load_blast_csv(
    filepath: str | Path,
    primer_direction: str,
    alignment_col: str,
    conditions: AssayConditions,
    operating_temp_c: float,
    n_bases: int,
    max_3prime_mismatches: int,
) -> dict[str, list[BlastHit]]:
    """Loads a single primer's BLAST output matrix and derives thermodynamic attributes,
    based on the sequence of the off-target duplex."""
    filepath = Path(filepath)
    df = pd.read_csv(filepath)
    hits_by_accession = defaultdict(list)
    c_template, c_primer = conditions.template_conc * 1e-6, conditions.oligo_conc * 1e-6

    if alignment_col not in df.columns:
        logger.warning(f"Target alignment column '{alignment_col}' not found in {filepath.name}. Available: {list(df.columns)}")

    # Locate tracking columns case-insensitively
    primer_pair_col = next((col for col in df.columns if col.lower() == "primerpair"), None)
    assay_name_col  = next((col for col in df.columns if col.lower() in {"assayname", "assay_name"}), None)

    for _, row in df.iterrows():
        # Extract and clean the alignment string from the CSV row
        raw_align = row.get(alignment_col)
        alignment = raw_align if isinstance(raw_align, str) and pd.notna(raw_align) else None

        if alignment:
            lines = alignment.splitlines()
            if lines and any(c.isdigit() for c in lines[0][:5]):
                import re
                lines[0] = re.sub(r'^\d+\s*', '', lines[0])
                alignment = "\n".join(lines)
    
        primer_pair = str(row[primer_pair_col]) if primer_pair_col and pd.notna(row[primer_pair_col]) else filepath.stem
        
        if assay_name_col and pd.notna(row[assay_name_col]):
            assay_name = str(row[assay_name_col])
        else:
            assay_name = _infer_assay_name(primer_pair, primer_pair)

        hit = BlastHit(
            seq_id=str(row["Background_SeqID"]), primer_direction=primer_direction,
            pct_identity=float(row["Percentage Identity"]), q_start=int(row["qStart"]), q_end=int(row["qEnd"]),
            back_seq_start=int(row["BackSeq_Start"]), back_seq_end=int(row["BackSeq_End"]),
            evalue=float(row["Evalue"]), num_gaps=int(row["Number Of Gaps"]), strand=str(row["Strand"]).strip().lower(),
            primer_pair=primer_pair, assay_name=assay_name, alignment=alignment,
        )

        parse_alignment_window(hit, n_bases, max_3prime_mismatches)

        if hit.query_sequence and hit.subject_sequence:
            thermo = calculate_thermodynamics(hit.query_sequence, hit.subject_sequence, conditions)
            if thermo is not None:
                hit.dG = thermo["Adjusted Gibbs"]
                dH_cal = thermo['Enthalpy'] * 1000.0
                dS_cal = thermo['Adjusted Entropy']
                
                K = calculate_equilibrium_constant(dH_cal, dS_cal, operating_temp_c)
                hit.theta = solve_fraction_bound(K, c_template, c_primer)

        hits_by_accession[hit.seq_id].append(hit)
    return dict(hits_by_accession)

def check_amplification_risk(
    fwd_hit: BlastHit, rev_hit: BlastHit, max_amplicon_size: int, theta_threshold: float, max_3prime_mismatches: int, api_key: Optional[str] = None
) -> Optional[PairRisk]:
    """
    Evaluates spatial/geometric constraints - convergent orientiation, amplicon size -
    and thermodynamic stability to categorise primer-pair off-target amplification risks.
    """
    convergent = False
    fwd_3p, rev_3p = fwd_hit.back_seq_end, rev_hit.back_seq_end

    if fwd_hit.strand == "plus" and rev_hit.strand == "minus" and fwd_3p < rev_3p:
        convergent = True
    elif fwd_hit.strand == "minus" and rev_hit.strand == "plus" and rev_3p < fwd_3p:
        convergent = True

    if not convergent:
        return None

    amp_left = min(fwd_hit.back_seq_start, fwd_hit.back_seq_end, rev_hit.back_seq_start, rev_hit.back_seq_end)
    amp_right = max(fwd_hit.back_seq_start, fwd_hit.back_seq_end, rev_hit.back_seq_start, rev_hit.back_seq_end)
    amplicon_size = amp_right - amp_left + 1
    
    if amplicon_size > max_amplicon_size:
        return None

    amplification_risk = bool(fwd_hit.three_prime_stable) and bool(rev_hit.three_prime_stable)
    fwd_t, rev_t = fwd_hit.theta or 0.0, rev_hit.theta or 0.0

    stable_risk = amplification_risk and fwd_t >= theta_threshold and rev_t >= theta_threshold

    # Elevated top-down gate: Tier 1 checks for absolute zero flawless matches first
    if (
        stable_risk 
        and fwd_hit.three_prime_mismatches <= 2 
        and rev_hit.three_prime_mismatches <= 2
    ):
        classification = "Tier 1: Severe Cross-Reactivity Risk"
        gene_symbol = lookup_intersecting_gene(accession=fwd_hit.seq_id, hit_start=amp_left, hit_end=amp_right, api_key=api_key)
    elif stable_risk:
        classification = "Tier 2: Plausible Amplification Risk"
        gene_symbol = "Skipped (Low Structural Threat Profile)"
    else:
        # Non-risk items return None directly so they never reach the output dataframe
        return None

    return PairRisk(
        primer_pair=fwd_hit.primer_pair, assay_name=fwd_hit.assay_name, accession=fwd_hit.seq_id, fwd_strand=fwd_hit.strand, 
        rev_strand=rev_hit.strand,fwd_back_start=fwd_hit.back_seq_start, fwd_back_end=fwd_hit.back_seq_end,
        rev_back_start=rev_hit.back_seq_start, rev_back_end=rev_hit.back_seq_end,
        amplicon_size_bp=amplicon_size, fwd_pct_identity=fwd_hit.pct_identity, rev_pct_identity=rev_hit.pct_identity,
        fwd_theta=fwd_hit.theta, rev_theta=rev_hit.theta, fwd_dG=fwd_hit.dG, rev_dG=rev_hit.dG,
        fwd_3prime_window=fwd_hit.three_prime_window, rev_3prime_window=rev_hit.three_prime_window,
        fwd_3prime_mismatches=fwd_hit.three_prime_mismatches, rev_3prime_mismatches=rev_hit.three_prime_mismatches,
        convergent_pair=True, amplification_risk=amplification_risk, stable_amplification_risk=stable_risk,
        risk_classification=classification, risk_basis="structural_and_thermodynamic", intersected_gene=gene_symbol
    )

def evaluate_primer_sequestration_risks(hits: dict[str, list[BlastHit]], theta_threshold: float) -> pd.DataFrame:
    """
    Identifies stable off-target monomer binding events that risk depleting free primers.
    
    This function scans individual primer alignments to capture high-occupancy 
    (high-theta) single-hybridisation events. 
    
    While these isolated alignments may not satisfy the geometric constraints required to generate a spurious amplicon,
    an excess of energetically stable primer-template interactions exhausts the available primer pool for the intended target, 
    thereby reducing overall assay sensitivity.
    
    Parameters
    ----------
    hits : dict[str, list[BlastHit]]
        A dictionary mapping background sequence accessions to lists of their 
        respective BLAST alignment metrics.
    theta_threshold : float
        The minimum fractional occupancy value (0.0 to 1.0) required to 
        classify an alignment as a significant sequestration hazard.
        
    Returns
    -------
    pd.DataFrame
        A compiled summary table containing accessions, structural coordinates, 
        and thermodynamic properties of individual primers flagged for consumption risk.
    """
    rows = []
    for accession, hit_list in hits.items():
        for hit in hit_list:
            if hit.theta is not None and hit.theta >= theta_threshold:
                rows.append({
                    "primer_pair": hit.primer_pair, "accession": accession, "primer_direction": hit.primer_direction,
                    "primer": hit.query_sequence, "template": hit.subject_sequence,
                    "back_seq_start": hit.back_seq_start, "back_seq_end": hit.back_seq_end,
                    "strand": hit.strand, "theta": hit.theta, "dG": hit.dG, "pct_identity": hit.pct_identity,
                    "three_prime_window": hit.three_prime_window, "three_prime_mismatches": hit.three_prime_mismatches,
                    "sequestration_risk": True,
                })
    return pd.DataFrame(rows)

def plot_specificity_landscape(pairwise_df: pd.DataFrame, output_png: str | Path) -> None:
    """
    Generates a modern, publication-quality scatter plot using Seaborn mapping 
    Predicted Amplicon Size vs. Maximum Primer Occupancy.
    """ 
    output_png = Path(output_png)
    output_png.parent.mkdir(parents=True, exist_ok=True)
    
    # 1. Set clean visual theme and context upfront
    sns.set_theme(style="whitegrid", palette="muted")
    fig, ax = plt.subplots(figsize=(9, 6))

    # Handle empty data gracefully
    if pairwise_df.empty or "risk_classification" not in pairwise_df.columns:
        ax.text(0.5, 0.5, "No Pairwise Amplification Risks Identified", 
                ha="center", va="center", transform=ax.transAxes, color="gray", fontsize=12)
        ax.set_xlim(0, 500)
        ax.set_ylim(-0.05, 1.05)
        fig.savefig(output_png, dpi=300, bbox_inches="tight")
        plt.close(fig)
        return

    # 2. Dynamically calculate the true maximum occupancy risk if individual columns exist
    # Fallback to a single 'theta' column if that's what the dataclass uses
    if "fwd_theta" in pairwise_df.columns and "rev_theta" in pairwise_df.columns:
        pairwise_df["max_theta_risk"] = pairwise_df[["fwd_theta", "rev_theta"]].max(axis=1)
    elif "theta" in pairwise_df.columns:
        pairwise_df["max_theta_risk"] = pairwise_df["theta"]
    else:
        pairwise_df["max_theta_risk"] = 1.0  # Absolute fallback indicator

    # 3. Custom color palette mapping directly to your pipeline's Tiers
    tier_colors = {
        "Tier 1: Severe Cross-Reactivity Risk": "#D62728",    # Crisp Crimson
        "Tier 2: Plausible Amplification Risk": "#FF7F0E"     # Modern Orange
    }

    # 4. Generate the plot using standard Seaborn objects
    sns.scatterplot(
        data=pairwise_df,
        x="amplicon_size_bp",
        y="max_theta_risk",
        hue="risk_classification",
        style="risk_classification",
        palette=tier_colors,
        s=100,              # Increased marker size for scannability
        alpha=0.85,         # Slight transparency for overlapping points
        edgecolor="w",      # Crisp white borders around markers
        linewidth=1.2,
        ax=ax
    )

    # 5. Clean layout aesthetics
    ax.set_title("RPA Assay Cross-Amplification Risk Profile", fontsize=14, fontweight="bold", pad=15)
    ax.set_xlabel("Predicted Off-Target Amplicon Size (bp)", fontsize=11, fontweight="semibold", labelpad=8)
    ax.set_ylabel(r"Maximum Fractional Occupancy ($\max(\\theta_{fwd}, \\theta_{rev})$)", fontsize=11, fontweight="semibold", labelpad=8)
    ax.set_ylabel(r"Maximum Fractional Occupancy ($\max(\theta_{fwd}, \theta_{rev})$)", fontsize=11, fontweight="semibold", labelpad=8)
    
    # Pad layout limits based on maximum detected amplicon size
    max_size = pairwise_df["amplicon_size_bp"].max()
    ax.set_xlim(0, max(max_size + 50, 500))
    ax.set_ylim(-0.05, 1.05)

    # Clean up legend placement
    ax.legend(title="Assay Risk Tiers", title_fontsize="10", loc="upper right", frameon=True, shadow=False)
    sns.despine(left=True, bottom=True) # Removes ugly black outer border box

    # 6. Save high-resolution file
    fig.tight_layout()
    fig.savefig(output_png, dpi=300, bbox_inches="tight")
    plt.close(fig)
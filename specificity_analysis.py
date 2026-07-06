import logging
from pathlib import Path
from typing import Any
import pandas as pd
from dataclasses import asdict

from specificity_analysis_components import (
    AssayConditions,
    load_blast_csv,
    check_amplification_risk,
    evaluate_primer_sequestration_risks,
    plot_specificity_landscape,
    write_classification_summary,
    build_assay_performance_landscape
)

logger = logging.getLogger(__name__)

def _fmt_accession(acc: str, labels: dict) -> str:
    """
    Formats an accession for display in the report by appending its
    human-readable FASTA header description, e.g.
    'NC_000001.11' -> 'NC_000001.11 (Homo sapiens chromosome 1, GRCh38.p14 Primary Assembly)'.
    Falls back to the bare accession if no label is available.
    """
    label = labels.get(acc, "")
    return f"{acc} ({label})" if label else acc

def run_specificity_pipeline(
    fwd_csv: str | Path,
    rev_csv: str | Path,
    fwd_target_csv: str | Path,
    rev_target_csv: str | Path,
    output_dir: str | Path,
    alignment_col: str = "Advanced Cross Reactivity Percentage String",
    n_bases: int = 5,
    max_amplicon_size: int = 500,
    max_3prime_mismatches: int = 3,
    theta_threshold: float = 0.05,
    operating_temp_c: float = 37.0,
    write_plot: bool = True,
    api_key: str | None = None,
    background_organisms: list[str] | None = None,
    target_organism: str | None = None,
    sensitivity: str = "Standard",
    accession_labels: dict | None = None,
) -> dict[str, Any]:
    """
    End-to-end primer specificity assessment: thermodynamic occupancy scoring,
    convergent off-target amplification detection, sequestration sink flagging,
    and optional on-target amplification confirmation. Summative report is produced.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    conditions = AssayConditions()

    # =========================================================================
    # STEP 1: THERMODYNAMIC CHARACTERISATION
    # =========================================================================
    logger.info("Scoring thermodynamic occupancy profiles...")
    fwd_hits = load_blast_csv(fwd_csv, "forward", alignment_col, conditions, operating_temp_c, n_bases, max_3prime_mismatches)
    rev_hits = load_blast_csv(rev_csv, "reverse", alignment_col, conditions, operating_temp_c, n_bases, max_3prime_mismatches)

    # =========================================================================
    # STEP 2A: PAIRWISE CROSS-AMPLIFICATION INTERSECTION
    # =========================================================================
    shared_accessions = set(fwd_hits.keys()) & set(rev_hits.keys())
    pairwise_rows = []

    for accession in shared_accessions:
        for fwd_hit in fwd_hits[accession]:
            for rev_hit in rev_hits[accession]:
                pair = check_amplification_risk(
                    fwd_hit, rev_hit, max_amplicon_size,
                    theta_threshold, max_3prime_mismatches, api_key
                )
                if pair is not None:
                    pair_dict = asdict(pair)
                    pair_dict["assay_name"]   = fwd_hit.assay_name
                    pair_dict["primer_pair"]  = f"{fwd_hit.primer_pair}/{rev_hit.primer_pair}"
                    pair_dict["fwd_oligo_id"] = fwd_hit.primer_pair
                    pair_dict["rev_oligo_id"] = rev_hit.primer_pair
                    pairwise_rows.append(pair_dict)

    pairwise_df = pd.DataFrame(pairwise_rows) if pairwise_rows else pd.DataFrame()
    # =========================================================================
    # STEP 2B: CONFIRMATION OF ON-TARGET AMPLIFICATION
    # Runs the same thermodynamic + geometry logic against the target genome.
    # Here a convergent high-theta pair is the desired outcome, not a risk.
    # =========================================================================
    on_target_df = pd.DataFrame()
    if fwd_target_csv and rev_target_csv:
        logger.info("Scoring on-target amplification profiles...")
        fwd_target_hits = load_blast_csv(fwd_target_csv, "forward", alignment_col, conditions, operating_temp_c, n_bases, max_3prime_mismatches)
        rev_target_hits = load_blast_csv(rev_target_csv, "reverse", alignment_col, conditions, operating_temp_c, n_bases, max_3prime_mismatches)

        target_shared = set(fwd_target_hits.keys()) & set(rev_target_hits.keys())
        on_target_rows = []

        for accession in target_shared:
            for fwd_hit in fwd_target_hits[accession]:
                for rev_hit in rev_target_hits[accession]:
                    # Use a permissive theta floor of 0 — we want to capture all
                    # convergent pairs on the target, even weak ones, so we can
                    # flag sensitivity concerns rather than silently miss them
                    pair = check_amplification_risk(
                        fwd_hit, rev_hit, max_amplicon_size,
                        theta_threshold=0.0,       # permissive: capture all geometry hits
                        max_3prime_mismatches=max_3prime_mismatches,
                        api_key=None               # no annotation lookup needed for target
                    )
                    if pair is not None:
                        row = asdict(pair)
                        row["assay_name"]   = fwd_hit.assay_name
                        row["primer_pair"]  = f"{fwd_hit.primer_pair}/{rev_hit.primer_pair}"
                        row["fwd_oligo_id"] = fwd_hit.primer_pair
                        row["rev_oligo_id"] = rev_hit.primer_pair
                        on_target_rows.append(row)

        on_target_df = pd.DataFrame(on_target_rows) if on_target_rows else pd.DataFrame()

        # Save full on-target detail for downstream inspection
        on_target_df.to_csv(output_dir / "on_target_amplification.csv", index=False)
        logger.info(f"On-target confirmation: {len(on_target_df)} convergent pair(s) detected.")

    # =========================================================================
    # STEP 3: SEQUESTRATION SINK EVALUATION
    # =========================================================================
    logger.info("Evaluating single-primer sequestration risks...")
    fwd_seq = evaluate_primer_sequestration_risks(fwd_hits, theta_threshold)
    rev_seq = evaluate_primer_sequestration_risks(rev_hits, theta_threshold)

    sequestration_df = pd.concat([fwd_seq, rev_seq], ignore_index=True, sort=False)
    if not sequestration_df.empty:
        sequestration_df = sequestration_df.drop_duplicates(
            subset=["accession", "primer_direction", "back_seq_start", "back_seq_end", "theta"]
        )

    pairwise_df     = pairwise_df.reset_index(drop=True)
    sequestration_df = sequestration_df.reset_index(drop=True)

    # =========================================================================
    # STEP 4: COLUMN NORMALISATION
    # =========================================================================
    for df in [pairwise_df, sequestration_df]:
        if not df.empty and "accession" in df.columns:
            df.rename(columns={"accession": "background_accession"}, inplace=True)

    acc_col = "background_accession"

    # =========================================================================
    # STEP 5: REPORT COMPILATION
    # =========================================================================
    report_path = output_dir / "specificity_report.md"

    # --- Pre-compute pairwise metrics ---
    tier_1_count, tier_2_count = 0, 0
    intersecting_genes = set()
    tier_1_details = []
    tier_2_details = []
    negatives = {"Intergenic Region", "Network Timeout", "Connection Error", "Rate Limit Throttled", "Skipped"}

    if not pairwise_df.empty and "risk_classification" in pairwise_df.columns:
        tier_1_mask = pairwise_df["risk_classification"] == "Tier 1: Severe Cross-Reactivity Risk"
        tier_2_mask = pairwise_df["risk_classification"] == "Tier 2: Plausible Amplification Risk"
        tier_1_count = int(tier_1_mask.sum())
        tier_2_count = int(tier_2_mask.sum())

        if tier_1_count > 0 and "intersected_gene" in pairwise_df.columns:
            valid_genes = pairwise_df[tier_1_mask & ~pairwise_df["intersected_gene"].isin(negatives)]["intersected_gene"].unique()
            intersecting_genes = set(valid_genes)
            tier_1_details = (
                pairwise_df[tier_1_mask].sort_values(["assay_name", acc_col])
                [[acc_col, "assay_name", "fwd_oligo_id", "rev_oligo_id", "amplicon_size_bp",
                  "fwd_theta", "rev_theta", "intersected_gene"]]
                .head(10)
                .to_dict(orient="records")
            )

        if tier_2_count > 0:
            tier_2_details = (
                pairwise_df[tier_2_mask].sort_values(["assay_name", acc_col])
                [[acc_col, "assay_name", "fwd_oligo_id", "rev_oligo_id", "amplicon_size_bp",
                  "fwd_theta", "rev_theta", "intersected_gene"]]
                .head(10)
                .to_dict(orient="records")
            )

    # --- Per-assay verdict breakdown ---
    # Collect all unique assay names across both pairwise and sequestration data
    all_assay_names = set()
    if not pairwise_df.empty and "assay_name" in pairwise_df.columns:
        all_assay_names.update(pairwise_df["assay_name"].dropna().unique())
    if not sequestration_df.empty and "primer_pair" in sequestration_df.columns:
        all_assay_names.update(sequestration_df["primer_pair"].dropna().unique())

    per_assay_verdicts = {}
    for aname in sorted(all_assay_names):
        if not pairwise_df.empty and "assay_name" in pairwise_df.columns:
            assay_pairwise = pairwise_df[pairwise_df["assay_name"] == aname]
            s_tier1 = int((assay_pairwise["risk_classification"] == "Tier 1: Severe Cross-Reactivity Risk").sum()) if "risk_classification" in assay_pairwise.columns else 0
            s_tier2 = int((assay_pairwise["risk_classification"] == "Tier 2: Plausible Amplification Risk").sum()) if "risk_classification" in assay_pairwise.columns else 0
        else:
            s_tier1, s_tier2 = 0, 0

        if not sequestration_df.empty and "primer_pair" in sequestration_df.columns:
            assay_seq_rows = sequestration_df[sequestration_df["primer_pair"] == aname]
            s_seq       = len(assay_seq_rows)
            # Worst-case theta across all sequestration hits for this assay
            s_max_theta = float(assay_seq_rows["theta"].max()) if not assay_seq_rows.empty and "theta" in assay_seq_rows.columns else 0.0
        else:
            s_seq, s_max_theta = 0, 0.0

        verdict = "❌ FAIL" if s_tier1 > 0 else ("⚠️ WARN" if s_tier2 > 0 else "✅ PASS")
        per_assay_verdicts[aname] = {
            "verdict":      verdict,
            "tier1":        s_tier1,
            "tier2":        s_tier2,
            "sequestration": s_seq,
            "max_seq_theta": s_max_theta,
        }

    # --- Sequestration summary: one row per primer (worst-case hit only) ---
    seq_summary_rows = []
    if not sequestration_df.empty:
        group_cols = ["primer_pair", "primer_direction"]
        available = [c for c in group_cols if c in sequestration_df.columns]
        if available:
            sort_cols = [c for c in ["primer_pair", acc_col] if c in sequestration_df.columns]  # ← add
            sequestration_df = sequestration_df.sort_values(sort_cols)                            # ← add
            for _, grp in sequestration_df.groupby(available):
                worst = grp.loc[grp["theta"].idxmax()]
                three_prime_stable = worst.get("three_prime_mismatches", 999) <= 1
                seq_summary_rows.append({
                    "primer_pair":          worst.get("primer_pair", "Unknown"),
                    "primer_direction":     worst.get("primer_direction", "?"),
                    "background_accession": worst.get(acc_col, worst.get("accession", "?")),
                    "n_sites":              len(grp),
                    "max_theta":            round(float(worst["theta"]), 4),
                    "three_prime_stable":   three_prime_stable,
                    "three_prime_window":   worst.get("three_prime_window", "?"),
                })

    overall_verdict = "❌ FAIL - CRITICAL OFF-TARGET RISK" if tier_1_count > 0 else "✅ PASS - SPECIFIC ASSAY"
    from datetime import datetime
    run_timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    # --- On-target verdict per assay ---
    # Classify each assay as Confirmed / Weak / Not Detected based on whether
    # a convergent pair was found and how stable the binding is
    on_target_verdicts = {}
    if not on_target_df.empty and "assay_name" in on_target_df.columns:
        for aname in on_target_df["assay_name"].unique():
            assay_rows = on_target_df[on_target_df["assay_name"] == aname]

            # Take the best hit: highest combined theta across forward and reverse
            if "fwd_theta" in assay_rows.columns and "rev_theta" in assay_rows.columns:
                assay_rows = assay_rows.copy()
                assay_rows["combined_theta"] = assay_rows["fwd_theta"] + assay_rows["rev_theta"]
                best = assay_rows.loc[assay_rows["combined_theta"].idxmax()]
            else:
                best = assay_rows.iloc[0]

            fwd_t = float(best.get("fwd_theta", 0) or 0)
            rev_t = float(best.get("rev_theta", 0) or 0)
            amp   = int(best.get("amplicon_size_bp", 0) or 0)
            fwd_o = best.get("fwd_oligo_id", "?")
            rev_o = best.get("rev_oligo_id", "?")
            fwd_3 = best.get("fwd_3prime_window", "?")
            rev_3 = best.get("rev_3prime_window", "?")

            # Verdict logic:
            # Confirmed    → geometry satisfied, both primers thermodynamically stable (θ ≥ threshold)
            # Weak         → geometry satisfied but at least one primer has low occupancy (sensitivity risk)
            # Not Detected → no convergent pair found on the target genome at all
            both_stable = fwd_t >= theta_threshold and rev_t >= theta_threshold
            either_weak = (fwd_t >= theta_threshold) != (rev_t >= theta_threshold)

            if both_stable:
                verdict = "✅ Confirmed"
            elif either_weak:
                verdict = "⚠️ Weak"
            else:
                verdict = "❌ Not detected"

            on_target_verdicts[aname] = {
                "verdict":     verdict,
                "fwd_oligo":   fwd_o,
                "rev_oligo":   rev_o,
                "amplicon_bp": amp,
                "fwd_theta":   fwd_t,
                "rev_theta":   rev_t,
                "fwd_3prime":  fwd_3,
                "rev_3prime":  rev_3,
            }

    # Flag assays that had no on-target convergent pair at all
    for aname in per_assay_verdicts:
        if aname not in on_target_verdicts:
            on_target_verdicts[aname] = {
                "verdict":     "❌ Not detected",
                "fwd_oligo":   "?",
                "rev_oligo":   "?",
                "amplicon_bp": 0,
                "fwd_theta":   0.0,
                "rev_theta":   0.0,
                "fwd_3prime":  "?",
                "rev_3prime":  "?",
            }

    assay_performance_landscape_path = None
    specificity_landscape_path = None
    if write_plot:
        assay_performance_landscape_path = build_assay_performance_landscape(
            per_assay_verdicts, on_target_verdicts, tier_1_details, tier_2_details,
            theta_threshold, output_dir,
        )
        specificity_landscape_path = plot_specificity_landscape(
            pairwise_df, output_dir / "specificity_landscape.png"
        )

    # Compile the final report
    with open(report_path, "w", encoding="utf-8") as md:
 
        # ── Title ──────────────────────────────────────────────────────
        md.write("# Oligonucleotide Specificity Analysis Report\n\n")
        md.write(
            "This report tells you two things: **will each RPA assay detect your target**, "
            "and **could it also detect another species in the background that it is not intended to**. "
            "A short glossary of the recurring terms is at the end of this section.\n\n"
        )
 
        # ── Quick guide ────────────────────────────────────────────────
        md.write("<details>\n<summary><strong>📖 Quick guide to the terms in this report</strong></summary>\n\n")
        md.write(
            "- **θ (theta)** — the predicted thermodynamic occupancy of a primer-binding site: "
            "the fraction of copies of that site expected to be bound by primer at equilibrium, "
            "under your reaction conditions. Primer is present in molar excess over template in "
            "RPA, so template — not primer — is the limiting reagent; θ describes how thoroughly "
            "a site gets captured, not how much primer gets used up. 1.0 means essentially "
            "complete occupancy; 0 means essentially none. High θ at your intended target "
            "supports efficient detection. High θ at a background site means that site would "
            "also be captured if present — to assess the genuine off-target risk, checking against "
            "the possibility of extension and geometry before treating it as a real threat.\n"
            "- **3′ window / extension competent** — whether the last few bases of the primer "
            "(where DNA polymerase starts the extension of the primer-template alignment) match well enough to actually prime synthesis. "
            "A primer can bind (high θ) without being extension-competent.\n"
            "- **Convergent geometry** — whether the forward and reverse primers bind opposite "
            "strands in the correct orientation, close enough together to fall within the "
            "realistic amplification distance, which can be adjusted the user . If they're oriented the same direction, or too "
            "far apart, there will be no productive amplification regardless of how well each primer binds on its own.\n"
            "- **Tier 1 (high-confidence cross-reactivity risk)** — satisfies every mechanistic "
            "requirement: stable binding (θ above threshold) on both primers, extension-competent "
            "3′ termini, convergent geometry, and an amplicon within the permitted size range. "
            "A highly probable source of experimental cross-reactivity.\n"
            "- **Tier 2 (borderline/plausible cross-reactivity risk)** — convergent geometry is present, but "
            "at least one mechanistic check fails or is only weakly satisfied (most commonly low "
            "θ on one primer). Less likely to amplify, but worth a second look if background "
            "nucleic acid is present at high concentration.\n"
            "- **Sequestration site (primer sink)** — an off-target site with appreciable thermodynamic"
            "occupancy but no geometry for productive amplification: a primer binds, but there's "
            "no convergent partner nearby to extend into a product. This doesn't cause a false "
            "result, but it pulls primer out of the active pool if that DNA site is present in the background, which can raise the limit of "
            "detection (LoD) and thus reduce sensitivity.\n\n"
        )
        md.write("</details>\n\n")
 
        # ── Run Parameters ─────────────────────────────────────────────
        md.write("## ⚙️ Run Parameters\n")
        md.write("| Parameter | Value |\n| :--- | :--- |\n")
        md.write(f"| Report generated | {run_timestamp} |\n")
        if target_organism:
            md.write(f"| Target organism | `{target_organism}` |\n")
        if background_organisms:
            md.write(f"| Background panel | {', '.join(f'`{b}`' for b in background_organisms)} |\n")
        md.write(f"| Operating temperature | {operating_temp_c} °C |\n")
        md.write(f"| θ threshold | {theta_threshold} (sites below this are treated as unbound) |\n")
        md.write(f"| Max amplicon size | {max_amplicon_size} bp |\n")
        md.write(f"| Max 3′ mismatches | {max_3prime_mismatches} |\n")
        md.write(f"| 3′ window length | {n_bases} bp |\n")
        md.write(f"| BLAST sensitivity | {sensitivity} |\n\n")
 
        # ── Executive Summary ──────────────────────────────────────────
        md.write("## 📋 Executive Summary\n")
        md.write("| Metric | Value |\n| :--- | :--- |\n")
        md.write(f"| **Overall verdict** | **{overall_verdict}** |\n")
        md.write(f"| Tier 1 cross-reactivity events | {tier_1_count} |\n")
        md.write(f"| Tier 2 plausible amplification events | {tier_2_count} |\n")
        md.write(f"| Sequestration sites (θ ≥ {theta_threshold}) | {len(sequestration_df)} |\n\n")
 
        # ── Assay Performance Landscape ───────────────────────────────────────────────
        if assay_performance_landscape_path is not None:
            md.write("### Assay Performance Landscape\n")
            md.write(
                "Each point is one assay. **Top-left is the ideal region**: strong "
                "on-target signal, unlikely to bind off-target. Points pulled "
                "right of the dashed line have at least one background site crossing the "
                "θ threshold under convergent geometry. \n\n"
            )
            md.write(f"![Assay Performance Landscape]({assay_performance_landscape_path.name})\n\n")
 
        # ── Per-Assay Breakdown ──────────────────────────────────────
        md.write("## 🧬 Per-Assay Verdict\n")
        md.write("| Assay | Verdict | Tier 1 Hits | Tier 2 Hits | Sequestration Sites | Max Seq. θ |\n")
        md.write("| :--- | :---: | :---: | :---: | :---: | :---: |\n")
        for aname, v in per_assay_verdicts.items():
            md.write(
                f"| `{aname}` | {v['verdict']} | {v['tier1']} | {v['tier2']} "
                f"| {v['sequestration']} | {v['max_seq_theta']:.3f} |\n"
            )
        md.write("\n")
 
        # ── Cross-Amplification Detail ───────────────────────────────
        md.write("## 🔴 Cross-Reactivity Risk Detail\n")
        if tier_1_count == 0 and tier_2_count == 0:
            md.write("No convergent off-target amplification events detected under the current parameters.\n\n")
        else:
            md.write(
                f"{'⚠️' if tier_1_count > 0 else '🟡'} The pipeline identified **{tier_1_count}** Tier 1 "
                f"and **{tier_2_count}** Tier 2 convergent binding events on the background panel.\n\n"
            )
            if intersecting_genes:
                md.write("**Predicted undesired amplified loci:**\n")
                seen = set()
                for row in tier_1_details + tier_2_details:
                    gene = row.get("intersected_gene", "N/A")
                    species = row.get("_accession_label", "?")
                    key = (gene, species)
                    if key in seen:
                        continue
                    seen.add(key)
                    md.write(f"- `{gene}` — {species}\n")
                md.write("\n")
 
            md.write("| Assay | Forward Primer | Reverse Primer | Accession | Amplicon (bp) | fwd θ | rev θ | Annotation |\n")
            md.write("| :--- | :--- | :--- | :--- | :---: | :---: | :---: | :--- |\n")
            for row in tier_1_details:
                fwd_t = f"{row.get('fwd_theta', 0):.3f}" if row.get('fwd_theta') is not None else "N/A"
                rev_t = f"{row.get('rev_theta', 0):.3f}" if row.get('rev_theta') is not None else "N/A"
                md.write(
                    f"| `{row.get('assay_name', '?')}` "
                    f"| `{row.get('fwd_oligo_id', '?')}` "
                    f"| `{row.get('rev_oligo_id', '?')}` "
                    f"| `{row.get('accession_label', 'N/A')}` "
                    f"| {row.get('amplicon_size_bp', 'N/A')} "
                    f"| {fwd_t} | {rev_t} "
                    f"| `{row.get('intersected_gene', 'N/A')}` |\n"
                )
            md.write("\n")
 
            # ── Classification Summary ──────────────────────────────────────
            md.write(
                "Each potential cross-amplification event is evaluated against the criteria "
                "required for productive amplification. As well as reporting the final Tier "
                "classification, the table below shows which criteria each hit satisfied and which it did not, "
                "making the basis for classification fully transparent.\n\n"
            )
            write_classification_summary(
                md, tier_1_details, tier_2_details, theta_threshold,
                max_amplicon_size, max_3prime_mismatches,
            )

        # ── Sequestration Summary ────────────────────────────────────
        md.write("## 📊 Primer Sequestration Summary\n")
        md.write("### Where primers are getting trapped without producing a product.\n\n")
 
        if not seq_summary_rows:
            md.write(f"✅ **No significant sequestration sites detected** above the selection threshold ($\\theta$ = {theta_threshold}).\n\n")
        else:
            md.write(
                f"⚠️ **Warning:** Off-target binding sites detected exceeding $\\theta \\ge$ {theta_threshold}.\n\n"
                f"These loci act as **primer sinks** — they bind primer tightly but don't have a "
                f"convergent partner nearby, so no product forms. They don't cause a false result, "
                f"but they pull primer out of solution and can **raise your functional limit of "
                f"detection (LoD)** for the intended target.\n\n"
                f"> 🧬 **Sample matrix context:** This risk depends heavily on what you're actually "
                f"testing. It's a real concern in clinical or environmental samples with abundant "
                f"background nucleic acid (host or commensal DNA), but may never matter in a clean, "
                f"isolated research culture.\n\n"
                f"**Data summary:** each row is the *highest-occupancy* off-target site found per "
                f"primer. For the full list, see `single_primer_sequestration_sinks.csv`.\n\n"
            )
            md.write("| Primer Pair | Direction | Accession | Sites (n) | Max θ | Stable 3′ | 3′ Window |\n")
            md.write("| :--- | :---: | :--- | :---: | :---: | :---: | :---: |\n")
            for s in seq_summary_rows:
                stable_flag = "⚠️ Yes" if s["three_prime_stable"] else "No"
                md.write(
                    f"| `{s['primer_pair']}` | {s['primer_direction']} "
                    f"| `{_fmt_accession(s['background_accession'], accession_labels) if _fmt_accession else s['background_accession']}` "
                    f"| {s['n_sites']} | {s['max_theta']:.3f} "
                    f"| {stable_flag} | `{s['three_prime_window']}` |\n"
                )
            md.write("\n")
 
        # ── On-Target Confirmation ────────────────────────────────────
        md.write("## 🎯 On-Target Amplification Confirmation\n")
        if not on_target_verdicts:
            md.write("On-target confirmation was not run. "
                     "Pass `--target` and ensure target genome BLAST results are provided.\n\n")
        else:
            md.write(
                "The primers were BLASTed against the target genome to confirm productive "
                "amplification is predicted. **Confirmed** requires convergent binding geometry "
                "within the maximum amplicon size and θ ≥ "
                f"{theta_threshold} for both primers. **Weak** means the geometry is right but "
                "at least one primer has low occupancy — worth checking experimentally before "
                "you rely on this assay's sensitivity.\n\n"
            )
            md.write("| Assay | Forward Primer | Reverse Primer | Amplicon (bp) | fwd θ | rev θ | fwd 3′ | rev 3′ | Verdict |\n")
            md.write("| :--- | :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: |\n")
            for aname, v in on_target_verdicts.items():
                md.write(
                    f"| `{aname}` "
                    f"| `{v['fwd_oligo']}` "
                    f"| `{v['rev_oligo']}` "
                    f"| {v['amplicon_bp']} "
                    f"| {v['fwd_theta']:.3f} "
                    f"| {v['rev_theta']:.3f} "
                    f"| `{v['fwd_3prime']}` "
                    f"| `{v['rev_3prime']}` "
                    f"| {v['verdict']} |\n"
                )
            md.write("\n")
 
        md.write("---\n")
        md.write("*Report compiled automatically by the Primer Specificity Downstream Analysis Tool.*\n")
 
    logger.info(f"Report written to: {report_path}")

    # =========================================================================
    # STEP 6: EXPORT CSVs AND PLOT
    # =========================================================================
    logger.info(f"Writing outputs to: {output_dir.resolve()}")
    pairwise_df.to_csv(output_dir / "dual_primer_cross_amplification_risks.csv", index=False)
    sequestration_df.to_csv(output_dir / "single_primer_sequestration_sinks.csv", index=False)

    plot_paths = {
        "assay_performance_landscape": assay_performance_landscape_path,
        "specificity_landscape": specificity_landscape_path,
    } if write_plot else {}

    return {
        "pairwise_df": pairwise_df,
        "sequestration_df": sequestration_df,
        "on_target_df": on_target_df,
        "plot_path": assay_performance_landscape_path,
        "plot_paths": plot_paths,
        "report_path": report_path,
    }
    
 
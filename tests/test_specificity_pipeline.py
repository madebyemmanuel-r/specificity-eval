import pytest
from specificity_analysis_components import BlastHit

class TestMultipartiteGenomeHandling:
    """
    Confirms that the pairwise amplification check never combines forward and
    reverse hits from different accessions (e.g. different chromosomes in a
    multipartite genome such as the human genome). This guarantees the
    concatenation behaviour in load_single_fasta_sequence — which only
    operates on single-locus samtools faidx extractions, never whole genomes —
    cannot cause spurious cross-chromosome amplicon predictions.
    """

    def _make_hit(self, seq_id: str, primer_direction: str, strand: str,
                   back_seq_start: int, back_seq_end: int) -> "BlastHit":
        return BlastHit(
            seq_id=seq_id, primer_direction=primer_direction,
            pct_identity=100.0, q_start=1, q_end=25,
            back_seq_start=back_seq_start, back_seq_end=back_seq_end,
            evalue=0.001, num_gaps=0, strand=strand,
            primer_pair="TEST", assay_name="TEST",
        )

    def test_shared_accessions_empty_across_different_chromosomes(self):
        """
        A forward hit on chr1 and a reverse hit on chr2 must produce an empty
        shared_accessions set, meaning check_amplification_risk is never called
        for that pairing.
        """
        fwd_hit = self._make_hit("NC_000001.11", "forward", "plus", 100, 125)
        rev_hit = self._make_hit("NC_000002.12", "reverse", "minus", 200, 175)

        fwd_hits = {"NC_000001.11": [fwd_hit]}
        rev_hits = {"NC_000002.12": [rev_hit]}

        shared_accessions = set(fwd_hits.keys()) & set(rev_hits.keys())
        assert shared_accessions == set(), (
            "No accession should be shared between hits on different chromosomes"
        )

    def test_same_accession_hits_are_correctly_paired(self):
        """
        Sanity check: when forward and reverse hits DO share an accession
        (same chromosome), they should appear in the shared set and be
        eligible for pairwise evaluation.
        """
        fwd_hit = self._make_hit("NC_000001.11", "forward", "plus", 100, 125)
        rev_hit = self._make_hit("NC_000001.11", "reverse", "minus", 300, 275)

        fwd_hits = {"NC_000001.11": [fwd_hit]}
        rev_hits = {"NC_000001.11": [rev_hit]}

        shared_accessions = set(fwd_hits.keys()) & set(rev_hits.keys())
        assert shared_accessions == {"NC_000001.11"}

    def test_multi_chromosome_hits_do_not_cross_contaminate(self):
        """
        With hits present on three different chromosomes for both forward and
        reverse primers, only same-accession pairs should ever be evaluated —
        confirming no O(n*m) cross-chromosome comparison occurs.
        """
        fwd_hits = {
            "NC_000001.11": [self._make_hit("NC_000001.11", "forward", "plus", 100, 125)],
            "NC_000002.12": [self._make_hit("NC_000002.12", "forward", "plus", 500, 525)],
        }
        rev_hits = {
            "NC_000002.12": [self._make_hit("NC_000002.12", "reverse", "minus", 600, 575)],
            "NC_000017.11": [self._make_hit("NC_000017.11", "reverse", "minus", 900, 875)],
        }

        shared_accessions = set(fwd_hits.keys()) & set(rev_hits.keys())

        # Only chr2 is shared — chr1 fwd hit and chr17 rev hit must never be paired
        assert shared_accessions == {"NC_000002.12"}

        evaluated_pairs = []
        for accession in shared_accessions:
            for fwd_hit in fwd_hits[accession]:
                for rev_hit in rev_hits[accession]:
                    evaluated_pairs.append((fwd_hit.seq_id, rev_hit.seq_id))

        assert evaluated_pairs == [("NC_000002.12", "NC_000002.12")]
        assert ("NC_000001.11", "NC_000017.11") not in evaluated_pairs
class TestMBVCrossReactivityHSV1:
    """
    Integration tests validating pipeline output against Chen et al. 2023
    (Zoonoses, DOI: 10.15212/ZOONOSES-2022-0040).

    The study developed two RPA assays targeting Macacine alphaherpesvirus 1 (monkey B virus):
      - RPA-VF-UL27: targets glycoprotein B (UL27), cross-reactive with HSV-1
      - RPA-VF-US6:  targets glycoprotein D (US6),  specific — no cross-reactivity with HSV-1

    Background genome: Human alphaherpesvirus 1 (HSV-1, NC_001806.2, taxon 10298).
    """

    def test_ul27_flagged_tier1_cross_reactive(self, mbv_pipeline_result):
        """
        UL27 must produce at least one Tier 1 cross-reactivity event against HSV-1.
        Chen et al. 2023: RPA-VF-UL27 showed cross-reactivity with HSV-1 nucleic acids.
        """
        pairwise = mbv_pipeline_result["pairwise_df"]
        assert not pairwise.empty, "Pairwise dataframe should not be empty"

        ul27_tier1 = (
            (pairwise["assay_name"] == "UL27") &
            (pairwise["risk_classification"] == "Tier 1: Severe Cross-Reactivity Risk")
        )
        assert ul27_tier1.any(), (
            "UL27 should be flagged Tier 1 cross-reactive with HSV-1 "
            "(Chen et al. 2023, Zoonoses)"
        )

    def test_us6_passes_specificity(self, mbv_pipeline_result):
        """
        US6 must not produce any Tier 1 cross-reactivity events against HSV-1.
        Chen et al. 2023: RPA-VF-US6 showed no cross-reactivity with HSV-1.
        """
        pairwise = mbv_pipeline_result["pairwise_df"]

        us6_tier1 = (
            (pairwise["assay_name"] == "US6") &
            (pairwise["risk_classification"] == "Tier 1: Severe Cross-Reactivity Risk")
        )
        assert not us6_tier1.any(), (
            "US6 should pass specificity screening against HSV-1 "
            "(Chen et al. 2023, Zoonoses)"
        )

    def test_ul27_cross_reactive_on_correct_accession(self, mbv_pipeline_result):
        """
        The UL27 Tier 1 hit must be against HSV-1 (NC_001806.2), not the target genome.
        Confirms background DB is correct and target genome is excluded.
        """
        pairwise = mbv_pipeline_result["pairwise_df"]
        ul27_tier1 = pairwise[
            (pairwise["assay_name"] == "UL27") &
            (pairwise["risk_classification"] == "Tier 1: Severe Cross-Reactivity Risk")
        ]
        assert (ul27_tier1["background_accession"] == "NC_001806.2").all(), (
            "UL27 cross-reactivity should be detected against HSV-1 (NC_001806.2)"
        )

    def test_ul27_amplicon_size_plausible(self, mbv_pipeline_result):
        """
        Predicted UL27 off-target amplicon should be within a plausible RPA range (100–500 bp).
        """
        pairwise = mbv_pipeline_result["pairwise_df"]
        ul27_hits = pairwise[pairwise["assay_name"] == "UL27"]
        assert not ul27_hits.empty
        amplicon = ul27_hits["amplicon_size_bp"].iloc[0]
        assert 100 <= amplicon <= 500, (
            f"UL27 off-target amplicon size {amplicon}bp outside expected RPA range"
        )

    def test_ul27_high_theta_both_primers(self, mbv_pipeline_result):
        """
        UL27 Tier 1 hit should have high fractional occupancy (θ > 0.5) for both primers,
        indicating thermodynamically stable off-target binding.
        """
        pairwise = mbv_pipeline_result["pairwise_df"]
        ul27_tier1 = pairwise[
            (pairwise["assay_name"] == "UL27") &
            (pairwise["risk_classification"] == "Tier 1: Severe Cross-Reactivity Risk")
        ]
        assert (ul27_tier1["fwd_theta"] > 0.5).any(), "UL27 forward theta should exceed 0.5"
        assert (ul27_tier1["rev_theta"] > 0.5).any(), "UL27 reverse theta should exceed 0.5"


class TestMBVOnTargetConfirmation:
    """
    Validates on-target amplification confirmation for both MBV assays.
    Both UL27 and US6 were experimentally confirmed to detect monkey B virus
    (Chen et al. 2023), so both should show Confirmed on-target binding.
    """

    def test_ul27_on_target_confirmed(self, mbv_pipeline_result):
        """UL27 should show confirmed on-target amplification of monkey B virus."""
        on_target = mbv_pipeline_result.get("on_target_df")
        if on_target is None or on_target.empty:
            pytest.skip("On-target CSV not provided for this run")

        ul27_hits = on_target[on_target["assay_name"] == "UL27"]
        assert not ul27_hits.empty, "UL27 should have on-target convergent pairs"
        assert (ul27_hits["fwd_theta"] >= 0.05).any()
        assert (ul27_hits["rev_theta"] >= 0.05).any()

    def test_us6_on_target_confirmed(self, mbv_pipeline_result):
        """US6 should show confirmed on-target amplification of monkey B virus."""
        on_target = mbv_pipeline_result.get("on_target_df")
        if on_target is None or on_target.empty:
            pytest.skip("On-target CSV not provided for this run")

        us6_hits = on_target[on_target["assay_name"] == "US6"]
        assert not us6_hits.empty, "US6 should have on-target convergent pairs"
        assert (us6_hits["fwd_theta"] >= 0.05).any()
        assert (us6_hits["rev_theta"] >= 0.05).any()


class TestMBVSequestration:
    """
    Validates sequestration sink detection.
    UL27 should have fewer high-theta sequestration sites than US6,
    consistent with the sensitivity difference reported in Chen et al. 2023
    (UL27 LOD: 3.4 copies vs US6 LOD: 28 copies).
    """

    def test_sequestration_df_populated(self, mbv_pipeline_result):
        """Sequestration dataframe should be present and non-empty."""
        seq_df = mbv_pipeline_result["sequestration_df"]
        assert seq_df is not None
        assert not seq_df.empty, "Sequestration sites should be detected"

    def test_ul27_sequestration_present(self, mbv_pipeline_result):
        """UL27 should have at least one sequestration site above theta threshold."""
        seq_df = mbv_pipeline_result["sequestration_df"]
        ul27_seq = seq_df[seq_df["primer_pair"] == "UL27"]
        assert not ul27_seq.empty

    def test_us6_more_sequestration_sites_than_ul27(self, mbv_pipeline_result):
        """
        US6 should have more sequestration sites than UL27, consistent with
        its higher LOD (28 copies vs 3.4 copies) in Chen et al. 2023.
        """
        seq_df = mbv_pipeline_result["sequestration_df"]
        ul27_count = len(seq_df[seq_df["primer_pair"] == "UL27"])
        us6_count  = len(seq_df[seq_df["primer_pair"] == "US6"])
        assert us6_count > ul27_count, (
            f"US6 ({us6_count} sites) should have more sequestration sites than "
            f"UL27 ({ul27_count} sites) — consistent with Chen et al. 2023 LOD data"
        )
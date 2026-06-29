import pytest


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
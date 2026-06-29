import pytest
from specificity_check_components import (
    calculate_thermodynamics,
    calculate_equilibrium_constant,
    solve_fraction_bound,
    parse_alignment_window,
    BlastHit,
    AssayConditions,
)


CONDITIONS = AssayConditions()


class TestCalculateThermodynamics:

    def test_returns_dict_for_valid_input(self):
        """Valid primer-template duplex should return a thermodynamics dict."""
        result = calculate_thermodynamics("ATCGATCG", "TAGCTAGC", CONDITIONS)
        assert result is not None
        assert "Enthalpy" in result
        assert "Adjusted Gibbs" in result

    def test_returns_none_for_short_sequence(self):
        """Sequences shorter than 2 bases should return None."""
        assert calculate_thermodynamics("A", "T", CONDITIONS) is None

    def test_gibbs_negative_for_stable_duplex(self):
        """A stable duplex should have negative adjusted Gibbs free energy."""
        result = calculate_thermodynamics("ATCGATCGATCG", "TAGCTAGCTAGC", CONDITIONS)
        assert result is not None
        assert result["Adjusted Gibbs"] < 0


class TestSolveFractionBound:

    def test_perfect_binding_approaches_one(self):
        """Very high K should give theta close to 1.0."""
        theta = solve_fraction_bound(K=1e10, C_A=1e-12, C_B=100e-6)
        assert theta > 0.99

    def test_zero_K_returns_zero(self):
        """K of zero should return theta of 0.0."""
        assert solve_fraction_bound(K=0.0, C_A=1e-12, C_B=100e-6) == 0.0

    def test_theta_bounded_between_zero_and_one(self):
        """Theta must always be in [0, 1]."""
        for K in [1e-5, 1.0, 1e5, 1e10]:
            theta = solve_fraction_bound(K=K, C_A=1e-12, C_B=100e-6)
            assert 0.0 <= theta <= 1.0


class TestCalculateEquilibriumConstant:

    def test_returns_float(self):
        """Should return a float for valid inputs."""
        K = calculate_equilibrium_constant(-40000, -100, 37.0)
        assert isinstance(K, float)

    def test_more_negative_dG_gives_larger_K(self):
        """More stable duplex (more negative dG) should give larger K."""
        K_stable   = calculate_equilibrium_constant(-60000, -150, 37.0)
        K_unstable = calculate_equilibrium_constant(-20000, -60, 37.0)
        assert K_stable > K_unstable


class TestParseAlignmentWindow:

    def _make_hit(self, alignment: str) -> BlastHit:
        return BlastHit(
            seq_id="NC_001806.2", primer_direction="forward",
            pct_identity=100.0, q_start=1, q_end=28,
            back_seq_start=55053, back_seq_end=55077,
            evalue=0.001, num_gaps=0, strand="plus",
            alignment=alignment,
        )

    def test_perfect_3prime_window(self):
        """All-match 3' window should give 0 mismatches and be flagged stable."""
        alignment = (
            "ATCGATCGATCGATCGATCGATCGATCG\n"
            "|||||||||||||||||||||||||||||\n"
            "ATCGATCGATCGATCGATCGATCGATCG"
        )
        hit = self._make_hit(alignment)
        parse_alignment_window(hit, n_bases=5, max_3prime_mismatches=3)
        assert hit.three_prime_mismatches == 0
        assert hit.three_prime_stable is True
        assert hit.three_prime_window == "|||||"

    def test_mismatched_3prime_window(self):
        """4 mismatches in 3' window should be flagged as unstable."""
        alignment = (
            "ATCGATCGATCGATCGATCGATCGATCG\n"
            "||||||||||||||||||||||||-----\n"
            "ATCGATCGATCGATCGATCGATCGATCG"
        )
        hit = self._make_hit(alignment)
        parse_alignment_window(hit, n_bases=5, max_3prime_mismatches=3)
        assert hit.three_prime_mismatches == 4
        assert hit.three_prime_stable is False

    def test_query_and_subject_assigned(self):
        """query_sequence and subject_sequence should be populated after parsing."""
        alignment = (
            "ATCGATCGATCG\n"
            "||||||||||||\n"
            "ATCGATCGATCG"
        )
        hit = self._make_hit(alignment)
        parse_alignment_window(hit, n_bases=5, max_3prime_mismatches=3)
        assert hit.query_sequence == "ATCGATCGATCG"
        assert hit.subject_sequence == "ATCGATCGATCG"

    def test_none_alignment_does_not_crash(self):
        """A hit with no alignment string should return without error."""
        hit = self._make_hit(None)
        parse_alignment_window(hit, n_bases=5, max_3prime_mismatches=3)
        assert hit.three_prime_window is None
import pytest
import pandas as pd
from pathlib import Path
import tempfile, csv

from oligo_handling import parse_and_validate_oligos, _infer_assay_name


# ---------------------------------------------------------------------------
# _infer_assay_name
# ---------------------------------------------------------------------------

class TestInferAssayName:
    """Unit tests for the assay name inference helper."""

    def test_space_delimited_strips_orientation(self):
        """'Eco fliC F1' should infer assay name 'Eco fliC'."""
        assert _infer_assay_name("Eco fliC F1", "Eco fliC F1") == "Eco fliC"

    def test_space_delimited_reverse(self):
        """'Eco fliC R1' should infer assay name 'Eco fliC'."""
        assert _infer_assay_name("Eco fliC R1", "Eco fliC R1") == "Eco fliC"

    def test_underscore_delimited(self):
        """'UL27_F' should infer assay name 'UL27'."""
        assert _infer_assay_name("UL27_F", "UL27") == "UL27"

    def test_underscore_delimited_reverse(self):
        """'UL27_R' should infer assay name 'UL27'."""
        assert _infer_assay_name("UL27_R", "UL27") == "UL27"

    def test_numbered_suffix(self):
        """'rfbE_F1' should infer assay name 'rfbE'."""
        assert _infer_assay_name("rfbE_F1", "rfbE") == "rfbE"

    def test_forward_written_out(self):
        """'SmMIT_Forward' should infer assay name 'SmMIT'."""
        assert _infer_assay_name("SmMIT_Forward", "SmMIT") == "SmMIT"

    def test_single_opaque_token_falls_back_to_pair(self):
        """A single unrecognised token should fall back to primer_pair."""
        assert _infer_assay_name("OPAQUE", "OPAQUE_PAIR") == "OPAQUE_PAIR"


# ---------------------------------------------------------------------------
# parse_and_validate_oligos
# ---------------------------------------------------------------------------

def _write_csv(rows: list[dict], path: Path):
    """Helper — writes a list of dicts as a CSV to path."""
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


class TestParseAndValidateOligos:

    def test_standard_headers_parsed_correctly(self, tmp_path):
        """Standard column names should parse without error."""
        csv_path = tmp_path / "oligos.csv"
        _write_csv([
            {"OligoID": "UL27_F", "Sequence": "ATCGATCGATCG", "PrimerPair": "UL27", "Orientation": "Forward", "AssayName": "UL27"},
            {"OligoID": "UL27_R", "Sequence": "CGATCGATCGAT", "PrimerPair": "UL27", "Orientation": "Reverse", "AssayName": "UL27"},
        ], csv_path)
        df = parse_and_validate_oligos(str(csv_path))
        assert len(df) == 2
        assert set(df["AssayName"]) == {"UL27"}
        assert set(df["PrimerPair"]) == {"UL27"}

    def test_flexible_header_matching(self, tmp_path):
        """Variant column names should be matched flexibly."""
        csv_path = tmp_path / "oligos.csv"
        _write_csv([
            {"name": "UL27_F", "sequence": "ATCGATCGATCG", "assay_id": "UL27", "direction": "Forward"},
            {"name": "UL27_R", "sequence": "CGATCGATCGAT", "assay_id": "UL27", "direction": "Reverse"},
        ], csv_path)
        df = parse_and_validate_oligos(str(csv_path))
        assert len(df) == 2

    def test_assay_name_inferred_when_absent(self, tmp_path):
        """AssayName should be inferred from OligoID when no AssayName column is present."""
        csv_path = tmp_path / "oligos.csv"
        _write_csv([
            {"OligoID": "Eco fliC F1", "Sequence": "ATCGATCGATCG", "PrimerPair": "Eco fliC F1", "Orientation": "Forward"},
            {"OligoID": "Eco fliC R1", "Sequence": "CGATCGATCGAT", "PrimerPair": "Eco fliC R1", "Orientation": "Reverse"},
        ], csv_path)
        df = parse_and_validate_oligos(str(csv_path))
        assert all(df["AssayName"] == "Eco fliC")

    def test_invalid_orientation_raises(self, tmp_path):
        """Invalid orientation value should raise ValueError."""
        csv_path = tmp_path / "oligos.csv"
        _write_csv([
            {"OligoID": "UL27_F", "Sequence": "ATCGATCGATCG", "PrimerPair": "UL27", "Orientation": "Sideways"},
        ], csv_path)
        with pytest.raises(ValueError, match="Invalid orientation"):
            parse_and_validate_oligos(str(csv_path))

    def test_invalid_sequence_characters_raises(self, tmp_path):
        """Non-IUPAC characters in sequence should raise ValueError."""
        csv_path = tmp_path / "oligos.csv"
        _write_csv([
            {"OligoID": "UL27_F", "Sequence": "ATCGXYZ", "PrimerPair": "UL27", "Orientation": "Forward"},
        ], csv_path)
        with pytest.raises(ValueError, match="invalid characters"):
            parse_and_validate_oligos(str(csv_path))

    def test_missing_mandatory_column_raises(self, tmp_path):
        """Missing a mandatory column should raise ValueError."""
        csv_path = tmp_path / "oligos.csv"
        _write_csv([
            {"OligoID": "UL27_F", "Sequence": "ATCGATCGATCG", "Orientation": "Forward"},
        ], csv_path)
        with pytest.raises(ValueError, match="Missing mapping"):
            parse_and_validate_oligos(str(csv_path))

    def test_degenerate_bases_accepted(self, tmp_path):
        """IUPAC degenerate bases should be accepted without error."""
        csv_path = tmp_path / "oligos.csv"
        _write_csv([
            {"OligoID": "UL27_F", "Sequence": "ATCGRYSWKMBDHVN", "PrimerPair": "UL27", "Orientation": "Forward"},
        ], csv_path)
        df = parse_and_validate_oligos(str(csv_path))
        assert len(df) == 1

    def test_returns_dataframe(self, tmp_path):
        """Return type should be a pandas DataFrame."""
        csv_path = tmp_path / "oligos.csv"
        _write_csv([
            {"OligoID": "UL27_F", "Sequence": "ATCGATCGATCG", "PrimerPair": "UL27", "Orientation": "Forward"},
        ], csv_path)
        result = parse_and_validate_oligos(str(csv_path))
        assert isinstance(result, pd.DataFrame)
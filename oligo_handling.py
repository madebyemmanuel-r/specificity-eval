import csv
import logging
from pathlib import Path
import pandas as pd

logger = logging.getLogger(__name__)

def parse_and_validate_oligos(csv_path: str) -> list:
    """
    Parses input CSV containing primer pairs with flexible header matching.
    """
    path = Path(csv_path)
    if not path.exists():
        raise FileNotFoundError(f"Oligo CSV file not found at: {csv_path}")

    # Map standardised internal keys to a list of acceptable user variations
    header_mapping = {
        "OligoID":    ["oligoid", "oligo_id", "oligo id", "primerid", "primer_id", "primer id", "name"],
        "Sequence":   ["sequence", "seq", "nucleotide", "nucleotide sequence", "sequence (5'-3')", "5'-3' sequence"],
        "PrimerPair": ["primerpair", "primer_pair", "primer pair", "pairid", "pair_id", "pair id", "assayid", "assay_id", "assay id"],
        "Orientation":["type", "orientation", "direction", "strand"],
        "AssayName":  ["assayname", "assay_name", "assay name", "setid", "set_id", "panel", "set", "gene", "locus"],
    }

    valid_orientations = {"Forward", "Reverse"}
    parsed_oligos = []

    with open(path, mode="r", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)

        actual_headers = reader.fieldnames if reader.fieldnames else []
        column_map = {}

        for standard_key, variations in header_mapping.items():
            norm_variations = {
                v.strip().lower().replace("_", "").replace("-", "").replace(" ", "")
                for v in variations
            }
            for header in actual_headers:
                norm_header = "".join(header.split()).lower().replace("_", "").replace("-", "")
                if norm_header in norm_variations:
                    column_map[standard_key] = header
                    break

        # AssayName is optional — it can be inferred from OligoID if absent
        mandatory = {"OligoID", "Sequence", "PrimerPair", "Orientation"}
        missing_cols = mandatory - set(column_map.keys())
        if missing_cols:
            raise ValueError(
                f"Could not map CSV columns automatically. Missing mapping for: {missing_cols}. "
                f"Please ensure columns are close to: OligoID, Sequence, PrimerPair, Orientation."
            )

        for row_idx, row in enumerate(reader, start=2):
            oligo_id    = row[column_map["OligoID"]].strip()
            sequence    = row[column_map["Sequence"]].strip().upper()
            primer_pair = row[column_map["PrimerPair"]].strip()
            oligo_type  = row[column_map["Orientation"]].strip().capitalize()

            if not all([oligo_id, sequence, primer_pair, oligo_type]):
                raise ValueError(f"QC Error [Row {row_idx}]: Found an incomplete or empty field.")

            if oligo_type not in valid_orientations:
                raise ValueError(
                    f"QC Error [Row {row_idx}]: Invalid orientation '{oligo_type}'. "
                    f"Must be close to 'Forward' or 'Reverse'."
                )

            invalid_chars = set(sequence) - {"A", "T", "C", "G", "N", "R", "Y", "S", "W", "K", "M", "B", "D", "H", "V"}
            if invalid_chars:
                raise ValueError(
                    f"QC Error [Row {row_idx}]: Sequence contains invalid characters: {invalid_chars}."
                )

            # AssayName: use explicit column if present, otherwise infer from OligoID
            if "AssayName" in column_map and row[column_map["AssayName"]].strip():
                assay_name = row[column_map["AssayName"]].strip()
            else:
                assay_name = _infer_assay_name(oligo_id, primer_pair)

            parsed_oligos.append({
                "OligoID":    oligo_id,
                "Sequence":   sequence,
                "PrimerPair": primer_pair,
                "Orientation": oligo_type,
                "AssayName":  assay_name,
            })

    logging.info(f"QC: Loaded and verified {len(parsed_oligos)} oligos using automated header matching.")
    return pd.DataFrame(parsed_oligos)


def _infer_assay_name(oligo_id: str, primer_pair: str) -> str:
    """
    Strips trailing orientation token (F/R/F1/R1/Fwd/Rev) to recover assay name.
    Falls back to primer_pair if the name is a single unrecognised token.
    e.g. 'Eco fliC F1' -> 'Eco fliC', 'UL27_F' -> 'UL27'
    """
    import re
    normalised = oligo_id.replace("_", " ").strip()
    parts = normalised.split()
    if len(parts) > 1 and re.fullmatch(r'[FRfr](?:wd|ev|orward|everse)?\d*', parts[-1], re.IGNORECASE):
        return " ".join(parts[:-1])
    return primer_pair
import os
import zipfile
import logging
import requests
from pathlib import Path
from NCBI import ncbi_session

logger = logging.getLogger(__name__)

def verify_ncbi_connectivity(timeout_seconds: int = 3) -> bool:
    """
    Quality Control: Pre-flight ping to confirm an active internet connection 
    and accessibility to the NCBI Datasets API.
    """

    logging.info("QC: Verifying internet connectivity to NCBI...")
    # We use the base health/version endpoint for a lightweight, instant check
    ping_url = ("https://api.ncbi.nlm.nih.gov/datasets/v2/version")

    try:

        response = ncbi_session.get(ping_url,timeout=timeout_seconds)

        response.raise_for_status()

        logging.info("QC: NCBI connection verified successfully.")

        return True

    except requests.exceptions.Timeout as exc:

        raise ConnectionError("Connection to the NCBI API timed out.") from exc

    except requests.exceptions.ConnectionError as exc:

        raise ConnectionError("Unable to establish a connection to the NCBI API.") from exc

    except requests.exceptions.HTTPError as exc:

        raise ConnectionError(f"NCBI API returned HTTP status "f"{response.status_code}.") from exc

    except requests.exceptions.RequestException as exc:

        raise ConnectionError(f"Unexpected network error while contacting "f"the NCBI API: {exc}") from exc


def validate_fasta_file(file_path: str) -> bool:
    """
    Quality Control: Checks that the file path exists, points to a valid file,
    and structurally conforms to the expected FASTA format.
    """
    path = Path(file_path)
    
    # QC Step 1: Validate Path and File Existence
    if not path.exists():
        raise FileNotFoundError(f"QC Error: The specified path does not exist: '{file_path}'")
    if not path.is_file():
        raise IsADirectoryError(f"QC Error: Expected a file, but found a directory: '{file_path}'")
    if os.path.getsize(path) == 0:
        raise ValueError(f"QC Error: The file is completely empty (0 bytes): '{file_path}'")
        
    # QC Step 2: Structural Integrity File Check
    with open(path, "r") as f:
        # Read lines until we hit a non-empty line (skips accidental leading whitespace)
        first_char = ""
        for line in f:
            cleaned_line = line.strip()
            if cleaned_line:
                first_char = cleaned_line[0]
                break
                
        if first_char != ">":
            raise ValueError(
                f"QC Error: File '{file_path}' is not a valid FASTA file. "
                f"Expected first non-empty character to be '>', but found '{first_char}'."
            )
            
    logging.info(f"QC: File passed validation successfully: {path.name}")
    return True
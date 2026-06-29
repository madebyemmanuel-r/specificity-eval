import os
import time
import logging
import shutil
import zipfile
import requests
from pathlib import Path
from urllib.parse import quote
from requests import Session
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

logger = logging.getLogger(__name__)

def create_ncbi_session() -> Session:
    """
    Creates a shared Requests session for all NCBI API communication.

    The session automatically retries transient failures (e.g. HTTP 429 or
    temporary server errors) using exponential backoff (takes twice as long each time)
    to improve pipeline robustness, without requiring retry logic at each request.
    """

    # Retry transient HTTP errors up to 5 times using exponential backoff
    # Delays are approximately 0 s, 2 s, 4 s, 8 s and 16 s
    retry = Retry(
        total=5,
        backoff_factor=2,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["GET"],
        raise_on_status=False,
    )

    # Attach the retry policy to all HTTP/HTTPS requests made through this session
    session = Session()
    adapter = HTTPAdapter(max_retries=retry)

    session.mount("https://", adapter)
    session.mount("http://", adapter)

    return session

# Shared session used for all NCBI API requests throughout the pipeline
# Centralising retry behaviour keeps network handling consistent
ncbi_session = create_ncbi_session()

def resolve_identity(organism_query: str, api_key: str | None = None) -> dict:
    """
    Takes a Scientific Name or Taxon ID and resolves it to the best NCBI 
    reference assembly, safely handling URL encoding and malformed API payloads.
    """
    #Rate-limit guard
    time.sleep(0.15 if api_key else 0.5)
    
    # URL-encode the query to handle spaces ("Homo sapiens") and special characters safely
    safe_query = quote(organism_query.strip())
    base_url = "https://api.ncbi.nlm.nih.gov/datasets/v2"
    search_url = f"{base_url}/genome/taxon/{safe_query}/dataset_report"
    
    params = {
        "filters.assembly_source": "refseq",
        "filters.has_annotation": "true"
    }
    if api_key:
        params["api_key"] = api_key
    
    logging.info(f"Resolving identity for: {organism_query}")
    try:
        response = ncbi_session.get(search_url, params=params, timeout=20.0)
        response.raise_for_status() # Catches 4xx and 5xx errors elegantly
    except requests.RequestException as e:
        raise RuntimeError(f"NCBI network request failed for '{organism_query}': {e}")
        
    try:
        data = response.json()
    except ValueError:
        raise RuntimeError(f"NCBI returned malformed, non-JSON data: {response.text[:200]}")
        
    # Quality hierarchy (lower number = higher priority)
    level_priority = {
        "Complete Genome": 0,
        "Chromosome":      1,
        "Scaffold":        2,
        "Contig":          3,
    }
    
    reports = data.get("reports", [])
    if not reports:
        raise ValueError(f"No reference genome found on NCBI for: {organism_query}")
    
    # Pass 1: Sort by release date descending (Newest files move to the front)
    reports.sort(
        key=lambda r: r.get("assembly_info", {}).get("release_date", "1900-01-01"), 
        reverse=True
    )
    
    # Pass 2: Sort by assembly level ascending (Best quality moves to the front)
    # The "newest first" order is preserved within each tier
    reports.sort(
        key=lambda r: level_priority.get(r.get("assembly_info", {}).get("assembly_level", "Contig"), 99)
    )
    
    # Safely return the gold-standard match
    return reports[0]


def fetch_reference_genome(
    best_match: dict, 
    output_dir: str = "data/genomes",
    api_key: str | None = None
) -> dict:
    """
    Queries NCBI Datasets v2 API for an organism, downloads NCBI reference genome assembly package, 
    extracts and caches the FASTA file, and saves the returned metadata.
    """
    #Rate-limit guard
    time.sleep(0.15 if api_key else 0.5)
    
    target_path = Path(output_dir)
    target_path.mkdir(parents=True, exist_ok=True)
    
    # Safe dictionary access
    accession = best_match.get("accession")
    if not accession:
        raise ValueError("Provided metadata dictionary is missing the critical 'accession' key.")
        
    organism_info = best_match.get("organism", {})
    assembly_info = best_match.get("assembly_info", {})

    # Fallback chain handles common naming fields returned by the NCBI Datasets API:
    # Bacterial/eukaryotic assemblies often use "sci_name", while viral assemblies often provide "organism_name" instead.
    # Both conventions are accepted and the org name is in the metadata
    organism_display_name = (
    organism_info.get("sci_name")
    or organism_info.get("organism_name")
    or "unknown"
    )
    
    scientific_name = organism_display_name.replace(" ", "_")
    final_fasta_path = target_path / f"{scientific_name}_{accession}.fna"
    
    genome_meta = {
        "organism":    organism_display_name if organism_display_name != "unknown" else None,
        "common_name": organism_info.get("common_name", "N/A"),
        "tax_id":      organism_info.get("tax_id"),
        "accession":   accession,
        "assembly_name":  assembly_info.get("assembly_name"),
        "assembly_level": assembly_info.get("assembly_level"),
        "seq_length":  int(best_match.get("assembly_stats", {}).get("total_sequence_length", 0)),
        "fasta_path":  str(final_fasta_path)
    }

    # Cache hit logic
    if final_fasta_path.exists() and final_fasta_path.stat().st_size > 0:
        logging.info(f"Cache Hit: {final_fasta_path.name} already exists. Skipping.")
        return genome_meta

    temp_zip_path = target_path / f"{accession}_tmp.zip"
    tmp_fasta_path = target_path / f"{final_fasta_path.name}.tmp"

    try:
        logging.info(f"Downloading data package zip for {accession}...")
        base_url = "https://api.ncbi.nlm.nih.gov/datasets/v2"
        download_url = f"{base_url}/genome/accession/{accession}/download"
        download_params = {"include_annotation_type": "GENOME_FASTA"}
        if api_key:
            download_params["api_key"] = api_key
        
        with ncbi_session.get(download_url, params=download_params, stream=True, timeout=60.0) as r:
            r.raise_for_status()
            with open(temp_zip_path, "wb") as f:
                for chunk in r.iter_content(chunk_size=8192):
                    f.write(chunk)

        # Extraction phase
        logging.info(f"Extracting FASTA file memory-safely...")
        with zipfile.ZipFile(temp_zip_path, "r") as zip_ref:
            for file_info in zip_ref.infolist():
                if file_info.filename.endswith(".fna"):
                    # Use shutil.copyfileobj to stream data without inflating RAM
                    with zip_ref.open(file_info) as source, open(tmp_fasta_path, "wb") as dest:
                        shutil.copyfileobj(source, dest)
                    break
                    
        if not tmp_fasta_path.exists() or tmp_fasta_path.stat().st_size == 0:
            raise FileNotFoundError("No valid .fna file found inside downloaded archive.")

        # Only create the final file if everything succeeded flawlessly
        tmp_fasta_path.rename(final_fasta_path)
        logging.info(f"Successfully saved FASTA to: {final_fasta_path}")

    except Exception as e:
        # Clean up partial work on error so future runs are not impacted
        if tmp_fasta_path.exists():
            tmp_fasta_path.unlink()
        raise RuntimeError(f"Failed to fetch reference genome for {accession}: {e}")
        
    finally:
        # Always clean up the temporary zip package
        temp_zip_path.unlink(missing_ok=True)
            
    return genome_meta
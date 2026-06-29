import argparse
import shutil
import subprocess
import logging
from pathlib import Path

logger = logging.getLogger(__name__)

import hashlib
import re
from pathlib import Path

def generate_db_prefix(fasta_inputs: list[str], base_dir: str = "outputs/blast_dbs") -> tuple[str, str]:
    """
    Generates a human-readable, deterministic prefix and title for the BLAST database 
    using a hybrid combination of species names and an MD5 integrity hash.
    """
    # 1. Sort inputs so order changes don't disrupt the cache signature
    sorted_inputs = sorted([str(Path(p).resolve()) for p in fasta_inputs])
    
    # 2. Compute the unique verification hash
    hash_obj = hashlib.md5()
    for path in sorted_inputs:
        hash_obj.update(path.encode("utf-8"))
    unique_hash = hash_obj.hexdigest()[:8]  # 8 characters is plenty for isolation
    
    # 3. Build a clean, human-readable species label string
    species_labels = []
    for path in sorted_inputs:
        # Extract the filename stem (e.g., "Homo_sapiens_genome" -> "Homo_sapiens_genome")
        stem = Path(path).stem
        # Strip out common technical suffixes to keep filenames concise
        clean_name = re.sub(r'(_genome|_assembly|_ref|data|v\d+)', '', stem, flags=re.IGNORECASE)
        # Convert to PascalCase/Alphanumeric only (e.g., "homo_sapiens" -> "HomoSapiens")
        clean_name = "".join([w.capitalize() for w in re.findall(r'[A-Za-z0-9]+', clean_name)])
        
        if clean_name:
            species_labels.append(clean_name[:12])  # Cap individual species names at 12 characters
            
    # Fallback if names are completely non-alphanumeric
    if not species_labels:
        species_labels = ["CompiledBackground"]
        
    # Combine up to the first 3 species names to keep the filename from growing too wide
    combined_labels = "_".join(species_labels[:3])
    if len(species_labels) > 3:
        combined_labels += "_and_others"
        
    # 4. Construct final system parameters
    db_title = f"RPA_Background_{combined_labels}_{unique_hash}"
    output_prefix = str(Path(base_dir) / f"bg_db_{combined_labels}_{unique_hash}")
    
    return output_prefix, db_title

def build_blast_db(fasta_inputs: list[str], output_prefix: str, db_title: str) -> str:
    """
    Builds a unified BLAST database and corresponding samtools index.

    Evaluates the existence and integrity of previously compiled database binaries 
    to leverage cached assets.
    
    If the files are not already there, it unifies multiplexed input files 
    (if there is more than one input fasta), executes makeblastdb and generates a 
    samtools index to enable downstream flanking-window context lookups.
    """
    # Ensure external system dependencies exist in the active environment PATH
    for cmd in ["makeblastdb", "samtools"]:
        if not shutil.which(cmd):
            raise EnvironmentError(f"Missing systemic dependency: {cmd} must be in PATH.")

    db_path = Path(output_prefix)
    db_dir = db_path.parent
    db_dir.mkdir(parents=True, exist_ok=True)

    # Managed workspace path for the local FASTA asset
    combined_fasta = db_path.with_suffix(".fna")

    # 0. If the files exist and aren't empty, a cache hit is recorded and the db build is skipped
    required_extensions = [".nsq", ".nin", ".nhr"]
    cache_valid = True
    
    # Verify BLAST binaries exist and are populated
    for ext in required_extensions:
        check_path = Path(f"{db_path}{ext}")
        if not check_path.exists() or check_path.stat().st_size == 0:
            cache_valid = False
            break
            
    # Verify accompanying samtools index file exists and is populated
    fai_path = Path(f"{combined_fasta}.fai")
    if not fai_path.exists() or fai_path.stat().st_size == 0:
        cache_valid = False

    if cache_valid:
        logging.info(f"Cache Hit: A valid BLAST database already exists at '{output_prefix}'. Skipping build.")
        return str(db_path)

    # 1. Merge all species FASTAs into one unified background file - if necessary
    # Otherwise, use the single input file directly to build the database
    if len(fasta_inputs) > 1:
        logging.info(f"Merging {len(fasta_inputs)} genome files into unified asset: {combined_fasta}")
        with open(combined_fasta, "w") as outfile:
            for fasta_file in fasta_inputs:
                with open(fasta_file, "r") as infile:
                    shutil.copyfileobj(infile, outfile)
    else:
        logging.info(f"Localising single genome file directly to workspace: {combined_fasta}")
        shutil.copy2(fasta_inputs[0], combined_fasta)

    # 2. Build the BLAST database (.ndb, .nsq, etc.)
    logging.info(f"Running makeblastdb for database title: {db_title}")
    blast_cmd = [
        "makeblastdb",
        "-in", str(combined_fasta),
        "-dbtype", "nucl",
        "-parse_seqids",
        "-out", str(db_path),
        "-title", db_title
    ]
    
    # Save the program's output text so we can read the exact error message if it crashes
    run_blast = subprocess.run(blast_cmd, capture_output=True, text=True)
    if run_blast.returncode != 0:
        raise RuntimeError(f"makeblastdb compilation failed:\n{run_blast.stderr}")

    # 3. Build the samtools index (.fai) for downstream tool use
    logging.info("Running samtools faidx context indexing...")
    run_samtools = subprocess.run(["samtools", "faidx", str(combined_fasta)], capture_output=True, text=True)
    if run_samtools.returncode != 0:
        raise RuntimeError(f"samtools faidx indexing failed:\n{run_samtools.stderr}")
    
    logging.info("Database creation complete and cached for future operational runs.")
    return str(db_path)
import os
import logging
import subprocess
import tempfile
import pandas as pd
import argparse
import re

logger = logging.getLogger(__name__)

def load_single_fasta_sequence(fasta_path: str) -> str:
    """
    Loads sequence content from a FASTA file.
    
    If multiple FASTA records are present, their sequences are concatenated in
    file order and a warning is logged. 
    
    This behaviour is compatible withfragmented or multipartite genome assemblies,
    but should be considered when analysing genomes comprising multiple replicons.
    """
    sequence = []
    header_count = 0

    with open(fasta_path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                header_count += 1
                if header_count > 1:
                    logging.warning(
                        f"Multiple FASTA records detected in '{fasta_path}' "
                        f"(record {header_count}: {line[:60]!r}). "
                        f"This is expected for fragmented or multipartite genome assemblies; "
                        f"records will be concatenated in file order."
                    )
                continue
                
            sequence.append(line)

    return "".join(sequence).upper()

def get_complement(seq: str, reverse: bool = False, rule: str = "N2N") -> str:
    """
    Obtains complement or reverse complement.
    """
    base_map = {
        "A": "T",
        "T": "A",
        "C": "G",
        "G": "C",
        "U": "A",
        "-": "N"
    }

    seq_comp = []

    for base in seq.upper():
        if base == "N":
            seq_comp.append("-" if rule == "N2-" else "N")
        else:
            seq_comp.append(base_map.get(base, "N"))

    result = "".join(seq_comp)

    return result[::-1] if reverse else result


def secondary_structure_check(seq_one: str, seq_two: str, reverse_orientation: bool, fixed_back: bool = False) -> tuple:
    nuc_dict = {"A": "T", "T": "A", "C": "G", "G": "C"}

    max_binding_sites = 0
    max_binding_string = ""
    max_binding_percentage = 0

    hard_fail = False
    hard_fail_string = ""

    max_length = len(seq_one) + len(seq_two) - 1

    seq_one_norm = seq_one + " " * (max_length - len(seq_one))
    seq_two_norm = seq_two + " " * (max_length - len(seq_two))
    seq_two_reversed = seq_two[::-1] + " " * (max_length - len(seq_two))

    combos = (
        [
            (seq_two_norm, seq_one_norm),
            (seq_one_norm, seq_two_norm),
            (seq_one_norm, seq_two_reversed),
            (seq_two_reversed, seq_one_norm),
        ]
        if reverse_orientation
        else [
            (seq_two_norm, seq_one_norm),
            (seq_one_norm, seq_two_norm),
        ]
    )

    for var_pair in combos:

        for i in range(var_pair[0].count(" ") + 1):

            dynamic_seq = (
                var_pair[0]
                if i == 0
                else " " * i + var_pair[0][:-i]
            )

            dynamic_list = list(dynamic_seq)
            fixed_list = list(var_pair[1])

            possible_binding_string = ""

            for pos in range(len(dynamic_list)):

                if (
                    dynamic_list[pos] in nuc_dict
                    and nuc_dict[dynamic_list[pos]] == fixed_list[pos]
                ):
                    possible_binding_string += "|"
                else:
                    possible_binding_string += "-"

            number_of_matches = possible_binding_string.count("|")

            complete_binding_string = (
                f"{dynamic_seq}\n"
                f"{possible_binding_string}\n"
                f"{var_pair[1]}"
            )

            if var_pair[1] == seq_one_norm:
                five_index = 0
                three_index = len(seq_one)
            else:
                five_index = i
                three_index = i + len(seq_one)

            upper_limit = five_index + 22
            lower_limit = max(0, three_index - 22)

            five_prime_counts = possible_binding_string[five_index:upper_limit]
            three_prime_counts = possible_binding_string[lower_limit:three_index]

            weightings = [3, 2, 1.5]

            adjusted_weights = []

            for orig_score, indexes in [
                (five_prime_counts, [0, 1, 2]),
                (three_prime_counts, [-1, -2, -3]),
            ]:

                temp_score = (
                    orig_score.count("|")
                    - orig_score.count("-")
                )

                for loc, weight in zip(indexes, weightings):

                    if -len(orig_score) <= loc < len(orig_score):

                        if orig_score[loc] == "|":
                            temp_score += weight - 1

                adjusted_weights.append(temp_score)

            if max(adjusted_weights) >= 21.5:
                hard_fail = True
                hard_fail_string = complete_binding_string

            if max_binding_sites < number_of_matches:

                max_binding_sites = number_of_matches
                max_binding_string = complete_binding_string

                denom = (
                    len(seq_one)
                    if fixed_back
                    else min(len(seq_one), len(seq_two))
                )

                max_binding_percentage = (
                    max_binding_sites / denom
                ) * 100

    return (
        max_binding_percentage,
        max_binding_string,
        hard_fail,
        hard_fail_string,
    )


def extract_background_sequence(db_path: str, seq_id: str, start: int, end: int, reverse: bool = False) -> str:
    """
    Extracts genomic sequence using samtools faidx against the compiled background reference database.
    """
    with tempfile.NamedTemporaryFile(suffix=".fa", delete=False) as tmp:
        tmp_path = tmp.name

    try:
        # 1. Use the explicit compiled asset path string directly 
        fasta_source = f"{db_path}.fna" if not db_path.endswith(".fna") else db_path
        
        # 2. Build the command string cleanly using single quotes around coordinates
        reverse_flag = "-i " if reverse else ""
        cmd = f"samtools faidx {reverse_flag}{fasta_source} '{seq_id}:{start}-{end}' -o {tmp_path}"

        logger.debug(f"Executing background extraction: {cmd}")

        run = subprocess.run(
            cmd,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )

        if run.returncode != 0:
            raise RuntimeError(
                f"samtools faidx failed:\n{run.stderr}\nCommand executed: {cmd}"
            )

        return load_single_fasta_sequence(tmp_path)

    finally:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)


def Blastn_Background_Check(oligo: dict, args: argparse.Namespace, output_dir: str, db_path: str, hit_prefix: str = "bg") -> pd.DataFrame:
    """
    Executes BLAST + secondary structure analysis for a single primer.
    """
    oligo_id = oligo["OligoID"]
    query_sequence = oligo["Sequence"]

    fasta_path = os.path.join(output_dir, f"{oligo_id}.fa")

    with open(fasta_path, "w") as handle:
        handle.write(f">{oligo_id}\n{query_sequence}\n")

    output_filename = os.path.join(output_dir, f"{oligo_id}_Blastn_Output.csv")

    # Establish sensitivity profiles
    penalty = "-2"
    word_size = "4"

    if args.sensitivity == "Advanced":
        penalty = "-1"
    elif args.sensitivity == "Fast":
        penalty = "-3"
        word_size = "7"
    elif args.sensitivity == "Very Advanced":
        penalty = "-1"
        word_size = "3"

    blast_cmd = (
        f"blastn "
        f"-word_size {word_size} "
        f"-gapopen 5 -gapextend 2 -reward 1 "
        f"-penalty {penalty} "
        f"-evalue {args.evalue} "
        f"-query {fasta_path} "
        f"-db {db_path} "  # Standardised to use the dynamic db_path
        f"-out {output_filename} "
        f'-outfmt "10 sseqid pident qstart qend sstart send evalue gaps sstrand"'
    )

    logger.info(f"Running BLAST for {oligo_id}")

    blast_run = subprocess.run(
        blast_cmd,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )

    if blast_run.returncode != 0:
        raise RuntimeError(f"BLAST failed for {oligo_id}:\n{blast_run.stderr}")

    columns = [
        "Background_SeqID", "Percentage Identity", "qStart", "qEnd",
        "BackSeq_Start", "BackSeq_End", "Evalue", "Number Of Gaps", "Strand"
    ]

    try:
        blast_df = pd.read_csv(output_filename, header=None)
    except pd.errors.EmptyDataError:
        fallback_cols = ["OligoID"] + columns
        blast_df = pd.DataFrame(columns=fallback_cols)

    if blast_df.empty:
        blast_df.columns = columns
        return blast_df

    blast_df.columns = columns
    blast_df.insert(0, "OligoID", oligo_id)
    blast_df["PrimerPair"] = oligo["PrimerPair"]
    blast_df["AssayName"]  = oligo["AssayName"]
    blast_df["Advanced Cross Reactivity Percentage String"] = ""
    blast_df["Cross Reactivity Hard Fail"] = False

    for idx in blast_df.index:
        row = blast_df.loc[idx]
        seq_id = row["Background_SeqID"]
        q_start = int(row["qStart"])
        q_end = int(row["qEnd"])
        back_start = int(row["BackSeq_Start"])
        back_end = int(row["BackSeq_End"])
        strand = row["Strand"]
            
        # Remove BLAST pipeline artifacts from the SeqID string
        # E.g. converts 'ref|NC_076164.1|' or 'gi|...|ref|NC_076164.1|' into just 'NC_076164.1'
        if "ref|" in seq_id:
            seq_id = seq_id.split("ref|")[-1].strip("|")
        elif "|" in seq_id:
            seq_id = seq_id.split("|")[-2]

        blast_df.loc[idx, "Background_SeqID"] = seq_id
            
        shorter_background = False

        if strand == "plus":
            upper_extension = (len(query_sequence) - q_end) + back_end
            lower_extension = back_start - q_start + 1

            if lower_extension < 0:
                lower_extension = 0
                shorter_background = True

            # Updated parameter routing
            extracted_seq = extract_background_sequence(
                db_path, seq_id, lower_extension, upper_extension, reverse=False
            )
        else:
            upper_extension = back_start + q_start - 1
            lower_extension = back_end - (len(query_sequence) - q_end)

            if lower_extension < 0:
                lower_extension = 0
                shorter_background = True

            # Updated parameter routing
            extracted_seq = extract_background_sequence(
                db_path, seq_id, lower_extension, upper_extension, reverse=True
            )

        complement_seq = get_complement(extracted_seq)

        if len(complement_seq) != len(query_sequence):
            blast_df.loc[idx, "Advanced Cross Reactivity Percentage String"] = "LENGTH_MISMATCH_SKIPPED"
            continue

        (_, binding_string, hard_fail, _) = secondary_structure_check(
            query_sequence, complement_seq, False, shorter_background
        )

        # Cleaning Step: Strip leading digits and tabs/spaces from the top line
        if binding_string and isinstance(binding_string, str):
                    lines = binding_string.splitlines()
                    if len(lines) >= 3:
                        import re
                        
                        # 1. Strip any coordinate digits, tabs and whitespaces from line 1
                        lines[0] = re.sub(r'^\d+\s*', '', lines[0]).rstrip()
                        lines[1] = re.sub(r'^\d+\s*', '', lines[1]).rstrip()
                        lines[2] = re.sub(r'^\d+\s*', '', lines[2]).rstrip()
                        
                        # 2. Trim trailing whitespaces from all lines safely
                        lines[0] = lines[0].rstrip()
                        lines[1] = lines[1].rstrip()
                        lines[2] = lines[2].rstrip()
                        
                        # 3. Truncate all lines to match the exact length of the primer query
                        target_len = len(query_sequence)
                        lines[0] = lines[0][:target_len]
                        lines[1] = lines[1][:target_len]
                        lines[2] = lines[2][:target_len]
                        
                        binding_string = "\n".join(lines)
            
        blast_df.loc[idx, "Advanced Cross Reactivity Percentage String"] = binding_string
        blast_df.loc[idx, "Cross Reactivity Hard Fail"] = hard_fail

    blast_df.to_csv(output_filename, index=False)
    logger.info(f"Completed analysis for {oligo_id}")
    return blast_df


def run_blast_searches(oligos_df: pd.DataFrame, args: argparse.Namespace, analysis_dir: str, db_path: str, hit_prefix: str = "bg") -> tuple[str, str]:
    """
    Orchestrates the individual background alignment and structural verification 
    pipelines for a collection of oligonucleotide primers.
    """
    fwd_dfs = []
    rev_dfs = []

    if isinstance(oligos_df, list):
        oligos_df = pd.DataFrame(oligos_df)
        
    for _, row in oligos_df.iterrows():
        oligo_dict = row.to_dict()
        orientation = str(oligo_dict.get("Orientation", "")).lower()
        
        # Clean execution step using the verified parameters
        result_df = Blastn_Background_Check(oligo_dict, args, analysis_dir, db_path, hit_prefix)
        
        if "fwd" in orientation or "forward" in orientation:
            fwd_dfs.append(result_df)
        else:
            rev_dfs.append(result_df)
            
    columns = [
        "OligoID", "Background_SeqID", "Percentage Identity", "qStart", "qEnd", "BackSeq_Start", 
        "BackSeq_End", "Evalue", "Number Of Gaps", "Strand", 
        "Advanced Cross Reactivity Percentage String", "Cross Reactivity Hard Fail"
    ]

    master_fwd = pd.concat(fwd_dfs, ignore_index=True) if fwd_dfs else pd.DataFrame(columns=columns)
    master_rev = pd.concat(rev_dfs, ignore_index=True) if rev_dfs else pd.DataFrame(columns=columns)
    
    # Use hit_prefix to dynamically differentiate filenames
    fwd_path = os.path.join(analysis_dir, f"{hit_prefix}_fwd_hits.csv")
    rev_path = os.path.join(analysis_dir, f"{hit_prefix}_rev_hits.csv")
    
    master_fwd.to_csv(fwd_path, index=False)
    master_rev.to_csv(rev_path, index=False)
    
    return fwd_path, rev_path
import logging, subprocess, os, sys, re, math
from utility import file_format, run_subprocess, make_directory, file_exists_and_not_empty
from Bio import SeqIO
import pandas as pd

#Authors: "Amrei Knuth", "Lars Gabriel"
#Credits: "Katharina Hoff"
#Email:"lars.gabriel@uni-greifswald.de"
#Date: "Janurary 2025"

logger = logging.getLogger(__name__)

def indexing(genome_fasta, hisat_path, threads, index_name="genome"):
    """
    Build a genome index using HISAT2.

    This function runs the `hisat2-build` command to create a genome index for use with HISAT2. 
    It ensures the indexing process is logged and handles any errors gracefully.

    :param genome_fasta: Path to the genome FASTA file to be indexed.
    :param hisat_path: Path to the HISAT2 executable directory.
    :param threads: Number of threads to use for indexing.
    :param index_name: Name for the generated index files (default: "genome").
    """
    hisat2_build = os.path.join(hisat_path, "hisat2-build")
    # Construct the command
    command = [
        hisat2_build,
        "--quiet",
        "-p",
        str(threads),
        genome_fasta,
        index_name
    ]

    # Run and log indexing
    logger.info(f"Building genome index from {genome_fasta} using {threads} threads...")
    error_msg  = f"Indexing failed. Check logs for details."
    run_subprocess(command, error_message=error_msg)
    logger.info("Indexing completed successfully.")


def mapping_short(rnaseq_paired_sets, rnaseq_single_sets, hisat_path, threads, 
                    output_dir='./', use_existing=False):
    """
    Map RNA-Seq short reads (single-end and paired-end) to a genome using HISAT2.

    This function maps single-end and paired-end RNA-Seq short reads to a reference genome 
    index using the HISAT2 aligner. It handles input file format validation, error handling, 
    and returns the list of generated alignment files.

    :param rnaseq_paired_sets: List of paired-end RNA-Seq read files. Files should be ordered as [read1_1, read1_2, read2_1, read2_2, ...].
    :param rnaseq_single_sets: List of single-end RNA-Seq read files.
    :param hisat_path: Path to the directory containing the HISAT2 executable.
    :param threads: Number of threads to use for HISAT2 (default: 4).
    :return: List of alignment file paths generated by HISAT2.
    """
    alignments_list = []
    hisat2_executable = f"{hisat_path}/hisat2"

    make_directory(output_dir)

    def determine_format_options(file_path):
        """Determine the file format (fasta, fastq, or unknown) of a given file."""
        format = file_format(file_path)
        if format == "unknown":
            error_msg = f"Error: Unknown file format for {file_path}. Please provide a FASTA or FASTQ file."
            logging.error(error_msg)
            raise ValueError(error_msg)
        return "-f" if format == "fasta" else ""

    # Process single-end reads
    if rnaseq_single_sets:
        format_option = determine_format_options(rnaseq_single_sets[0])
        string_with_sets = ",".join(rnaseq_single_sets)
        output_single = f"{output_dir}/alignment_single_rnaseq.sam"        
        
        hisat2_single_command = [
            hisat2_executable,
            format_option,
            "-x", "genome",
            "-U", string_with_sets,
            "--dta",
            "-p", str(threads),
            "-S", output_single
        ]

        logging.info("Mapping single-end short reads to the genome...")
        if use_existing and file_exists_and_not_empty(output_single):
            logging.info(f"Using existing file instead of rerunning program: {output_single}")
        else:
            error_msg  = f"Error during mapping of single-end short reads."
            run_subprocess(hisat2_single_command, error_message=error_msg)
            logger.info("Mapping of single-end short reads completed successfully.")
        alignments_list.append(output_single)

    # Process paired-end reads
    if rnaseq_paired_sets:
        format_option = determine_format_options(rnaseq_paired_sets[0])
        # !THIS HAS TO BE IMPROVED
        string_with_first = ",".join(rnaseq_paired_sets[0::2])  # Read 1 files
        string_with_second = ",".join(rnaseq_paired_sets[1::2])  # Read 2 files
        output_paired = f"{output_dir}/alignment_paired_rnaseq.sam"
        
        hisat2_paired_command = [
            hisat2_executable,
            format_option,
            "-x", "genome",
            "-1", string_with_first,
            "-2", string_with_second,
            "--dta",
            "-p", str(threads),
            "-S", output_paired
        ]

        logging.info("Mapping paired-end short reads to the genome...")
        if use_existing and file_exists_and_not_empty(output_paired):
            logging.info(f"Using existing file instead of rerunning program: {output_paired}")
        else:
            error_msg  = f"Error during mapping of paired-end short reads."
            run_subprocess(hisat2_paired_command, error_message=error_msg)
            logging.info("Mapping of paired-end short reads completed successfully.")
        
        alignments_list.append(output_paired)
        
    return alignments_list

def mapping_long(genome, isoseq_sets, minimap_path, threads=4, output_dir="./", use_existing=False):
    """Map long reads (IsoSeq data) to a reference genome using Minimap2.

    This function uses Minimap2 to map IsoSeq long reads to a reference genome. The output
    is a SAM file containing the alignments.

    :param genome: Path to the reference genome FASTA file.
    :param isoseq_sets: List of paths to IsoSeq read files.
    :param minimap_path: Path to the directory containing the Minimap2 executable.
    :param threads: Number of threads to use for Minimap2 (default: 4).
    :param output_dir: Path to the output directory (default: "./").
    :return: Path to the generated SAM file.
    """

    # Construct the path to the Minimap2 executable
    minimap_executable = os.path.join(minimap_path, "minimap2")

    make_directory(output_dir)
    output_sam = f"{output_dir}/alignment_isoseq.sam"
    # Construct the Minimap2 command
    command = [
        minimap_executable,
        "-ax", "splice:hq",  # Alignment preset for high-quality long reads with spliced alignment
        "-uf",  # Unspliced mode for full-length cDNA alignment
        genome,
        "-t", str(threads)
    ] + isoseq_sets + ["-o", output_sam]

    logging.info("Mapping long reads to the genome using Minimap2...")
    if use_existing and file_exists_and_not_empty(output_sam):
        logging.info(f"Using existing file instead of rerunning program: {output_sam}")
    else:
        error_msg  = f"Minimap2 mapping failed. Check logs for details."
        run_subprocess(command, error_message=error_msg)
        logging.info(f"Mapping of long reads completed successfully. Output written to {output_sam}")
    
    return output_sam


def assembling(alignment_rnaseq, alignment_isoseq, stringtie_path, 
            threads=4, output_dir="./", mode="mixed", use_existing=False):
    """
    Assemble reads into transcripts using StringTie.

    This function assembles RNA-Seq and/or IsoSeq alignments into transcripts using StringTie,
    based on the specified mode. The output is a GTF file containing the assembled transcripts.

    :param alignment_rnaseq: Path to the RNA-Seq alignment BAM file.
    :param alignment_isoseq: Path to the IsoSeq alignment BAM file.
    :param stringtie_path: Path to the directory containing the StringTie executable.
    :param threads: Number of threads to use for StringTie (default: 4).
    :param output_dir: Path to the output directory (default: "./").
    :param mode: Assembly mode ("mixed", "rnaseq", or "isoseq", default: "mixed").
    :return: Path to the generated GTF file.
    """
    # Construct the path to the StringTie executable
    stringtie_executable = os.path.join(stringtie_path, "stringtie")

    # Ensure the output directory exists
    make_directory(output_dir)

    # Define the output GTF file
    output_gtf = os.path.join(output_dir, "transcripts.gtf")

    # Define base command
    command = [
        stringtie_executable,
        "-p", str(threads),
        "-o", output_gtf
    ]

    # Add mode-specific options
    mode_options = {
        "mixed": ["--mix", alignment_rnaseq, alignment_isoseq],
        "rnaseq": [alignment_rnaseq],
        "isoseq": ["-L", alignment_isoseq]
    }

    if mode not in mode_options:
        raise ValueError(f"Invalid mode: {mode}. Choose from 'mixed', 'rnaseq', or 'isoseq'.")

    # Extend the base command with mode-specific options
    command += mode_options[mode]
    log_message = f"Assembling short and long reads ({mode} mode)..."

    # Log the process
    logging.info(log_message)
    if use_existing and file_exists_and_not_empty(output_gtf):
        logging.info(f"Using existing file instead of rerunning program: {output_gtf}")
    else:
        error_msg = f"StringTie assembly failed in {mode} mode. Check logs for details."
        # Run the subprocess
        run_subprocess(command, error_message=error_msg)
        logging.info(f"Assembly completed successfully. Output written to {output_gtf}")

    return output_gtf

def orfsearching(genome_fa, transcripts_gtf, transdecoder_path, transdecoder_util_path,
                        output_dir="./", use_existing=False):
    """
    Search for open reading frames (ORFs) using TransDecoder.

    This function uses TransDecoder modules to process a genome and transcripts GTF file, 
    generate a FASTA file, extract long ORFs, and predict likely coding regions.

    :param genome_fa: Path to the genome FASTA file.
    :param transcripts_gtf: Path to the transcripts GTF file.
    :param transdecoder_path: Path to the directory containing TransDecoder executables.
    :param transdecoder_util_path: Path to the directory containing TransDecoder utility scripts.
    :param output_dir: Path to the output directory (default: "./").
    :param use_existing: Whether to use existing files instead of rerunning (default: False).
    :return: Path to the predicted coding regions directory.
    """
    # Ensure the output directory exists
    make_directory(output_dir)

    # Step 1: Generate the transcripts FASTA file
    transcripts_fasta = os.path.join(output_dir, "transcripts.fasta")
    gtf_to_fasta_script = os.path.join(transdecoder_util_path, "gtf_genome_to_cdna_fasta.pl")
    fasta_file_command = [
        gtf_to_fasta_script,
        transcripts_gtf,
        genome_fa
    ]

    logging.info("Creating a FASTA file using the genome and transcripts.gtf files...")
    if use_existing and file_exists_and_not_empty(transcripts_fasta):
        logging.info(f"Using existing FASTA file: {transcripts_fasta}")
    else:
        error_msg = "Failed to generate transcripts FASTA file. Check logs for details."
        with open(transcripts_fasta, "w") as output_file:
            run_subprocess(fasta_file_command, stdout=output_file, error_message=error_msg,
                                capture_output=False)
        logging.info(f"Created transcripts FASTA file successfully: {transcripts_fasta}")

    # Step 2: Extract long ORFs
    transdecoder_longorfs = os.path.join(transdecoder_path, "TransDecoder.LongOrfs")
    long_orf_command = [
        transdecoder_longorfs,
        "-O",
        output_dir,
        "-t",
        transcripts_fasta
    ]

    logging.info("Extracting long open reading frames (ORFs)...")
    transcripts_pep = os.path.join(output_dir, "transcripts.fasta.transdecoder.pep")
    if use_existing and file_exists_and_not_empty(transcripts_pep):
        logging.info(f"Using existing ORFs extraction results in: {output_dir}")
    else:
        error_msg = "Failed to extract long ORFs. Check logs for details."
        run_subprocess(long_orf_command, error_message=error_msg)
        logging.info(f"Extracted long ORFs successfully in: {output_dir}")

    # Step 3: Predict likely coding regions
    transdecoder_predict = os.path.join(transdecoder_path, "TransDecoder.Predict")
    predict_command = [
        transdecoder_predict,
        "-O",
        output_dir,
        "-t",
        transcripts_fasta
    ]

    logging.info("Predicting likely coding regions...")
    transcripts_gff = os.path.join(output_dir, "transcripts.fasta.transdecoder.gff3")
    if use_existing and file_exists_and_not_empty(transcripts_gff):
        logging.info(f"Using existing prediction results in: {output_dir}")
    else:
        error_msg = "Failed to predict likely coding regions. Check logs for details."
        run_subprocess(predict_command, error_message=error_msg)
        logging.info(f"Predicted likely coding regions successfully in: {output_dir}")

    return transcripts_pep, transcripts_fasta

def shorten_incomplete_orfs(transdecoder_pep, output_path="shortened_candidates.pep"):
    """
    Shorten incomplete ORFs in a TransDecoder peptide file by trimming sequences 
    to start at the first methionine (M) codon, and update their descriptions.

    ORFs with the type `5prime_partial` or `internal` are considered incomplete.
    The function calculates new coordinates and lengths, and outputs shortened ORFs.

    :param transdecoder_pep: Path to the TransDecoder peptide FASTA file.
    :param output_dir: Directory to store the output file (default: "./").
    :param output_file: Name of the output file (default: "shortened_candidates.pep").
    :return: Path to the output file containing shortened ORFs.
    """
    # Open the output file for writing
    with open(output_path, "w") as output:
        # Parse the input peptide file
        for record in SeqIO.parse(transdecoder_pep, "fasta"):
            # Identify incomplete ORFs (lacking a start codon)
            if "type:5prime_partial" in record.description or "type:internal" in record.description:
                # Find the position of the first methionine (M)
                m_position = record.seq.find("M")

                # If an M is found, shorten the sequence
                if m_position != -1:
                    # Shorten the ORF to start at the first M
                    record.seq = record.seq[m_position:]

                    # Extract coordinates from the description
                    description_parts = record.description.split(" ")
                    coords_match = re.search(r":(\d+)-(\d+)\([\+\-]\)", description_parts[7])

                    if coords_match:
                        old_start = int(coords_match.group(1))
                        stop = int(coords_match.group(2))

                        # Calculate the new start position
                        new_start = old_start + m_position * 3
                        new_length = len(record.seq) - 1  # Subtract 1 to exclude the stop codon (*)

                        # Update the description with the new start and length
                        record.description = re.sub(f"{old_start}-{stop}", f"{new_start}-{stop}", record.description)
                        record.description = re.sub(r"len:\d+", f"len:{new_length}", record.description)

                        # Write the updated record to the output file
                        SeqIO.write(record, output, "fasta")

    logging.info(f"Shortened ORFs successfully written to {output_path}")
    return output_path


def make_diamond_db(protein_file, diamond_path, output_dir="./", db_name="protein_db"):
    """
    Create a DIAMOND protein database from a given protein file.

    This function uses DIAMOND to create a protein database from the specified input file.
    If the database already exists and `use_existing` is True, it skips re-creating the database.

    :param protein_file: Path to the input protein FASTA file.
    :param diamond_path: Path to the directory containing the DIAMOND executable.
    :param output_dir: Directory to store the database (default: "./").
    :param db_name: Name of the database to create (default: "protein_db").
    :return: Path to the created or existing database.
    """
    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Define the path to the database
    db_path = os.path.join(output_dir, db_name)

    # Construct the path to the DIAMOND executable
    diamond_executable = os.path.join(diamond_path, "diamond")

    # Construct the DIAMOND makedb command
    command = [
        diamond_executable,
        "makedb",
        "--in", protein_file,
        "-d", db_path
    ]

    # Log the process
    logging.info("Creating a DIAMOND protein database...")
    error_msg = "Error creating DIAMOND database. Check logs for details."
    run_subprocess(command, error_message=error_msg)
    logging.info(f"Protein database created successfully: {db_path}.dmnd")

    return db_path

def validating_orfs(transdecoder_pep, protein_db_path, diamond_path, 
            output_path="output.tsv", use_existing=False):
    """
    Validate ORFs by searching for matches in a protein database using DIAMOND blastp.

    This function uses DIAMOND blastp to compare predicted ORFs against a protein database and outputs
    the results in a TSV file.

    :param transdecoder_pep: Path to the TransDecoder peptide FASTA file.
    :param protein_db_path: Path to the DIAMOND protein database.
    :param diamond_path: Path to the directory containing the DIAMOND executable.
    :param output_path: Path to the output TSV file (default: "output.tsv").
    :param use_existing: Whether to use an existing TSV file if it already exists (default: False).
    :return: Path to the output TSV file.
    """

    # Check if the output TSV file already exists
    if use_existing and os.path.isfile(output_path):
        logging.info(f"Using existing DIAMOND output file: {output_path}")
        return output_path

    # Construct the path to the DIAMOND executable
    diamond_executable = os.path.join(diamond_path, "diamond")

    # Construct the DIAMOND blastp command
    command = [
        diamond_executable,
        "blastp",
        "-d", protein_db_path,
        "-q", transdecoder_pep,
        "-o", output_path
    ]

    # Log the process
    logging.info(f"Searching for matches in the protein database for {os.path.basename(transdecoder_pep)}...")
    error_msg = "DIAMOND blastp search failed. Check logs for details."
    run_subprocess(command, error_message=error_msg)
    logging.info(f"Blastp search completed successfully. Output written to {output_path}")

    return output_path



def get_optimized_pep_file(normal_pep, shortened_pep, classifications, output_path="revised_candidates.pep"):
    """
    Create a revised PEP file by incorporating classifications and shortened ORFs.

    This function reads a normal PEP file and a shortened PEP file, updates records based on
    classifications, and writes the revised records to a new file.

    :param normal_pep: Path to the original PEP file with all ORFs.
    :param shortened_pep: Path to the PEP file containing shortened ORFs.
    :param classifications: Dictionary of classifications for ORFs.
    :param output_path: Path to the output revised PEP file (default: "revised_candidates.pep").
    :return: Path to the output PEP file.
    """
    # Ensure the output directory exists
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    # Parse the shortened ORFs into a dictionary for quick access
    shortened_pep_dict = {
        record.id: (record.seq, record.description)
        for record in SeqIO.parse(shortened_pep, "fasta")
    }

    logging.info("Creating a revised PEP file using classifications...")

    with open(output_path, "w") as output:
        # Process each record in the normal PEP file
        for record in SeqIO.parse(normal_pep, "fasta"):
            orf_id = record.id

            # Determine the classification type from the record's description
            if "type:complete" in record.description:
                classification = "complete"
            elif "type:5prime_partial" in record.description:
                classification = "5prime_partial"
            elif "type:internal" in record.description:
                classification = "internal"
            elif "type:3prime_partial" in record.description:
                classification = "3prime_partial"
            else:
                classification = None

            # Handle incomplete ORFs (5prime_partial or internal)
            if classification in {"5prime_partial", "internal"}:
                if orf_id in classifications:
                    # Write the original record if classified as incomplete
                    if classifications[orf_id] == "incomplete":
                        SeqIO.write(record, output, "fasta")
                    else:
                        # Update the record with shortened sequence and adjusted description
                        seq, description = shortened_pep_dict[orf_id]
                        record.seq = seq
                        if "type:5prime_partial" in description:
                            description = description.replace("type:5prime_partial", "type:complete")
                        elif "type:internal" in description:
                            description = description.replace("type:internal", "type:3prime_partial")
                        record.description = description
                        SeqIO.write(record, output, "fasta")
                else:
                    # Write the original record if no classification is found
                    SeqIO.write(record, output, "fasta")
            else:
                # Write original records for complete or 3prime_partial ORFs
                SeqIO.write(record, output, "fasta")

    logging.info(f"Revised PEP file created successfully at: {output_path}")
    return output_path

def get_cds_classification(normal_tsv, shortened_tsv):
    """
    Compare DIAMOND search results to classify CDS candidates as complete or incomplete.

    This function parses the DIAMOND results from normal and shortened ORFs, calculates 
    support scores, and determines whether each CDS candidate should be classified as 
    "complete" or "incomplete" based on its best protein alignments.

    :param normal_tsv: Path to the DIAMOND results for the original TransDecoder output.
    :param shortened_tsv: Path to the DIAMOND results for the shortened ORFs.
    :return: Dictionary with CDS IDs as keys and classifications ("complete" or "incomplete") as values.
    """
    logging.info("Analyzing incomplete ORF predictions...")

    # Define DIAMOND output column names
    header_list = ["cdsID", "proteinID", "percIdentMatches", "alignLength", "mismatches", "gapOpenings",
                   "queryStart", "queryEnd", "targetStart", "targetEnd", "eValue", "bitScore"]

    # Load DIAMOND results into pandas DataFrames
    df_normal = pd.read_csv(normal_tsv, delimiter='\t', header=None, names=header_list)
    df_shortened = pd.read_csv(shortened_tsv, delimiter='\t', header=None, names=header_list)

    # Merge results on "cdsID" and "proteinID" to compare normal and shortened ORFs
    merged_df = pd.merge(df_shortened, df_normal, on=["cdsID", "proteinID"], suffixes=('_short', '_normal'))

    # Drop unnecessary columns
    columns_to_drop = ["alignLength_short", "alignLength_normal", "mismatches_short", "mismatches_normal",
                       "gapOpenings_short", "gapOpenings_normal", "queryEnd_short", "queryEnd_normal",
                       "targetEnd_short", "targetEnd_normal", "eValue_short", "eValue_normal"]
    merged_df.drop(columns=columns_to_drop, inplace=True)

    # Add support score column
    merged_df["supportScore"] = None

    # Compute support scores
    for i, row in merged_df.iterrows():
        q_incomplete_start = row["queryStart_normal"]
        t_incomplete_start = row["targetStart_normal"]
        t_complete_start = row["targetStart_short"]
        aai_incomplete = row["percIdentMatches_normal"]
        aai_complete = row["percIdentMatches_short"]

        # Prevent division by zero by assigning a small value if AAI is zero
        if aai_complete == 0:
            aai_complete = 0.0001

        match_log = math.log(aai_incomplete / aai_complete)
        support_score = (t_complete_start - t_incomplete_start) - (q_incomplete_start - 1) + match_log**1000
        merged_df.at[i, "supportScore"] = support_score

    # Compute maximum bit score for each candidate
    merged_df["bitScore_max"] = merged_df[["bitScore_short", "bitScore_normal"]].max(axis=1)

    # Group by cdsID and analyze the best 25 alignments
    classifications = {}
    for cds_id, group in merged_df.groupby("cdsID"):
        sorted_group = group.sort_values(by="bitScore_max", ascending=False)

        # Determine if any of the top 25 alignments has a positive support score
        incomplete = any(sorted_group.head(25)["supportScore"] > 0)

        classifications[cds_id] = "incomplete" if incomplete else "complete"

    logging.info("CDS classification completed successfully.")
    return classifications
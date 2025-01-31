import logging, sys, os, re, subprocess
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from utility import make_directory, run_subprocess, file_exists_and_not_empty

logger = logging.getLogger(__name__)

def getting_hc_supported_by_proteins(diamond_tsv, transdecoder_pep, protein_file, 
        output_path="hc_genes.pep"):
    """
    Select high-confidence genes supported by protein evidence.

    This function filters CDS candidates based on protein alignment criteria using DIAMOND results.
    Only complete CDS candidates that meet specific criteria are classified as high-confidence.

    :param diamond_tsv: Path to the DIAMOND BLASTp results file.
    :param transdecoder_pep: Path to the TransDecoder peptide FASTA file.
    :param protein_file: Path to the protein FASTA file used as the reference database.
    :param output_path: Path to the output file for high-confidence genes (default: "hc_genes.pep").
    :return: Dictionary of remaining CDS candidates that are not classified as high-confidence.
    """
    logging.info("Selecting high-confidence genes supported by protein evidence...")

    # Ensure the output directory exists
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    # Dictionary to store protein sequence lengths from the reference database
    t_length_dict = {
        record.id: len(record.seq)
        for record in SeqIO.parse(protein_file, "fasta")
    }

    # Dictionary to store query protein sequences and descriptions
    q_dict = {
        record.id: [record.description, record.seq]
        for record in SeqIO.parse(transdecoder_pep, "fasta")
    }

    # List to store the IDs of high-confidence genes
    already_hc_genes = set()

    # Process DIAMOND results and identify high-confidence genes
    with open(diamond_tsv, "r") as tsv, open(output_path, "w") as output:
        for line in tsv:
            part = line.strip().split("\t")
            query_id, protein_id = part[0], part[1]
            aaident = float(part[2])  # Percentage of identical amino acids
            align_length = int(part[3])  # Length of the alignment
            q_start, q_end = int(part[6]), int(part[7])  # Query alignment region
            t_start, t_end = int(part[8]), int(part[9])  # Target alignment region

            # Check if the query ID exists in q_dict
            if query_id in q_dict:
                record = SeqIO.SeqRecord(seq=q_dict[query_id][1], id=query_id, description=q_dict[query_id][0])

                # Calculate sequence lengths
                q_length = len(record.seq) - 1  # Exclude stop codon (*)
                t_length = t_length_dict.get(protein_id, 0)

                # Consider only complete candidates for high-confidence classification
                if "type:complete" in record.description:
                    if (q_start - t_start) < 6 and (t_length - align_length) < 15 and aaident > 95:
                        SeqIO.write(record, output, "fasta")
                        gene_id = ".".join(query_id.split(".")[:2])
                        already_hc_genes.add(gene_id)
                        del q_dict[query_id]  # Remove from further analysis

    # Remove CDS candidates belonging to genes already classified as high-confidence
    q_dict = {
        id: value
        for id, value in q_dict.items()
        if ".".join(id.split(".")[:2]) not in already_hc_genes
    }

    logging.info("High-confidence gene selection completed successfully.")
    return q_dict


def align_proteins(genome, protein, miniprot_path, miniprot_scorer_path, miniprothint_path,
        alignment_scoring=None, threads=4, output_dir="./", use_existing=False):
    """
    Align proteins to a reference genome using Miniprot and score alignments with Miniprot Boundary Scorer.

    This function first aligns protein sequences to a reference genome using Miniprot. Then, 
    it scores the alignment using Miniprot Boundary Scorer with the specified scoring matrix.

    :param genome: Path to the reference genome FASTA file.
    :param protein: Path to the protein FASTA file.
    :param miniprot_path: Path to the directory containing the Miniprot executable.
    :param miniprot_scorer_path: Path to the directory containing the Miniprot Boundary Scorer executable.
    :param miniprothint_path: Path to the directory containing the MiniprotHint executable.
    :param alignment_scoring: Path to the alignment scoring matrix.
    :param threads: Number of threads to use for Miniprot (default: 4).
    :param output_dir: Path to the output directory (default: "./").
    :param use_existing: Whether to use existing output files instead of rerunning (default: False).
    :return: Path to the scored alignment GFF file.
    """

    # Ensure output directory exists
    make_directory(output_dir)

    # Define output file paths
    output_aln = os.path.join(output_dir, "miniprot.aln")
    output_gff = os.path.join(output_dir, "miniprot_parsed.gff")
    output_gtf = os.path.join(output_dir, "miniprot.gtf")

    # Construct the Miniprot executable path
    miniprot_executable = os.path.join(miniprot_path, "miniprot")

    # Construct the Miniprot alignment command
    miniprot_command = [
        miniprot_executable,
        "-t", str(threads),
        "--genome", genome,
        "--protein", protein,
        "--aln"
    ]

    # miniprot_command = [
    #     miniprot_executable,
    #     "-t", str(threads),
    #     "-d", os.path.join(output_dir, "genome.mpi"), genome,
    # ]
    # run_subprocess(miniprot_command,  error_message='',)
    # miniprot_command = [
    #     miniprot_executable,'-I',
    #     "-ut", str(threads), "--outn=1", "--aln", os.path.join(output_dir, "genome.mpi"), 
    #     protein,        
    # ]
    logging.info("Aligning proteins to the genome using Miniprot...")
    if use_existing and file_exists_and_not_empty(output_aln):
        logging.info(f"Using existing alignment file: {output_aln}")
    else:
        error_msg = "Miniprot alignment failed. Check logs for details."
        with open(output_aln, "w") as output_file:
            run_subprocess(miniprot_command, stdout=output_file, error_message=error_msg,
                        capture_output=False)
        logging.info(f"Protein alignment completed successfully. Output written to {output_aln}")

    # Set default alignment scoring file if not provided
    if alignment_scoring is None:
        alignment_scoring = os.path.join(miniprot_scorer_path, "blosum62.csv")

    # Construct the Miniprot Boundary Scorer executable path
    miniprot_scorer_executable = os.path.join(miniprot_scorer_path, "miniprot_boundary_scorer")

    scorer_command = [miniprot_scorer_executable,
            "-o", output_gff, "-s", alignment_scoring]

    logging.info("Scoring the protein-to-genome alignment using Miniprot Boundary Scorer...")
    if use_existing and file_exists_and_not_empty(output_gff):
        logging.info(f"Using existing scored alignment file: {output_gff}")
    else:
        error_msg = "Miniprot Boundary Scorer failed. Check logs for details."
        with open(output_aln, 'r') as stdin_file:
            run_subprocess(scorer_command, error_message=error_msg, stdin=stdin_file)
            logging.info(f"Alignment scoring completed successfully. Output written to {output_gff}")


    # Convert scorer output to GTF using miniprothint
    miniprothint_executable = os.path.join(miniprothint_path, "miniprothint.py")

    cmd = [
        miniprothint_executable, output_gff, "--workdir", output_dir,
        "--ignoreCoverage", "--topNperSeed", "10", "--minScoreFraction", "0.5"
    ]
    
    logging.info("Convert scorer output to GTF using miniprothint...")
    if use_existing and file_exists_and_not_empty(output_gtf):
        logging.info(f"Using existing file: {output_gtf}")
    else:
        error_msg = "Miniprot failed. Check logs for details."
        with open(output_aln, 'r') as stdin_file:
            run_subprocess(cmd, error_message=error_msg)
            logging.info(f"Conversion completed successfully. Output written to {output_gtf}")

    return output_gff, output_gtf


def from_dict_to_pep_file(input_dict, output_path="output.pep"):
    """
    Convert a dictionary of sequences into a peptide FASTA file.

    This function takes a dictionary where keys are sequence IDs, and values are tuples
    containing sequence descriptions and Seq objects. It writes the sequences to a 
    FASTA file.

    :param input_dict: Dictionary containing sequence data {id: (description, SeqObject)}
    :param output_path: Path to the output peptide FASTA file (default: "output.pep").
    :return: Path to the generated FASTA file.
    """
    # Ensure the output directory exists
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    logging.info(f"Writing sequences to peptide FASTA file: {output_path}")

    # Convert dictionary to FASTA and write to output file
    with open(output_path, "w") as output_file:
        records = [
            SeqRecord(seq=seq, id=seq_id, description=desc)
            for seq_id, (desc, seq) in input_dict.items()
        ]
        SeqIO.write(records, output_file, "fasta")

    logging.info(f"Successfully created peptide FASTA file: {output_path}")
    return output_path

def choose_one_isoform(input_pep, output_pep):
    """
    Selects one isoform per gene based on length.

    If multiple isoforms exist, the longest is selected. If ties occur, the first 
    encountered isoform is chosen.

    :param input_pep: Path to the input PEP file.
    :param output_pep: Path to the output PEP file with one isoform per gene.
    """
    logging.info("Selecting one isoform per gene...")

    isoform_dict = {}

    # Parse sequences and group by gene ID
    for record in SeqIO.parse(input_pep, "fasta"):
        gene_id = record.id.split(".p")[0]
        cds_coords = re.search(r":(\d+)-(\d+)\((\+|-)\)", record.description)
        start, stop = int(cds_coords.group(1)), int(cds_coords.group(2))

        if gene_id not in isoform_dict:
            isoform_dict[gene_id] = [(start, stop, record)]
        else:
            isoform_dict[gene_id].append((start, stop, record))

    # Select and write longest isoform per gene
    with open(output_pep, "w") as output:
        for gene_id, isoforms in isoform_dict.items():
            record = max(isoforms, key=lambda x: x[1] - x[0])[2]  # Select longest
            SeqIO.write(record, output, "fasta")

    logging.info("Isoform selection complete. Output saved to %s", output_pep)

    return output_pep


def from_transcript_to_genome(cds_gff3, transcripts_gff3, transcripts_fasta, output_path, 
                        transdecoder_util):
    """
    Convert transcript coordinates to genome coordinates using TransDecoder.

    :param cds_gff3: Path to CDS annotation in GFF3 format in transcript coordinates.
    :param transcripts_gff3: Path to transcript annotation in GFF3 format in genome coordinates.
    :param transcripts_fasta: Path to the transcript sequences in FASTA format.
    :param output_path: Path to the output genome-coordinate GFF3 file.
    :param transdecoder_util: Path to Transdecoder utility scripts.
    """
    logging.info("Converting transcript coordinates to genome coordinates...")

    transdecoder_executable = os.path.join(transdecoder_util, "cdna_alignment_orf_to_genome_orf.pl")
    command = [
        transdecoder_executable,
        cds_gff3,
        transcripts_gff3,
        transcripts_fasta
    ]

    with open(output_path, "w") as output:
        run_subprocess(command, stdout=output, capture_output=False,
                error_message="TransDecoder coordinate transformation failed.")

    logging.info("Transcript-to-genome coordinate transformation complete.")

def preparing_candidates_for_conflict_comparison(candidates_gff3, transcripts_gtf, output_bed="candidates.bed"):
    """
    Prepare a BED file from CDS candidates to compare against Miniprot predictions.

    This function extracts chromosome information from the StringTie GTF file and converts 
    CDS candidate regions from the GFF3 file into a BED format.

    :param candidates_gff3: Path to the CDS candidates in GFF3 format.
    :param transcripts_gtf: Path to the StringTie transcript annotations (GTF).
    :param output_bed: Path to the output BED file (default: "candidates.bed").
    :return: Path to the generated BED file.
    """
    logging.info("Preparing BED file for CDS candidates conflict comparison...")

    # Extract chromosome information from the GTF file
    chromosome_dict = {}
    with open(transcripts_gtf, "r") as transcripts_file:
        for line in transcripts_file:
            if line.startswith("#"):
                continue
            parts = line.split("\t")
            if parts[2] == "transcript":
                gene_id_match = re.search(r'gene_id "([^"]+)"', line)
                if gene_id_match:
                    gene_id = gene_id_match.group(1)
                    chromosome_dict[gene_id] = parts[0]

    # Convert GFF3 CDS records to BED format
    with open(output_bed, "w") as output, open(candidates_gff3, "r") as candidates_file:
        for line in candidates_file:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split("\t")
            if parts[2] == "CDS":
                start, stop = parts[3], parts[4]
                gene_id = re.search(r"(STRG\.\d+)", line).group(1)
                transcript_id = re.search(r'Parent=([^\s;]+)', line).group(1)
                chromosome = chromosome_dict.get(gene_id, "unknown")
                strand = parts[6]

                output.write(f"{chromosome}\t{start}\t{stop}\t{transcript_id}\t.\t{strand}\n")

    logging.info("BED file for conflict comparison created: %s", output_bed)
    return output_bed


def finding_protein_conflicts(candidates_bed, reference_bed, bedtools_path, output_conflicts="conflicts.bed"):
    """
    Identify overlaps between CDS candidates and Miniprot gene predictions using Bedtools.

    :param candidates_bed: Path to the BED file of candidate CDS regions.
    :param reference_bed: Path to the BED file of Miniprot-predicted regions.
    :param bedtools_path: Path to the directory containing the Bedtools executable.
    :param output_conflicts: Path to the output BED file with conflict results (default: "conflicts.bed").
    :return: Path to the generated conflict BED file.
    """
    logging.info("Identifying overlaps between CDS candidates and Miniprot predictions using Bedtools...")

    # Construct the Bedtools command
    bedtools_executable = os.path.join(bedtools_path, "bedtools")
    command = [
        bedtools_executable,
        "coverage",
        "-a", candidates_bed,
        "-b", reference_bed,
        "-s"
    ]

    with open(output_conflicts, "w") as output_file:
        run_subprocess(command, stdout=output_file, error_message="Bedtools conflict detection failed.",
            capture_output=False)    

    logging.info("Conflicts identified successfully. Output saved to %s", output_conflicts)
    return output_conflicts

def finding_stop_in_utr(transcripts_fasta, intrinsic_candidates_gff3):
    """
    Identify CDS candidates with stop codons in the 5' UTR.

    :param transcripts_fasta: Path to the transcript sequences in FASTA format.
    :param intrinsic_candidates_gff3: Path to the candidate CDS regions in GFF3 format.
    :return: Dictionary {transcript_id: True} for candidates with stop codons in the 5' UTR.
    """
    logging.info("Searching for stop codons in the 5' UTR of candidates...")

    stop_codons = {"TAA", "TAG", "TGA"}
    stop_in_utr_dict = {}

    # Load transcript sequences into a dictionary
    transcripts_dict = {
        record.id: str(record.seq) for record in SeqIO.parse(transcripts_fasta, "fasta")
    }

    # Process GFF3 file to find 5' UTRs
    with open(intrinsic_candidates_gff3, "r") as candidates_file:
        for line in candidates_file:
            if line.startswith("#") or not line.strip():
                continue

            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue  # Skip malformed lines

            feature_type, start, end, attributes = parts[2], int(parts[3]), int(parts[4]), parts[8]
            transcript_match = re.search(r"ID=([^;]+?)\.p\d+\b", attributes)
            if transcript_match:
                transcript_id = transcript_match.group(1)

            # Check for stop codons in 5' UTR
            if feature_type == "five_prime_UTR" and transcript_id in transcripts_dict:
                utr_sequence = transcripts_dict[transcript_id][: (end - start + 1)]
                if any(codon in utr_sequence for codon in stop_codons):
                    stop_in_utr_dict[transcript_id] = True

    logging.info("Identified stop codons in the 5' UTR for %d candidates.", len(stop_in_utr_dict))
    return stop_in_utr_dict


def from_pep_file_to_gff3(orf_pep, transcript_gtf, output_gff3):
    """
    Convert a TransDecoder peptide FASTA file and StringTie transcript annotations (GTF) into a GFF3 file.

    This function calculates transcript lengths from the GTF file and generates a GFF3 file containing 
    mRNA, exon, CDS, and UTR annotations.

    :param orf_pep: Path to the TransDecoder peptide FASTA file.
    :param transcript_gtf: Path to the StringTie transcript annotation file (GTF).
    :param output_gff3: Path to the output GFF3 file.
    :return: Path to the generated GFF3 file.
    """
    logging.info("Generating GFF3 file from TransDecoder PEP and StringTie GTF...")

    # Step 1: Parse transcript lengths from GTF
    transcript_lengths = {}

    with open(transcript_gtf, "r") as transcript_file:
        for line in transcript_file:
            if line.startswith("#"):
                continue

            parts = line.strip().split("\t")
            if parts[2] == "exon":
                start, stop = int(parts[3]), int(parts[4])
                length = stop - start + 1
                transcript_id_match = re.search(r'transcript_id "([^"]+)"', parts[8])

                if transcript_id_match:
                    transcript_id = transcript_id_match.group(1)
                    transcript_lengths[transcript_id] = transcript_lengths.get(transcript_id, 0) + length

    # Step 2: Create the GFF3 file
    tool_stringtie = "StringTie"
    tool_transdecoder = "TransDecoder"

    with open(output_gff3, "w") as output:
        for record in SeqIO.parse(orf_pep, "fasta"):
            transcript_id = record.id.split(".p")[0]

            # if transcript_id not in transcript_lengths:
            #     logging.warning(f"Transcript ID {transcript_id} not found in GTF, skipping.")
            #     continue

            transcript_length = transcript_lengths[transcript_id]
            description = record.description

            # Extract gene ID, gene name, ORF coordinates, and strand
            # gene_id_match = re.search(r"gene=([^ ]+)", description)
            # gene_name_match = re.search(r"Name=([^ ]+)", description)
            orf_coords_match = re.search(r":(\d+)-(\d+)\((\+|-)\)", description)

            if not orf_coords_match:
                logging.warning(f"Skipping {record.id} due to missing annotation.")
                continue

            description_parts = description.split()
            gene_id = description_parts[0]
            gene_name = description_parts[1]
            # gene_name = gene_name_match.group(1) if gene_name_match else "unknown_gene"
            orf_start, orf_stop, strand = int(orf_coords_match.group(1)), int(orf_coords_match.group(2)), orf_coords_match.group(3)

            # Write GFF3 records
            output.write(f"{transcript_id}\t{tool_stringtie}\tgene\t1\t{transcript_length}\t.\t{strand}\t.\tID={gene_id};Name={gene_name}\n")
            output.write(f"{transcript_id}\t{tool_stringtie}\tmRNA\t1\t{transcript_length}\t.\t{strand}\t.\tID={record.id};Parent={gene_id};Name={gene_name}\n")
            output.write(f"{transcript_id}\t{tool_stringtie}\texon\t1\t{transcript_length}\t.\t{strand}\t.\tID={record.id}.exon1;Parent={record.id}\n")
            output.write(f"{transcript_id}\t{tool_transdecoder}\tCDS\t{orf_start}\t{orf_stop}\t.\t{strand}\t0\tID=cds.{record.id};Parent={record.id}\n")

            # Add UTR annotations if applicable
            if strand == "+":
                if orf_start > 1:
                    output.write(f"{transcript_id}\t{tool_transdecoder}\tfive_prime_UTR\t1\t{orf_start-1}\t.\t{strand}\t.\tID={record.id}.utr5p1;Parent={record.id}\n")
                if orf_stop < transcript_length:
                    output.write(f"{transcript_id}\t{tool_transdecoder}\tthree_prime_UTR\t{orf_stop+1}\t{transcript_length}\t.\t{strand}\t.\tID={record.id}.utr3p1;Parent={record.id}\n")
            elif strand == "-":
                if orf_stop < transcript_length:
                    output.write(f"{transcript_id}\t{tool_transdecoder}\tfive_prime_UTR\t{orf_stop+1}\t{transcript_length}\t.\t{strand}\t.\tID={record.id}.utr5p1;Parent={record.id}\n")
                if orf_start > 1:
                    output.write(f"{transcript_id}\t{tool_transdecoder}\tthree_prime_UTR\t1\t{orf_start-1}\t.\t{strand}\t.\tID={record.id}.utr3p1;Parent={record.id}\n")

            output.write("\n")

    logging.info("GFF3 file created successfully: %s", output_gff3)
    return output_gff3

def preparing_miniprot_gff_for_conflict_comparison(miniprot_gff, output_bed="reference.bed"):
    """
    Convert Miniprot gene predictions (GFF) into a BED format for conflict comparison.

    This function extracts CDS entries from the Miniprot GFF file and converts them into a BED format 
    for comparison with TransDecoder predictions.

    :param miniprot_gff: Path to the Miniprot GFF file.
    :param output_bed: Path to the output BED file (default: "reference.bed").
    :return: Path to the generated BED file.
    """
    logging.info("Converting Miniprot GFF to BED format for conflict comparison...")

    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_bed), exist_ok=True)

    # Process GFF and extract CDS features
    with open(miniprot_gff, "r") as gff_file, open(output_bed, "w") as output:
        for line in gff_file:
            if line.startswith("#") or not line.strip():
                continue

            parts = line.strip().split("\t")
            if len(parts) < 9 or parts[2] != "CDS":
                continue  # Skip non-CDS entries

            chromosome, start, stop, strand = parts[0], parts[3], parts[4], parts[6]
            prot_id_match = re.search(r"prot=(prot\d+)", parts[8])

            if prot_id_match:
                prot_id = prot_id_match.group(1)
                output.write(f"{chromosome}\t{start}\t{stop}\t{prot_id}\t.\t{strand}\n")
            else:
                logging.warning(f"Skipping line due to missing protein ID: {line.strip()}")

    logging.info("Miniprot BED file created successfully: %s", output_bed)
    return output_bed


def getting_hc_supported_by_intrinsic(q_dict, stringtie_gtf, stringtie_gff3, 
                    transcript_fasta, miniprot_gff, transdecoder_util_path, bedtools_path,
                    output_dir="./", min_cds_len=300, score_threshold=50):
    """
    Identify and classify high-confidence (HC) genes supported by intrinsic criteria.

    This function processes candidate coding sequences (CDS) and applies a series of filtering 
    criteria to classify them as either high-confidence (HC) or low-confidence (LC) genes. The 
    filtering is based on:
      - Minimum coding sequence (CDS) length
      - Alignment score threshold
      - Conflict status with Miniprot predictions
      - Presence of a stop codon in the 5' UTR
      - Completion status of the CDS (only complete ORFs are considered)

    The function generates multiple intermediate files including PEP, GFF3, and BED files for 
    conflict detection. The final output consists of two PEP files:
      - `hc_genes.pep`: High-confidence genes
      - `lc_genes.pep`: Low-confidence genes

    :param q_dict: Dictionary containing candidate ORFs.
    :param stringtie_gtf: Path to the StringTie transcript annotation file.
    :param transdecoder_gff: Path to the TransDecoder prediction of transcriptomic data.
    :param transcript_fasta: Path to the transcripts FASTA file.
    :param miniprot_gff: Path to the Miniprot predictions file.
    :param transdecoder_util_path: Path to TransDecoder utilities for coordinate transformations.
    :param bedtools_path: Path to the Bedtools executable directory for conflict detection.
    :param output_dir: Directory to store the output files (default: "./").
    :param min_cds_len: Minimum length required for a coding sequence (default: 300).
    :param score_threshold: Minimum required alignment score for high-confidence classification (default: 50).
    
    :return: Paths to the generated high-confidence (HC) and low-confidence (LC) gene PEP files.
    """
    make_directory(output_dir)

    logging.info("Selecting high-confidence genes supported by intrinsic criteria...")

    # files to be created
    intrinsic_pep = os.path.join(output_dir, "intrinsic_candidates.pep")
    intrinsic_one_isoform = os.path.join(output_dir, "intrinsic_one_isoform.pep")
    intrinsic_gff3 = os.path.join(output_dir, "intrinsic_candidates.gff3")
    intrinsic_genome_gff3 = os.path.join(output_dir, "intrinsic_candidates_genome.gff3")
    candidates_bed = os.path.join(output_dir, "candidates.bed")
    reference_bed = os.path.join(output_dir, "reference.bed")
    conflicts_bed = os.path.join(output_dir, "conflicts.bed")
    hc_genes_pep = os.path.join(output_dir, "hc_genes.pep")
    lc_genes_pep = os.path.join(output_dir, "lc_genes.pep")

    # Generate candidate PEP file
    from_dict_to_pep_file(q_dict, intrinsic_pep)

    # Select one isoform per gene
    choose_one_isoform(intrinsic_pep, intrinsic_one_isoform)

    # Convert PEP to GFF3
    from_pep_file_to_gff3(intrinsic_one_isoform, stringtie_gtf, intrinsic_gff3)

    # Convert transcript coordinates to genome coordinates
    from_transcript_to_genome(intrinsic_gff3, stringtie_gff3, transcript_fasta, 
        intrinsic_genome_gff3, transdecoder_util_path)

    # Prepare conflict comparison files
    preparing_candidates_for_conflict_comparison(intrinsic_genome_gff3, stringtie_gtf, candidates_bed)
    preparing_miniprot_gff_for_conflict_comparison(miniprot_gff, reference_bed)

    #Find conflicts between the candidates and the miniprot predictions.
    finding_protein_conflicts(candidates_bed, reference_bed, bedtools_path, conflicts_bed)

    #Find candidates with a stop codon in the 5' UTR.
    stop_in_utr_dict = finding_stop_in_utr(transcript_fasta, intrinsic_gff3)

    # Process conflicts
    no_conflicts_dict = {}
    with open(conflicts_bed, "r") as conflicts_file:
        for line in conflicts_file:
            parts = line.strip().split("\t")
            transcript_id = parts[3]
            no_conflicts_dict[transcript_id] = int(parts[6]) == 0  # True if no conflicts

    # Filter high-confidence candidates
    already_hc_genes = set()
    with open(hc_genes_pep, "a") as output_hc, open(intrinsic_one_isoform, "r") as candidates_file:
        for record in SeqIO.parse(candidates_file, "fasta"):
            transcript_id = record.id
            length_pep = int(re.search(r"len:(\d+)", record.description).group(1))
            length_cds = length_pep * 3
            score = float(re.search(r"score=(-?[\d.]+)", record.description).group(1))
            id_stringtie = transcript_id.split(".p")[0]

            conditions = [
                length_cds >= min_cds_len,  # Min length
                score > score_threshold,  # Score threshold
                no_conflicts_dict.get(transcript_id, False),  # No conflicts
                "type:complete" in record.description,  # Complete ORFs
                stop_in_utr_dict.get(id_stringtie, False)  # Stop codon in 5' UTR
            ]

            if all(conditions):
                already_hc_genes.add(id_stringtie)
                SeqIO.write(record, output_hc, "fasta")

    logging.info("High-confidence genes identified and saved.")

    # Identify low-confidence genes
    with open(lc_genes_pep, "w") as output_lc, open(intrinsic_one_isoform, "r") as candidates_file:
        for record in SeqIO.parse(candidates_file, "fasta"):
            transcript_id = record.id
            id_stringtie = transcript_id.split(".p")[0]
            if id_stringtie not in already_hc_genes:
                SeqIO.write(record, output_lc, "fasta")

    logging.info("Low-confidence genes identified and saved.")

    return hc_genes_pep, lc_genes_pep
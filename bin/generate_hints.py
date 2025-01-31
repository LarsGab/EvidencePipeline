import os, logging, shutil
from utility import run_subprocess

logger = logging.getLogger(__name__)

def aln2hints(aln_file, final_out_file, scripts_path, prg="miniprot", priority=4):
    """
    Convert protein alignment file to AUGUSTUS-compatible hints.

    This function runs the `aln2hints.pl` script to convert alignments into
    hints for gene prediction with AUGUSTUS.

    :param aln_file: Path to the input alignment file.
    :param final_out_file: Path to the final hints output file.
    :param scripts_path: Path to the scripts directory with aln2hints.pl.
    :param prg: Alignment program used, either "gth" or "miniprot" (default: "miniprot").
    :param priority: Priority level for the hints (default: 4).
    """

    logging.info(f"Converting alignments from {aln_file} to AUGUSTUS hints...")

    # Check if the alignment file is empty
    if os.path.isfile(aln_file) and os.path.getsize(aln_file) > 0:
        out_file_name = os.path.join(os.path.dirname(final_out_file), "prot_hintsfile.aln2hints.temp.gff")
        aln2hints_executable = os.path.join(scripts_path, "aln2hints.pl")

        # Construct the command
        command = [
            "perl", aln2hints_executable,
            f"--in={aln_file}",
            f"--out={out_file_name}",
            f"--prg={prg}",
            f"--priority={priority}"
        ]
        # Run the aln2hints command
        run_subprocess(command, error_message=f"Failed to convert {aln_file} to hints.")

        # Concatenate the temporary hints file to the final output
        with open(final_out_file, 'a+') as out_file, open(out_file_name, "r") as in_file:
            for line in in_file.readlines():
                out_file.write(line)

        logging.info(f"Protein hints successfully generated and saved to {final_out_file}.")

    else:
        logging.warning(f"Alignment file {aln_file} was empty. No hints generated.")

def bam2hints(bam_file, genome, out_file, hints_file, bam2hints_path, scripts_dir):
    """
    Convert bam file to AUGUSTUS hints.

    :param bam_file: Path to the input BAM file.
    :param out_file: Path to the final hints output file.
    :param bam2hints_path: Path to the scripts directory with bam2hints.
    """

    logging.info(f"Converting alignments from {bam_file} to AUGUSTUS hints...")

    # Check if the alignment file is empty
    if os.path.isfile(bam_file) and os.path.getsize(bam_file) > 0:        
        bam2hints_executable = os.path.join(bam2hints_path, "bam2hints")

        # Construct the command
        command = [
            bam2hints_executable,
            '--intronsonly',
            f"--in={bam_file}",
            f"--out={out_file}.temp",
        ]

        # Run the bam2hints command
        run_subprocess(command, error_message=f"Failed to convert {bam_file} to hints.")

        logging.info(f"Transcriptome hints successfully generated and saved to {out_file}.temp.")

        # filterIntronsFindStrand.pl
        logging.info(f"Filter intron hints and add strand: {bam_file}.temp .")
        filter_executable = os.path.join(scripts_dir, "filterIntronsFindStrand.pl")

        # Construct the command
        command = [
            filter_executable,
            genome,
            f"{out_file}.temp",
            '--score'
        ]
        with open(out_file, 'w+') as out_f:
            run_subprocess(command, error_message=f"Failed to add strand to {out_file}.temp.",
                    stdout=out_f,  capture_output=False,)

        # Concatenate the temporary hints file to the final output
        with open(hints_file, 'a+') as h_file, open(out_file, "r") as in_file:
            for line in in_file.readlines():
                h_file.write(line)
    else:
        logging.warning(f"BAM file {bam_file} was empty. No hints generated.")


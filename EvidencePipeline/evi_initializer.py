import argparse, logging, shutil, os
from EvidencePipeline.utility import read_yaml, find_executable, file_exists_and_not_empty, make_absolute_path

#Authors: "Amrei Knuth", "Lars Gabriel"
#Credits: "Katharina Hoff"
#Email:"lars.gabriel@uni-greifswald.de"
#Date: "Janurary 2025"

logger = logging.getLogger(__name__)


def get_config(args):
    """Read argument and generate a config dictionary"""
    if args.config:
        cfg = read_yaml(args.config)
    else:
        cfg = {}
    # Overwrite values in dict with args
    for key in vars(args):
        value = getattr(args, key)
        if value is not None:  # Only overwrite if the argument is provided            
            if key in ["rnaseq_paired", "rnaseq_single", "isoseq"]:
                values = value.split(",") 
            cfg[key] = value
            if value:
                logging.info(f'Setting {key} to {value}.')

    for key in ["genome", "proteins", "scoring_matrix", "output_path"]:
        if key in cfg and cfg[key]:
            cfg[key] = make_absolute_path(cfg[key])
    for key in ["rnaseq_paired", "rnaseq_single", "isoseq"]:
        if key in cfg and cfg[key]:
            cfg[key] = [make_absolute_path(v) for v in cfg[key]]
    return cfg

def determine_mode(cfg):
    """Determine the mode based on the given input data.
    
    :cfg: Dictionary containing input configuration data.
    :return: The determined mode as a string, or raises an error if requirements are not met.
    """
    # Define mode requirements
    mode_requirements = {'mixed' : [['rnaseq_paired', 'rnaseq_single'], ['isoseq', 'proteins']],
                'isoseq': [[], ['isoseq', 'proteins']], 
                'rnaseq': [['rnaseq_paired', 'rnaseq_single'], ['proteins']],
                'proteins': [[], ['proteins']]}
    for mode, (any_keys, all_keys) in mode_requirements.items():
        # Check if at least one key from `any_keys` exists and has a non-None value
        any_valid = not any_keys or any(
            key in cfg and cfg[key] not in [None, []] for key in any_keys
        )
        
        # Check if all keys from `all_keys` exist and have non-None values
        all_valid = all(
            key in cfg and cfg[key] not in [None, []] for key in all_keys
        )

        # If both conditions are satisfied, return the mode
        if any_valid and all_valid:
            logging.info(f'Running in {mode} mode.')
            return mode

    # If no mode matches, raise an error
    error_message = (
        "Input missing! You need to specify at least 'proteins' and either Isoseq or RNA-Seq data."
    )
    logging.error(error_message)
    raise ValueError(error_message)


def check_dependencies(cfg):
    """Try to find required dependencies for mode"""
    cfg["GALBA_SCRIPTS"] = find_executable(
                tool_path=cfg['GALBA_SCRIPTS'] if cfg['GALBA_SCRIPTS'] else os.path.dirname(os.path.realpath(__file__)), 
                executable_name="aln2hints.pl", tool_name='GALBA_SCRIPTS')
    cfg["DIAMOND"] = find_executable(tool_path=cfg['DIAMOND'], 
                executable_name="diamond", tool_name='DIAMOND')
    cfg["BEDTOOLS"] = find_executable(tool_path=cfg['BEDTOOLS'], 
                executable_name="bedtools", tool_name='BEDTOOLS')
    cfg["MINIPROT"] = find_executable(tool_path=cfg['MINIPROT'], 
                executable_name="miniprot", tool_name='MINIPROT')
    cfg["MINIPROT_HINT"] = find_executable(tool_path=cfg['MINIPROT_HINT'], 
                executable_name="miniprothint.py", tool_name='MINIPROT_HINT')
    cfg["MINIPROT_BOUNDARY_SCORER"] = find_executable(tool_path=cfg['MINIPROT_BOUNDARY_SCORER'], 
                executable_name="miniprot_boundary_scorer", tool_name='MINIPROT_BOUNDARY_SCORER')

    # add default scoring matrix for miniprot-boundary-scorer
    if not "scoring_matrix" in cfg or cfg["scoring_matrix"] is None:
        cfg["scoring_matrix"] = cfg["MINIPROT_BOUNDARY_SCORER"] + "/blosum62.csv"
        if not file_exists_and_not_empty(cfg["scoring_matrix"]):
            error_message = f"Couldn't find scoring matrix for miniprot-boundary-scorer. It should be located in the directory of miniprot-boundary-scorer or specified as argument."
            logging.error(error_message)
            raise ValueError(error_message)

    if cfg["mode"] in ['mixed', 'rnaseq']:
        cfg["HISAT"] = find_executable(tool_path=cfg['HISAT'], 
                    executable_name="hisat2", tool_name='HISAT')

    if cfg["mode"] in ['mixed', 'isoseq']:
        cfg["MINIMAP"] = find_executable(tool_path=cfg['MINIMAP'], 
                    executable_name="minimap2", tool_name='MINIMAP')

    if cfg["mode"] in ['mixed', 'isoseq', 'rnaseq']:
        cfg["STRINGTIE"] = find_executable(tool_path=cfg['STRINGTIE'], 
                    executable_name="stringtie", tool_name='STRINGTIE')
        cfg["SAMTOOLS"] = find_executable(tool_path=cfg['SAMTOOLS'], 
                    executable_name="samtools", tool_name='SAMTOOLS')
        cfg["TRANSDECODER"] = find_executable(tool_path=cfg['TRANSDECODER'], 
                    executable_name="TransDecoder.LongOrfs", tool_name='TRANSDECODER')
        cfg["TRANSDECODER_UTIL"] = find_executable(tool_path=f"{cfg['TRANSDECODER']}/util", 
                    executable_name="gtf_genome_to_cdna_fasta.pl", tool_name='TRANSDECODER')
        cfg["BAM2HINTS"] = find_executable(tool_path=cfg['BAM2HINTS'], 
                    executable_name="bam2hints", tool_name='BAM2HINTS')
    return cfg


def setup_logger(logfile=None):
    logging.basicConfig(
        level=logging.DEBUG,  # Change to INFO, WARNING, etc., as needed
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        handlers=[
            logging.FileHandler(logfile),  # Log to a file
            logging.StreamHandler(),        # Log to the console
        ],
    )

def get_parser():
    """Define the command-line arguments accepted by the script."""
    
    parser = argparse.ArgumentParser(description='Genome annotation with protein and transcriptomic evidence auch as RNA-Seq and Iso-Seq reads')  
    parser.add_argument('-t', '--threads', default=4, help='Number of threads (default=4)', required=False)
    parser.add_argument('-y', '--config', help='Config file path', metavar='<config.yaml>', required=False) 
    parser.add_argument('-g', '--genome', help='Genome file path', metavar='<genome.fasta>', required=False)
    parser.add_argument('-p', '--proteins', help='Protein file path', metavar='<proteins.fasta>', required=False)
    parser.add_argument('-sp', '--rnaseq_paired', help='Comma separated paired-end short read file paths', 
            metavar='<set1.fasta,set2.fasta>', required=False, default=None)
    parser.add_argument('-ss', '--rnaseq_single', help='Comma separated single-end short read file paths', 
            metavar='<set1.fasta,set2.fasta>', required=False, default=None)
    parser.add_argument('-l', '--isoseq', help='Comma separated long read file paths', 
            metavar='<set1.fasta,set2.fasta>', required=False, default=None)
    parser.add_argument('-m', '--scoring_matrix', help='Alignment scoring matrix for amino acids', metavar='<blosum62.csv>', required=False)
    parser.add_argument('--output_path',  default='./pregalba/',
            help='Specify the path you want the output folder to be created in if it should differ from the current working directory')
    parser.add_argument('--restart', action='store_true', 
            help='Use results from a previous run to skip steps. The results have to be located in --output_path.')
    parser.add_argument('--HISAT', help='Provide the path to the HISAT2 directory', required=False, default='')
    parser.add_argument('--MINIMAP', help='Provide the path to the Minimap2 directory', required=False, default='')
    parser.add_argument('--STRINGTIE', help='Provide the path to the StringTie2 directory', required=False, default='')
    parser.add_argument('--SAMTOOLS', help='Provide the path to the samtools directory', required=False, default='')
    parser.add_argument('--TRANSDECODER', help='Provide the path to the TransDecoder directory', required=False, default='')
    parser.add_argument('--DIAMOND', help='Provide the path to the DIAMOND directory', required=False, default='')
    parser.add_argument('--BEDTOOLS', help='Provide the path to the bedtools2 directory', required=False, default='')
    parser.add_argument('--MINIPROT', help='Provide the path to the miniprot directory', required=False, default='')
    parser.add_argument('--MINIPROT_HINT', help='Provide the path to the miniprot_hint directory', required=False, default='')
    parser.add_argument('--MINIPROT_BOUNDARY_SCORER', help='Provide the path to the miniprot-boundary-scorer directory', required=False, default='')
    parser.add_argument('--GALBA_SCRIPTS', help='Provide the path to galba scripts like aln2hints.pl directory. Default is the same directory as pregalba.py', required=False, default='')
    parser.add_argument('--BAM2HINTS', help='Provide the path to bam2hints directory.', required=False, default='')
    args = parser.parse_args()
    return args
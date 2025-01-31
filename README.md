# EvidencePipeline

EvidencePipeline is a tool for generating evidence for genome annotation using protein and transcriptomic data (RNA-Seq and Iso-Seq reads). It automates alignment, evidence generation, and high-confidence gene selection.

## Installation

### Prerequisites
Ensure that you have the following dependencies installed:
- [HISAT2](https://daehwankimlab.github.io/hisat2/)
- [Minimap2](https://github.com/lh3/minimap2)
- [StringTie](https://ccb.jhu.edu/software/stringtie/)
- [SAMtools](http://www.htslib.org/)
- [TransDecoder](https://github.com/TransDecoder/TransDecoder)
- [DIAMOND](https://github.com/bbuchfink/diamond)
- [Bedtools](https://bedtools.readthedocs.io/)
- [Biopython](https://biopython.org/)

The relevant scripts should be in your PATH or specified as command line arguments.

### Install Python Dependencies
```bash
pip install -r requirements.txt
```

## Usage

### Running EvidencePipeline

To run EvidencePipeline, use the following command:
```bash
python evidence_pipeline.py -g <genome.fasta> -p <proteins.fasta> [-sp <paired_reads>] [-ss <single_reads>] [-l <long_reads>] -y <config.yaml> --output_path <output_dir>
```

### Required Arguments:
- `-g, --genome` : Path to the genome FASTA file.
- `-p, --proteins` : Path to the protein FASTA file.
OR
- `-y, --config` : Path to a configuration YAML file that includes genome and proteins path. If config and command-line arguments are both provided, command-line arguments will take precedence.

### Optional Arguments:
- `-sp, --rnaseq_paired` : Comma-separated paired-end RNA-Seq read files.
- `-ss, --rnaseq_single` : Comma-separated single-end RNA-Seq read files.
- `-l, --isoseq` : Comma-separated long read files.
- `--output_path` : Directory where output files will be stored.
- `--restart` : Resume from previous results if available.
- `--threads` : Number of threads to use (default: 4).

### Example Run
```bash
python evidence_pipeline.py -g genome.fasta -p proteins.fasta -sp rnaseq_1.fq,rnaseq_2.fq -l isoseq.fq -y config.yaml --output_path results/
```

## Pipeline Overview
1. **Initialization:**
   - Reads configuration parameters.
   - Determines the mode (RNA-Seq, Iso-Seq, or Mixed).
   - Checks dependencies.

2. **Alignment:**
   - Align proteins to the genome using Miniprot.
   - Align RNA-Seq reads using HISAT2.
   - Align Iso-Seq reads using Minimap2.

3. **Transcript Assembly & ORF Prediction:**
   - Uses StringTie to assemble transcripts.
   - Uses TransDecoder to predict open reading frames (ORFs).

4. **Evidence Generation:**
   - Converts protein and transcript alignments into AUGUSTUS hints.
   - Generates high-confidence genes.

5. **Output Generation:**
   - Produces GFF3, FASTA, and other annotation files.

## Output Files
- `hintsfile.gff` : AUGUSTUS-compatible hints file.
- `training.gff` : High-confidence genes that can be used for training.

## Logging
All logs are saved in `LOG.out` inside the output directory.

## License
This project is licensed under the MIT License.

## Contact
For questions or issues, please open an issue on the repository or contact the developer.


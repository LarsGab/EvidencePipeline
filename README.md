# EvidencePipeline

EvidencePipeline is a Nextflow workflow for generating high-confidence genes and genome annotation hints from extrinsic evidence.

Three sources of extrinsic evidence can be provided:
- Protein sequences
- RNA-Seq short reads
- Iso-Seq long reads

Additionallty, Tiberius can be used to make *ab initio* predictions. These are then merged with the high-confidence genes for a genome-wide genome annotation. 

## Requirements
- [Nextflow](https://www.nextflow.io/) and Java 11+.
- Singularity in order to not require any other dependencies.
- A nextflow config for your HPC scheduler or local setting (documentaion and examples in `conf/`).


## Installation without Singularity Container

**We strongly recommend to use the available Singularity container.** Is this not possible, you can install all required tools manually on your system. Following tools are required and there executables must be in your PATH variable (or you have to set their paths in `conf/base.config`):

| Tool             | Description                        | Installation Page                                                                            |
| ---------------- | ---------------------------------- | -------------------------------------------------------------------------------------------- |
| **Minimap2**     | Used for alignment of Iso seq        | [https://github.com/lh3/minimap2](https://github.com/lh3/minimap2)                           |
| **StringTie**    | Transcript assembler               | [https://github.com/gpertea/stringtie](https://github.com/gpertea/stringtie)                 |
| **Samtools**     | BAM utilities                 | [http://www.htslib.org/](http://www.htslib.org/)                                             |
| **TransDecoder** | ORF prediction                     | [https://github.com/TransDecoder/TransDecoder](https://github.com/TransDecoder/TransDecoder) |
| **DIAMOND**      | Fast protein aligner               | [https://github.com/bbuchfink/diamond](https://github.com/bbuchfink/diamond)                 |
| **BEDTools**     | Genomic interval operations        | [https://github.com/arq5x/bedtools2](https://github.com/arq5x/bedtools2)                     |
| **BAM2HINTS**     | Generates genome annotation hints       | [https://github.com/Gaius-Augustus/Augustus](https://github.com/Gaius-Augustus/Augustus)     |
| **MiniProt**                 | Protein-to-genome aligner               | [https://github.com/lh3/miniprot](https://github.com/lh3/miniprot)                                               |
| **MiniProtHint**             | Improved protein hints for AUGUSTUS     | [https://github.com/tomasbruna/miniprothint](https://github.com/tomasbruna/miniprothint)                         |
| **MiniProt Boundary Scorer** | Boundary scoring for protein alignments | [https://github.com/tomasbruna/miniprot-boundary-scorer](https://github.com/tomasbruna/miniprot-boundary-scorer) |
| **HISAT2** (v2.2.1)          | Splice-aware short-read aligner         | [https://daehwankimlab.github.io/hisat2/](https://daehwankimlab.github.io/hisat2/)                               |
| **BEDTools**                 | Interval operations                     | [https://github.com/arq5x/bedtools2](https://github.com/arq5x/bedtools2)                                         |
| **Prefetch / fasterq-dump**  | Part of SRA-Tools                       | [https://github.com/ncbi/sra-tools](https://github.com/ncbi/sra-tools)                                           |



## Configuration

See `conf/README.md`


## Running the Evidence Pipeline

You can start the nextflow pipeline with `evidence_pipeline.py`:

```bash
cd EvidencePipeline
./evidence_pipeline.py \
  --params_yaml conf/params.yaml \
  --nf_config conf/slurm_generic.config 
```

The wrapper checks that
1. the provided params YAML and Nextflow config exist,
2. all referenced input data files are available, and
3. the required executables (Nextflow, Singularity, HISAT2, Minimap2, TransDecoder, MiniProt, etc.) can be located on your PATH.

Use `--dry_run` to only perform validations, `--check_tools` when running with native tool installations instead of the Singularity container, and pass additional Nextflow options after `--`, e.g. `./run_evidence_pipeline.py ... -- -with_report`.

## Outputs

| Path (under `${params.outdir}`) | Contents |
| --- | --- |
| `hintsfile.gff` | Combined genome annotation hints from extrinsic evidence. |
| `training.gff` | Filtered high-confidence gene set after stop-codon checks. |
| `tiberius.gff3` | Raw merged Tiberius genome annotation (if `tiberius.run: true`). |
| `tiberius_train.gff3` | Union of Tiberius predictions and HC training genes. |
| `tiberius_train_prio.gff3` | Same as above but prioritising HC genes on conflicting gene loci. |
| `sra_downloads/*.fastq.gz` | FASTQs produced by `DOWNLOAD_SRA_*` processes. |

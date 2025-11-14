process DOWNLOAD_SRA_PAIRED {
  // container params.sra_containers

  publishDir "${params.outdir}/sra_downloads", mode: 'copy'

  input:
    val acc

  output:
    tuple val(acc), path("${acc}_1.fastq.gz"), path("${acc}_2.fastq.gz")

  script:
  """
  fasterq-dump --split-files --threads ${params.threads} ${acc}
  gzip ${acc}_1.fastq ${acc}_2.fastq
  """
}

process DOWNLOAD_SRA_SINGLE {
  // container params.sra_containers

  publishDir "${params.outdir}/sra_downloads", mode: 'copy'

  input:
    val acc

  output:
    tuple val(acc), path("${acc}.fastq.gz")

  script:
  """
  fasterq-dump --threads ${params.threads} ${acc}
  gzip ${acc}.fastq
  """
}

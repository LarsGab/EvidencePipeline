nextflow.enable.dsl=2

process CONCAT_HINTS {
  publishDir "${params.outdir}", mode: 'copy'   

  input:
    path(prot)
    path(rnaseq), stageAs: 'rnaseq_hints.gff'
    path(isoseq), stageAs: 'isoseq_hints.gff'

  output:
    path "hintsfile.gff", emit: hints

  script:
  """
  cat ${prot} ${ rnaseq ?: "/dev/null" } ${ isoseq ?: "/dev/null" } > hintsfile.gff
  """
}

process EMPTY_FILE {
  output:
    path 'empty.txt'

  script:
  """
    touch empty.txt
  """
}

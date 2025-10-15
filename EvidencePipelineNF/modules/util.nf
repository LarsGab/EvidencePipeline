nextflow.enable.dsl=2

process CONCAT_HINTS {
  publishDir "${params.outdir}/hints", mode: 'copy'   

  input:
    path(prot)
    path(rnaseq)
    path(isoseq)

  output:
    path "hints/all_hints.gff", emit: hints

  script:
  """
  mkdir -p hints
  cat ${prot} ${ rnaseq ?: "/dev/null" } ${ isoseq ?: "/dev/null" } > hints/all_hints.gff
  """
}

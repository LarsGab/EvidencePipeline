nextflow.enable.dsl=2

process STRINGTIE_ASSEMBLE_RNA {
  input:
    path genome
    path rnabam
  output:
    path "stringtie/transcripts.gtf", emit: gtf
  script:
  """
  mkdir -p stringtie
  ${params.tools.stringtie} -p ${task.cpus} -o stringtie/transcripts.gtf ${rnabam}
  """
}

process STRINGTIE_ASSEMBLE_ISO {
  input:
    path genome
    path isobam
  output:
    path "stringtie/transcripts.gtf", emit: gtf
  script:
  """
  mkdir -p stringtie
  ${params.tools.stringtie} -p ${task.cpus} -o stringtie/transcripts.gtf -L ${isobam}
  """
}

process STRINGTIE_ASSEMBLE_MIX {
  input:
    path genome
    path rnabam
    path isobam
  output:
    path "stringtie/transcripts.gtf", emit: gtf
  script:
  """
  mkdir -p stringtie
  ${params.tools.stringtie} -p ${task.cpus} -o stringtie/transcripts.gtf --mix ${rnabam} ${isobam}
  """
}

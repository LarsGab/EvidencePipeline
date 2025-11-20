nextflow.enable.dsl=2

process STRINGTIE_ASSEMBLE_RNA {
  label 'container'
  input:
    path genome
    path rnabam
  output:
    path "stringtie/stringtie.gtf",  emit: gtf
    path "stringtie/stringtie.gff3", emit: gff3
  script:
  """
  mkdir -p stringtie
  ${params.tools.stringtie} -p ${params.threads} -o stringtie/stringtie.gtf ${rnabam}
  ${params.tools.transdecoder_gtf2gff} stringtie/stringtie.gtf > stringtie/stringtie.gff3 
  """
}

process STRINGTIE_ASSEMBLE_ISO {
  label 'container'
  input:
    path genome
    path isobam
  output:
    path "stringtie/stringtie.gtf",  emit: gtf
    path "stringtie/stringtie.gff3", emit: gff3
  script:
  """
  mkdir -p stringtie
  ${params.tools.stringtie} -p ${params.threads} -o stringtie/stringtie.gtf -L ${isobam}
  ${params.tools.transdecoder_gtf2gff} stringtie/stringtie.gtf > stringtie/stringtie.gff3 
  """
}

process STRINGTIE_ASSEMBLE_MIX {
  label 'container'
  input:
    path genome
    path rnabam,  stageAs: 'rna.bam'
    path isobam,  stageAs: 'isoseq.bam'
  output:
    path "stringtie/stringtie.gtf",  emit: gtf
    path "stringtie/stringtie.gff3", emit: gff3
  script:
  """
  mkdir -p stringtie
  ${params.tools.stringtie} -p ${params.threads} -o stringtie/stringtie.gtf --mix rna.bam isoseq.bam
  ${params.tools.transdecoder_gtf2gff} stringtie/stringtie.gtf > stringtie/stringtie.gff3 
  """
}

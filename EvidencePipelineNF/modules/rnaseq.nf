nextflow.enable.dsl=2
process HISAT2_BUILD {
  input:
    path genome
  output:
    path "hisat2", emit: idxdir

  script:
  """
  mkdir -p hisat2
  ${params.tools.hisat2_build} -p ${task.cpus} ${genome} hisat2/genome
  """
}

process HISAT2_MAP_SINGLE {
  input:
    path idxdir         
    path reads           
  output:
    path "rnaseq/${reads.simpleName}.sam", emit: sam

  script:
  """
  mkdir -p rnaseq
  ${params.tools.hisat2} -x ${idxdir}/genome -U ${reads} --dta -p ${task.cpus} -S rnaseq/${reads.simpleName}.sam
  """
}

process HISAT2_MAP_PAIRED {
  input:
    path idxdir
    tuple val(sample), path(read1), path(read2)
  output:
    path "rnaseq/${sample}.sam", emit: sam

  script:
  """
  mkdir -p rnaseq
  ${params.tools.hisat2} -x ${idxdir}/genome -1 ${read1} -2 ${read2} --dta -p ${task.cpus} -S rnaseq/${sample}.sam
  """
}




process SAMTOOLS_VIEW_SORT {
  input: path sam
  output: path "${sam.baseName}.bam", emit: bam
  script: """
  ${params.tools.samtools} view -bS ${sam} | ${params.tools.samtools} sort -@ ${params.threads} -o ${sam.baseName}.bam
  """
}

process SAMTOOLS_MERGE {
  input:
    path bams

  output:
    path "rnaseq/merged.bam", emit: bam

  script:
  """
  mkdir -p rnaseq
  ${params.tools.samtools} merge -f -@ ${task.cpus} -o rnaseq/merged.bam ${bams}
  """
}


process BAM2HINTS {
  input: path bam; path genome
  output: path "${bam.simpleName}.hints.gff", emit: hints
  script: """
  ${params.tools.bam2hints} --intronsonly --in=${bam} --out=${bam}.temp
  ${params.tools.filterIntronsFindStrand} ${genome} ${bam}.temp --score > ${bam.simpleName}.hints.gff
  """
}

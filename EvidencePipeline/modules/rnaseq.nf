nextflow.enable.dsl=2
process HISAT2_BUILD {
  label 'container'
  input:
    path genome

  output:
    path "hisat2_idx", emit: idxdir

  script:
  """
  mkdir -p hisat2_idx
  ${params.tools.hisat2_build} -p ${task.cpus} ${genome} hisat2_idx/genome
  """
}

process HISAT2_MAP_SINGLE {
  label 'container'
  input:
    path idxdir         
    path reads           
  output:
    path "${reads.simpleName}.sam", emit: sam

  script:
  """
  ${params.tools.hisat2} -x ${idxdir}/genome -U ${reads} --dta -p ${task.cpus} -S ${reads.simpleName}.sam
  """
}

process HISAT2_MAP_PAIRED {
  label 'container'
  input:
    path idxdir
    tuple val(sample), path(read1), path(read2)
  output:
    path "${sample}.sam", emit: sam

  script:
  """
  ${params.tools.hisat2} -x ${idxdir}/genome -1 ${read1} -2 ${read2} --dta -p ${task.cpus} -S ${sample}.sam
  """
}

process SAMTOOLS_VIEW_SORT {
  label 'container'
  input: path sam
  output: path "${sam.baseName}.bam", emit: bam
  script: """
  ${params.tools.samtools} view -bS ${sam} | ${params.tools.samtools} sort -@ ${params.threads} -o ${sam.baseName}.bam
  """
}

process SAMTOOLS_MERGE {
  label 'container'
  input:
    path bams

  output:
    path "merged.bam", emit: bam

  script:
  """
  ${params.tools.samtools} merge -f -@ ${task.cpus} -o merged.bam ${bams}
  """
}


process BAM2HINTS {
  label 'container'
  input: path bam; path genome
  output: path "${bam.simpleName}.hints.gff", emit: hints
  script: """
  ${params.tools.bam2hints} --intronsonly --in=${bam} --out=${bam}.temp
  ${projectDir}/scripts/filterIntronsFindStrand.pl ${genome} ${bam}.temp --score > ${bam.simpleName}.hints.gff
  """
}

nextflow.enable.dsl=2

process HISAT2_BUILD {
  input: path genome
  output: path "hisat2/genome.*.ht2", emit: idx
  script: """
  mkdir -p hisat2
  ${params.tools.hisat2_build} -p ${params.threads} ${genome} hisat2/genome
  """
}

process HISAT2_MAP_SINGLE {
  input: path idx; path reads
  output: path "rnaseq/${reads.simpleName}.sam", emit: sam
  script: """
  mkdir -p rnaseq
  ${params.tools.hisat2} -x hisat2/genome -U ${reads} --dta -p ${params.threads} -S rnaseq/${reads.simpleName}.sam
  """
}

process HISAT2_MAP_PAIRED {
  input: path idx; path read
  output: path "rnaseq/${pair}.sam", emit: sam
  when: read.name ==~ /.*(_1|R1)\..*/
  script:
  def pair = read.baseName.replaceAll(/(_1|R1)$/, "")
  def mate2 = file(read.name.replaceFirst(/(_1|R1)/, "_2").replaceFirst(/R1/, "R2"))
  """
  mkdir -p rnaseq
  ${params.tools.hisat2} -x hisat2/genome -1 ${read} -2 ${mate2} --dta -p ${params.threads} -S rnaseq/${pair}.sam
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
  input: path(bams)
  output: path "rnaseq/merged.bam", emit: bam
  script: """
  mkdir -p rnaseq
  ${params.tools.samtools} merge -f -o rnaseq/merged.bam ${" ".join(bams.collect{it.getName()})}
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

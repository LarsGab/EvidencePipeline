nextflow.enable.dsl=2

process MINIPROT_ALIGN {
  input: path genome; path proteins
  output: 
    path "miniprot/miniprot.aln", emit: aln
  script: """
  mkdir -p miniprot
  ${params.tools.miniprot} -t ${params.threads} --aln ${genome} ${proteins} > miniprot/miniprot.aln
  """
}

process MINIPROT_BOUNDARY_SCORE {
  input:
    path aln
    path score_matrix

  output:
    path "miniprot/miniprot_parsed.gff", emit: gff

  script:
  """
  mkdir -p miniprot
  ${params.tools.miniprot_boundary_scorer} \
    -s ${score_matrix} \
    -o miniprot/miniprot_parsed.gff \
    < ${aln}
  """
}

process MINIPROTHINT_CONVERT {
  input: path gff
  output: 
    path "miniprot/miniprot.gtf", emit: gtf
    path "miniprot/miniprot_trainingGenes.gff", emit: traingff
  script: """
  mkdir -p miniprot
  ${params.tools.miniprothint} ${gff} --workdir miniprot --ignoreCoverage --topNperSeed 10 --minScoreFraction 0.5
  """
}

process ALN2HINTS {
  input: path gtf
  output: path "hints/hints_protein.gff", emit: hints
  script: """
  mkdir -p hints
  ${projectDir}/scripts/aln2hints.pl --in=${gtf} --out=hints/prot_hintsfile.aln2hints.temp.gff --prg=miniprot --priority=4
  cp hints/prot_hintsfile.aln2hints.temp.gff hints/hints_protein.gff
  """
}

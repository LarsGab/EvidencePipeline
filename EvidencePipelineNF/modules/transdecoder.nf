nextflow.enable.dsl=2

process TD_GTF_TO_FA {
  input: path gtf; path genome
  output: path "transdecoder/transcripts.fasta", emit: cdna
  script: """
  mkdir -p transdecoder
  ${params.tools.transdecoder_util_gtf2fa} ${gtf} ${genome} > transdecoder/transcripts.fasta
  """
}
process TD_LONGORFS {
  input: path cdna
  output: path "transdecoder/longorfs", emit: longdir
  script: """
  mkdir -p transdecoder/longorfs
  ${params.tools.transdecoder_longorfs} -O transdecoder -t ${cdna}
  """
}
process TD_PREDICT {
  input: path cdna
  output: path "transdecoder/transcripts.fasta.transdecoder.pep", emit: pep
  script: """
  ${params.tools.transdecoder_predict} -O transdecoder -t ${cdna}
  """
}
process PY_SHORTEN_INCOMPLETE_ORFS {
  input: path pep
  output: path "transdecoder/shortened_candidates.pep", emit: pep_short
  script: """
  python - << 'PY'
from EvidencePipeline.prep_evidence import shorten_incomplete_orfs
shorten_incomplete_orfs("${pep}", "transdecoder/shortened_candidates.pep")
PY
  """
}
process CDS_CLASSIFY_AND_REVISE {
  input: path diamond_normal; path diamond_short; path transdecoder_pep; path shortened_pep
  output: path "transdecoder/revised_candidates.pep", emit: revised_pep
  output: path "transdecoder/classifications.json", emit: classes
  script: """
  python wrappers/classify_and_revise.py     --diamond_normal ${diamond_normal}     --diamond_short ${diamond_short}     --transdecoder_pep ${transdecoder_pep}     --shortened_pep ${shortened_pep}     --revised_pep transdecoder/revised_candidates.pep     --classifications_json transdecoder/classifications.json
  """
}

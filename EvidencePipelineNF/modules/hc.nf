nextflow.enable.dsl=2

process HC_SUPPORTED_BY_PROTEINS {
  input:
    path diamond_revised
    path revised_pep
    path proteins

  output:
    path "hc/hc_genes.protein.pep", emit: hc_prot_pep
    path "hc/q_dict.json",         emit: qdict

  script:
  """
  mkdir -p hc
  python wrappers/hc_by_proteins.py \
    --diamond_revised ${diamond_revised} \
    --revised_pep     ${revised_pep} \
    --proteins        ${proteins} \
    --hc_pep_out      hc/hc_genes.protein.pep \
    --qdict_json      hc/q_dict.json
  """
}

process HC_SUPPORTED_BY_INTRINSIC_AND_TRAINING {
  publishDir "${params.outdir}", pattern: "hc/training.gff", mode: 'copy'
  
  input:
    path qdict
    path stringtie_gtf
    path stringtie_gff3
    path transcripts_fasta
    path miniprot_gff

  output:
    path "hc/training.gff", emit: training_gff
    path "hc/hc_genes.pep", emit: hc_pep
    path "hc/lc_genes.pep", emit: lc_pep

  script:
  """
  mkdir -p hc
  python wrappers/hc_by_intrinsic_and_training.py \
    --qdict_json        ${qdict} \
    --stringtie_gtf     ${stringtie_gtf} \
    --stringtie_gff3    ${stringtie_gff3} \
    --transcripts_fasta ${transcripts_fasta} \
    --miniprot_gff      ${miniprot_gff} \
    --transdecoder_util ${params.tools.transdecoder_util_orf2genome} \
    --bedtools_path     ${params.tools.bedtools} \
    --outdir            hc \
    --training_gff      hc/training.gff
  """
}

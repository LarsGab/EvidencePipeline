nextflow.enable.dsl=2

process MINIPROT_ALIGN {
  label 'container', 'bigmem'

  input: path genome; path proteins

  output: 
    path "miniprot/miniprot.aln", emit: aln

  script: """
  mkdir -p miniprot
  ${params.tools.miniprot} -t ${params.threads} --aln ${genome} ${proteins} > miniprot/miniprot.aln
  """
}

process MINIPROT_BOUNDARY_SCORE {
  label 'container'

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
  label 'container'

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

  output: path "hints_protein.gff", emit: hints
  
  script: """
  ${projectDir}/scripts/aln2hints.pl --in=${gtf} --out=prot_hintsfile.aln2hints.temp.gff --prg=miniprot --priority=4
  cp prot_hintsfile.aln2hints.temp.gff hints_protein.gff
  """
}

process PREPROCESS_PROTEINDB {
  label 'container', 'bigmem'

  input: 
    path proteinDB 
    path tiberius_prot

  output: path "protein_preprocessed.fa"

  script: """
    # Count protein sequences (FASTA headers start with '>')
    N_PROT=\$(grep -c '^>' "${proteinDB}" || echo 0)

    echo "[PREPROCESS_PROTEINDB] Number of proteins in input: \$N_PROT" >&2

    if [[ "\$N_PROT" -le 1000000 ]]; then
        echo "[PREPROCESS_PROTEINDB] <= 1,000,000 proteins – using full DB." >&2
        ln -s "${proteinDB}" protein_preprocessed.fa
    else
        echo "[PREPROCESS_PROTEINDB] > 1,000,000 proteins – running DIAMOND soft filter." >&2

        diamond makedb \
          --in "${tiberius_prot}" \
          --db tib_db

        diamond blastp \
          --query "${proteinDB}" \
          --db tib_db \
          --out diamond_hits.tsv \
          --outfmt 6 qseqid sseqid pident length evalue bitscore \
          --evalue 1e-3 \
          --id 20 \
          --min-orf 40 \
          --max-target-seqs 5 \
          --threads ${params.threads} \
          --very-sensitive

        # Collect unique query IDs (proteins from proteinDB that hit something)
        cut -f1 diamond_hits.tsv | sort -u > keep_ids.txt

        N_KEEP=\$(wc -l < keep_ids.txt || echo 0)
        echo "[PREPROCESS_PROTEINDB] Proteins kept after DIAMOND filter: \$N_KEEP" >&2

        # Extract only those proteins
        seqkit grep -f keep_ids.txt "${proteinDB}" > protein_preprocessed.fa
    fi
  """
}
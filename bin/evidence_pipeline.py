import logging, os, time
import evi_initializer as p_i
import utility as util
import prep_evidence as p_e
import generate_hc_genes as hc
import generate_hints as hints


def main():
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s')
    args = p_i.get_parser()
    start_time = time.time()

    util.make_directory(args.output_path)    
    log_file = args.output_path + '/LOG.out' if args.output_path else 'LOG.out'
    p_i.setup_logger(log_file)
    logging.info('Running EvidencePipeline for generating evidence for genome annotation.')
    logging.info(f'The log messages will be saved in {log_file} .')
    logging.info
    
    cfg = p_i.get_config(args) 
    cfg["mode"] = p_i.determine_mode(cfg)
    cfg = p_i.check_dependencies(cfg)

    cfg["output_path"] = util.make_absolute_path(cfg["output_path"])
    os.chdir(cfg["output_path"])

    # init empty hintsfile
    hints_file = os.path.join(cfg["output_path"], "hintsfile.gff")    
    with open(hints_file, 'w+') as h_f:
        pass

    prot_aln_gff, prot_aln_gtf = hc.align_proteins(cfg["genome"], cfg["proteins"], 
            cfg["MINIPROT"], cfg["MINIPROT_BOUNDARY_SCORER"], cfg["MINIPROT_HINT"],
            alignment_scoring=cfg["scoring_matrix"], threads=cfg["threads"], 
            output_dir=f"{cfg['output_path']}/miniprot/", use_existing=cfg["restart"])

    hints.aln2hints(prot_aln_gtf, hints_file, cfg["GALBA_SCRIPTS"])

    alignment_rnaseq = None
    alignment_isoseq = None
    # align rnaseq
    if cfg["mode"] in ['mixed', 'rnaseq']:
        p_e.indexing(cfg["genome"], cfg["HISAT"], cfg["threads"])
        alignments_list = p_e.mapping_short(cfg["rnaseq_paired"], cfg["rnaseq_single"], 
                cfg["HISAT"], cfg["threads"], f'{cfg["output_path"]}/rnaseq_mapping',
                use_existing=cfg["restart"])
        if len(alignments_list) > 1:
            alignment_rnaseq = util.sam_to_bam(alignments_list, 
                    merge_bam=f'{cfg["output_path"]}/rnaseq_mapping/merged.bam', use_existing=cfg["restart"])
        else:
            alignment_rnaseq = util.sam_to_bam(alignments_list, samtools_path=cfg["SAMTOOLS"],
                        use_existing=cfg["restart"])[0]
        # generate hints from bam
        rnaseq_hints = hints.bam2hints(alignment_rnaseq, cfg["genome"],
            f"{cfg['output_path']}/rnaseq_mapping/hintsfile_rnaseq.gff", hints_file,
            cfg["BAM2HINTS"], cfg["GALBA_SCRIPTS"])

    # align isoseq
    if cfg["mode"] in ['mixed', 'isoseq']:
        minimap_isoseq = p_e.mapping_long(cfg["genome"], cfg["isoseq"], cfg["MINIMAP"], 
                threads=cfg["threads"], output_dir=f"{cfg['output_path']}/isoseq_mapping/",
                use_existing=cfg["restart"])
        alignment_isoseq = util.sam_to_bam([minimap_isoseq], samtools_path=cfg["SAMTOOLS"],
                use_existing=cfg["restart"])[0]
        # generate hints from bam
        isoseq_hints = hints.bam2hints(alignment_isoseq, cfg["genome"],
            f"{cfg['output_path']}/isoseq_mapping/hintsfile_isoseq.gff", hints_file,
            cfg["BAM2HINTS"], cfg["GALBA_SCRIPTS"])

    # assemble reads
    if cfg["mode"] in ['mixed', 'rnaseq', 'isoseq']:
        stringtie_gtf = p_e.assembling(alignment_rnaseq, alignment_isoseq, cfg['STRINGTIE'], 
                threads=cfg["threads"], output_dir=f"{cfg['output_path']}/stringtie/", 
                mode=cfg['mode'], use_existing=cfg["restart"])
        stringtie_gff3 = util.convert_gtf_to_gff3(stringtie_gtf, cfg['TRANSDECODER_UTIL'],
                output_gff3=f"{cfg['output_path']}/stringtie/transcripts.gff3")
        transdecoder_pep, transdecoder_fasta = p_e.orfsearching(cfg["genome"], 
                stringtie_gtf, cfg["TRANSDECODER"], 
                output_dir=f"{cfg['output_path']}/transdecoder/", use_existing=cfg["restart"],
                transdecoder_util_path=cfg['TRANSDECODER_UTIL'])

        shortened_pep = p_e.shorten_incomplete_orfs(transdecoder_pep, 
                output_path=f"{cfg['output_path']}/transdecoder/shortened_candidates.pep")
        protein_db_path = p_e.make_diamond_db(cfg["proteins"], cfg["DIAMOND"], 
                output_dir=f"{cfg['output_path']}/diamond/", db_name="protein_db")
        diamond_shortened = p_e.validating_orfs(shortened_pep, protein_db_path, cfg["DIAMOND"], 
                output_path=f"{cfg['output_path']}/diamond/diamond_shortened.tsv", use_existing=cfg["restart"])

        protein_db_path = p_e.make_diamond_db(cfg["proteins"], cfg["DIAMOND"], 
                output_dir=f"{cfg['output_path']}/diamond/", db_name="protein_db")
        diamond_normal = p_e.validating_orfs(transdecoder_pep, protein_db_path, cfg["DIAMOND"],
                output_path=f"{cfg['output_path']}/diamond/diamond_normal.tsv", use_existing=cfg["restart"])

        classifications_dict = p_e.get_cds_classification(diamond_shortened, diamond_normal)
        revised_pep = p_e.get_optimized_pep_file(transdecoder_pep, shortened_pep, 
                classifications_dict, output_path=f"{cfg['output_path']}/transdecoder/revised_candidates.pep")
        protein_db_path = p_e.make_diamond_db(cfg["proteins"], cfg["DIAMOND"], 
                output_dir=f"{cfg['output_path']}/diamond/", db_name="protein_db")
        diamond_revised = p_e.validating_orfs(revised_pep, protein_db_path, cfg["DIAMOND"],
                output_path=f"{cfg['output_path']}/diamond/diamond_revised.tsv", use_existing=cfg["restart"])

        # generate hc genes
        hc_out_dir = f"{cfg['output_path']}/hc_genes/"
        util.make_directory(hc_out_dir)

        q_dict = hc.getting_hc_supported_by_proteins(diamond_revised, revised_pep, cfg["proteins"], 
                output_path=f"{hc_out_dir}/hc_genes.pep")       
        hc_genes, lc_genes = hc.getting_hc_supported_by_intrinsic(q_dict, stringtie_gtf, stringtie_gff3, 
                    transdecoder_fasta, prot_aln_gff, cfg['TRANSDECODER_UTIL'], cfg['BEDTOOLS'],
                output_dir=hc_out_dir, min_cds_len=300, score_threshold=50)
        hc_single_pep = hc.choose_one_isoform(hc_genes, f"{hc_out_dir}/hc_one_isoform.pep")
        hc_single_gff = hc.from_pep_file_to_gff3(hc_single_pep, stringtie_gtf, 
                f"{hc_out_dir}/hc_one_isoform.gff3")
        hc.from_transcript_to_genome(hc_single_gff, stringtie_gff3, transdecoder_fasta, 
                f"{output_path}/training.gff", cfg['TRANSDECODER_UTIL'])
    else:
        logging.info("Using Miniprot training genes as training genes.")
        util.copy_file(f"{cfg['output_path']}/miniprot/miniprot_trainingGenes.gff", f"{output_path}/training.gff")


    
    logging.info(f"------> Training genes file {output_path}/training.gff3 was generated successfully!")

    duration = time.time()  - start_time
    logging.info(f"Finished! EvidencePipeline took {duration/60:.4f} minutes to execute.")

if __name__ == "__main__":
    main()
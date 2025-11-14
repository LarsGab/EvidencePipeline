nextflow.enable.dsl=2

include { MINIPROT_ALIGN; MINIPROT_BOUNDARY_SCORE; MINIPROTHINT_CONVERT; ALN2HINTS } from './modules/proteins.nf'
include {
  HISAT2_BUILD
  HISAT2_MAP_SINGLE
  HISAT2_MAP_PAIRED
  SAMTOOLS_VIEW_SORT as SAMTOOLS_VIEW_SORT_SINGLE
  SAMTOOLS_VIEW_SORT as SAMTOOLS_VIEW_SORT_PAIRED
  SAMTOOLS_VIEW_SORT as SAMTOOLS_VIEW_SORT_ISO
  SAMTOOLS_MERGE as SAMTOOLS_MERGE_RNA
  SAMTOOLS_MERGE as SAMTOOLS_MERGE_ISO
  BAM2HINTS as BAM2HINTS_RNA
} from './modules/rnaseq.nf'

include {
  MINIMAP2_MAP
} from './modules/isoseq.nf'

include {
  BAM2HINTS                 as BAM2HINTS_ISO
} from './modules/rnaseq.nf'

include {
  DIAMOND_MAKEDB
  DIAMOND_BLASTP as DIAMOND_BLASTP_NORM
  DIAMOND_BLASTP as DIAMOND_BLASTP_SHORT
  DIAMOND_BLASTP as DIAMOND_BLASTP_REV
} from './modules/diamond.nf'
include {
    DOWNLOAD_SRA_PAIRED
    DOWNLOAD_SRA_SINGLE as DOWNLOAD_SRA_SINGLE
    DOWNLOAD_SRA_SINGLE as DOWNLOAD_SRA_ISO
} from './modules/download.nf'
include { RUN_TIBERIUS; SPLIT_GENOME; MERGE_TIBERIUS; MERGE_TIBERIUS_TRAIN_PRIO; MERGE_TIBERIUS_TRAIN} from './modules/tiberius.nf'
include { CONCAT_HINTS; EMPTY_FILE } from './modules/util.nf'
include { STRINGTIE_ASSEMBLE_RNA; STRINGTIE_ASSEMBLE_ISO; STRINGTIE_ASSEMBLE_MIX } from './modules/assembly.nf'
include { TD_ALL; SHORTEN_INCOMPLETE_ORFS; CDS_CLASSIFY_AND_REVISE;} from './modules/transdecoder.nf'
include { HC_SUPPORTED; HC_FORMAT_FILTER } from './modules/hc.nf'
include { REMOVE_CONFLICTING_WITH_PROTHINT } from './modules/conflicts.nf'

log.info "DEBUG params: genome=${params.genome} proteins=${params.proteins} outdir=${params.outdir}"


def inferMode(boolean hasPaired, boolean hasSingle, boolean hasIso, boolean hasProteins) {
    if (!hasProteins) error 'params.proteins is required'
    if ((hasPaired || hasSingle) && hasIso) return 'mixed'
    if (hasIso) return 'isoseq'
    if (hasPaired || hasSingle) return 'rnaseq'
    return 'proteins'
}

workflow {
    main:
    if( !params.genome )   error "params.genome is required"
    if( !params.proteins ) error "params.proteins is required"

    def genomeFile   = file(params.genome)
    def proteinsFile = file(params.proteins)

    if( !genomeFile.exists() )
        error "Genome file not found on disk: ${genomeFile.toString()}"

    if( !proteinsFile.exists() )
        error "Proteins file not found on disk: ${proteinsFile.toString()}"

    CH_GENOME   = Channel.value(genomeFile)
    CH_PROTEINS = Channel.value(proteinsFile)

    def DO_SE_LOCAL = params.rnaseq_single && params.rnaseq_single.size() > 0
    def DO_PE_LOCAL = params.rnaseq_paired && params.rnaseq_paired.size() > 0
    def DO_ISO_LOCAL = params.isoseq && params.isoseq.size() > 0

    def DO_SE = DO_SE_LOCAL ||
        (params.rnaseq_sra_single && params.rnaseq_sra_single.size() > 0)
    def DO_PE = DO_PE_LOCAL ||
        (params.rnaseq_sra_paired && params.rnaseq_sra_paired.size() > 0)
    def DO_ISO = DO_ISO_LOCAL ||
        (params.isoseq_sra && params.isoseq_sra.size() > 0)    

    

    if( !params.scoring_matrix )
        error "params.scoring_matrix is required by miniprot_boundary_scorer"

    def SCORE = file(params.scoring_matrix)
    if( !SCORE.exists() )
        error "Score matrix not found: ${SCORE.toString()}"

    CH_SCORE = Channel.value(SCORE)

    def hasPaired   = params.rnaseq_paired && params.rnaseq_paired.size() > 0
    def hasSingle   = params.rnaseq_single && params.rnaseq_single.size() > 0
    def hasIso      = params.isoseq && params.isoseq.size() > 0
    def hasProteins = params.proteins != null

    def MODE = params.mode ?: inferMode(hasPaired, hasSingle, hasIso, hasProteins)
    log.info "Running mode: ${MODE}"
    outdir = params.outdir
    file(outdir).mkdirs()

    empty_file = EMPTY_FILE()

    // RNA-Seq PE channel 
    if (DO_PE) {
        CH_PAIRED_LOCAL = DO_PE_LOCAL ?
            Channel.fromFilePairs(params.rnaseq_paired, flat: true, checkIfExists: true)
            : Channel.empty()
        if (params.rnaseq_sra_paired) {
            CH_RNASEQ_SRA_IDS_PAIRED = params.rnaseq_sra_paired ?
                Channel.from(params.rnaseq_sra_paired)
                : Channel.empty()

            CH_RNASEQ_PAIRED_SRA = DOWNLOAD_SRA_PAIRED(CH_RNASEQ_SRA_IDS_PAIRED)

            CH_PAIRED_SRA = CH_RNASEQ_PAIRED_SRA.map { acc, r1, r2 ->
                tuple(acc, [r1, r2])
            }

            CH_PAIRED = CH_PAIRED_LOCAL.mix(CH_PAIRED_SRA)
        } else {
            CH_PAIRED = CH_PAIRED_LOCAL
        }
    }

    // RNA-Seq SE channel
    if (DO_SE) {
        CH_SINGLE_LOCAL = DO_SE_LOCAL
            ? Channel.fromPath(params.rnaseq_single, checkIfExists: true)
            : Channel.empty()

        if (params.rnaseq_sra_single) {
            CH_RNASEQ_SRA_IDS_SINGLE = params.rnaseq_sra_single ?
                Channel.from(params.rnaseq_sra_single)
                : Channel.empty()

            CH_RNASEQ_SINGLE_SRA = DOWNLOAD_SRA_SINGLE(CH_RNASEQ_SRA_IDS_SINGLE)

            CH_SINGLE = CH_SINGLE_LOCAL.mix(CH_RNASEQ_SINGLE_SRA.map { acc, f -> f })
        } else {
            CH_SINGLE = CH_SINGLE_LOCAL
        }
    }

    // Isoseq channel
    if (DO_ISO) {
        CH_ISO_LOCAL = DO_ISO_LOCAL
            ? Channel.fromPath(params.isoseq, checkIfExists: true)
            : Channel.empty()
        
        if (params.isoseq_sra) {
            CH_ISO_SRA_IDS = params.isoseq_sra ?
                Channel.from(params.isoseq_sra)
                : Channel.empty()

            CH_RNASEQ_ISO_SRA = DOWNLOAD_SRA_ISO(CH_ISO_SRA_IDS)

            CH_ISO = CH_ISO_LOCAL.mix(CH_RNASEQ_ISO_SRA.map { acc, f -> f })
        } else {
            CH_ISO = CH_ISO_LOCAL
        }
    }
  


    // Tiberius
    if (params.tiberius.run){
        genome_split = SPLIT_GENOME(CH_GENOME)
        chunks_ch      = genome_split.chunks.flatten()
        tiberius_split = RUN_TIBERIUS(chunks_ch, params.tiberius.model_cfg)
        tiberius = MERGE_TIBERIUS(tiberius_split.toList())
    }

    // Protein
    prot_aln = MINIPROT_ALIGN(CH_GENOME, CH_PROTEINS)
    scored   = MINIPROT_BOUNDARY_SCORE(prot_aln.aln, CH_SCORE)
    prot_gtf = MINIPROTHINT_CONVERT(scored.gff)
    prot_hints = ALN2HINTS(prot_gtf.gtf)
    //miniprot_gff_for_conf = file("miniprot/miniprot_parsed.gff")

    // RNA-seq
    if( MODE in ['mixed','rnaseq'] && (params.rnaseq_paired.size()>0 || params.rnaseq_single.size()>0) ){
        index        = HISAT2_BUILD(CH_GENOME)
        rnaseq_bams  = Channel.empty()

        if (DO_SE) {
            map_se     = HISAT2_MAP_SINGLE(index.idxdir, CH_SINGLE)
            sort_se    = SAMTOOLS_VIEW_SORT_SINGLE(map_se.sam)
            rnaseq_bams = rnaseq_bams.mix(sort_se.bam)
        }

        if (DO_PE) {
            map_pe     = HISAT2_MAP_PAIRED(index.idxdir, CH_PAIRED)
            sort_pe    = SAMTOOLS_VIEW_SORT_PAIRED(map_pe.sam)
            rnaseq_bams = rnaseq_bams.mix(sort_pe.bam)
        }

        rnaseq_merged = SAMTOOLS_MERGE_RNA(rnaseq_bams.collect())
        rnaseq_hints  = BAM2HINTS_RNA(rnaseq_merged.bam, CH_GENOME)
        HINTS_RNASEQ  = rnaseq_hints.hints
        BAM_FOR_ASM   = rnaseq_merged.bam
    } else {
        HINTS_RNASEQ = empty_file
        BAM_FOR_ASM  = empty_file
    }

    //  Iso-Seq
    if( MODE in ['mixed','isoseq'] && DO_ISO ){
        iso_bams = Channel.empty()
        iso_sam   = MINIMAP2_MAP(CH_GENOME, CH_ISO)
        iso_sort   = SAMTOOLS_VIEW_SORT_ISO(iso_sam.sam)
        iso_bams = iso_bams.mix(iso_sort.bam)
        iso_merged = SAMTOOLS_MERGE_ISO(iso_bams.collect())
        iso_hints  = BAM2HINTS_ISO(iso_merged.bam, CH_GENOME)
        HINTS_ISO    = iso_hints.hints
        ISO_BAM      = iso_merged.bam
    } else {
        HINTS_ISO = empty_file
        ISO_BAM   = empty_file
    }

    // Hints concat
    all_hints = CONCAT_HINTS(prot_hints.hints, HINTS_RNASEQ, HINTS_ISO )

    // Assembly + ORFs + DIAMOND + HC + Training
    if( MODE in ['mixed','rnaseq','isoseq'] ){
        def asm

        if( MODE == 'mixed' ) {
        asm = STRINGTIE_ASSEMBLE_MIX( CH_GENOME, BAM_FOR_ASM, ISO_BAM )
        }
        else if( MODE == 'rnaseq' ) {
        asm = STRINGTIE_ASSEMBLE_RNA( CH_GENOME, BAM_FOR_ASM )
        }
        else if( MODE == 'isoseq' ) {
        asm = STRINGTIE_ASSEMBLE_ISO( CH_GENOME, ISO_BAM )
        }
        td_all   = TD_ALL( asm.gtf, CH_GENOME )

        pep_short = SHORTEN_INCOMPLETE_ORFS( td_all.pep )
        db = DIAMOND_MAKEDB(CH_PROTEINS)
        dia_norm = DIAMOND_BLASTP_NORM( td_all.pep, db.db )
        dia_short = DIAMOND_BLASTP_SHORT( pep_short.pep_short,  db.db )
        rev = CDS_CLASSIFY_AND_REVISE( dia_norm.tsv, dia_short.tsv, td_all.pep, pep_short.pep_short )
        dia_rev = DIAMOND_BLASTP_REV( rev.revised_pep, db.db )

        hc = HC_SUPPORTED(
            dia_rev.tsv,
            rev.revised_pep, 
            CH_PROTEINS,
            asm.gtf,
            asm.gff3,
            td_all.cdna,
            scored.gff
            )
        traingff = hc.training_gff
    } else {
        traingff = prot_gtf.traingff

    }
    train_final = HC_FORMAT_FILTER(traingff, params.genome)

    if (params.tiberius.run) {
        MERGE_TIBERIUS_TRAIN(tiberius,train_final)
        MERGE_TIBERIUS_TRAIN_PRIO(tiberius,train_final)
    }
    emit: all_hints.hints
}

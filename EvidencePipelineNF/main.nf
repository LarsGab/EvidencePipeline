nextflow.enable.dsl=2

include { MINIPROT_ALIGN; MINIPROT_BOUNDARY_SCORE; MINIPROTHINT_CONVERT; ALN2HINTS } from './modules/proteins.nf'
include {
  HISAT2_BUILD
  HISAT2_MAP_SINGLE
  HISAT2_MAP_PAIRED
  SAMTOOLS_VIEW_SORT as SAMTOOLS_VIEW_SORT_SINGLE
  SAMTOOLS_VIEW_SORT as SAMTOOLS_VIEW_SORT_PAIRED
  SAMTOOLS_VIEW_SORT as SAMTOOLS_VIEW_SORT_ISO
  SAMTOOLS_MERGE
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
include { CONCAT_HINTS } from './modules/util.nf'
include { STRINGTIE_ASSEMBLE_RNA; STRINGTIE_ASSEMBLE_ISO; STRINGTIE_ASSEMBLE_MIX } from './modules/assembly.nf'
include { TD_GTF_TO_FA; TD_LONGORFS; TD_PREDICT; PY_SHORTEN_INCOMPLETE_ORFS; CDS_CLASSIFY_AND_REVISE } from './modules/transdecoder.nf'
include { HC_SUPPORTED_BY_PROTEINS; HC_SUPPORTED_BY_INTRINSIC_AND_TRAINING } from './modules/hc.nf'
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

    def genomeFile   = file(params.genome)      // does NOT do globbing; just wraps the path
    def proteinsFile = file(params.proteins)

    if( !genomeFile.exists() )
        error "Genome file not found on disk: ${genomeFile.toString()}"

    if( !proteinsFile.exists() )
        error "Proteins file not found on disk: ${proteinsFile.toString()}"

    // Create single-value channels
    CH_GENOME   = Channel.value(genomeFile)
    CH_PROTEINS = Channel.value(proteinsFile)

    CH_PAIRED = Channel.fromList( params.rnaseq_paired )
    CH_SINGLE = Channel.fromList( params.rnaseq_single )
    CH_ISO    = Channel.fromList( params.isoseq )

    def hasPaired   = params.rnaseq_paired && params.rnaseq_paired.size() > 0
    def hasSingle   = params.rnaseq_single && params.rnaseq_single.size() > 0
    def hasIso      = params.isoseq        && params.isoseq.size()        > 0
    def hasProteins = params.proteins != null

    def MODE = params.mode ?: inferMode(hasPaired, hasSingle, hasIso, hasProteins)
    log.info "Running mode: ${MODE}"
    outdir = params.outdir
    file(outdir).mkdirs()

    // Protein branch
    prot_aln = MINIPROT_ALIGN(CH_GENOME, CH_PROTEINS)
    scored   = MINIPROT_BOUNDARY_SCORE(prot_aln.aln)
    prot_gtf = MINIPROTHINT_CONVERT(scored.gff)
    prot_hints = ALN2HINTS(prot_gtf.gtf)
    miniprot_gff_for_conf = file("miniprot/miniprot_parsed.gff")

    // RNA-seq
    if( MODE in ['mixed','rnaseq'] && (params.rnaseq_paired.size()>0 || params.rnaseq_single.size()>0) ){
        index        = HISAT2_BUILD(CH_GENOME)
        rnaseq_bams  = Channel.empty()

        if (params.rnaseq_single.size()>0){
            map_single   = HISAT2_MAP_SINGLE(index.idx, CH_SINGLE)
            sort_single  = SAMTOOLS_VIEW_SORT_SINGLE(map_single.sam)
            rnaseq_bams  = rnaseq_bams.mix( sort_single.bam )
        }
        if (params.rnaseq_paired.size()>0){
            map_paired   = HISAT2_MAP_PAIRED(index.idx, CH_PAIRED)
            sort_paired  = SAMTOOLS_VIEW_SORT_PAIRED(map_paired.sam)
            rnaseq_bams  = rnaseq_bams.mix( sort_paired.bam )
        }

        rnaseq_merged = SAMTOOLS_MERGE(rnaseq_bams.collect())
        rnaseq_hints  = BAM2HINTS_RNA(rnaseq_merged.bam, CH_GENOME)
        HINTS_RNASEQ  = rnaseq_hints.hints
        BAM_FOR_ASM   = rnaseq_merged.bam
    } else {
        HINTS_RNASEQ = Channel.empty()
        BAM_FOR_ASM  = Channel.empty()
    }

    // Iso-Seq
    if( MODE in ['mixed','isoseq'] && params.isoseq.size()>0 ){
        iso_sam   = MINIMAP2_MAP(CH_GENOME, CH_ISO)
        iso_bam   = SAMTOOLS_VIEW_SORT_ISO(iso_sam.sam)   
        isoseq_hints = BAM2HINTS_ISO( iso_bam.bam, CH_GENOME )
        HINTS_ISO    = isoseq_hints.hints
        ISO_BAM      = iso_bam.bam
    } else {
        HINTS_ISO = Channel.empty()
        ISO_BAM   = Channel.empty()
    }

    // Hints concat
    all_hints = CONCAT_HINTS( prot_hints.hints, HINTS_RNASEQ, HINTS_ISO )

    // Assembly + ORFs + DIAMOND + HC + Training
    if( MODE in ['mixed','rnaseq','isoseq'] ){
        def mix = ""
        if (MODE in ['mixed','rnaseq'] && BAM_FOR_ASM)            mix = BAM_FOR_ASM.toString()
        if (MODE in ['mixed','isoseq'] && ISO_BAM)                mix = "-L ${ISO_BAM}"
        if (MODE in ['mixed']        && BAM_FOR_ASM && ISO_BAM)   mix = "--mix ${BAM_FOR_ASM} ${ISO_BAM}"

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
        td_fa = TD_GTF_TO_FA( asm.gtf, CH_GENOME )
        td_long = TD_LONGORFS( td_fa.cdna )
        td_pred = TD_PREDICT( td_fa.cdna )
        pep_short = PY_SHORTEN_INCOMPLETE_ORFS( td_pred.pep )
        db        = DIAMOND_MAKEDB(CH_PROTEINS)
        dia_norm  = DIAMOND_BLASTP_NORM( td_pred.pep,           db.db )
        dia_short = DIAMOND_BLASTP_SHORT( pep_short.pep_short,  db.db )
        rev       = CDS_CLASSIFY_AND_REVISE( dia_norm.tsv, dia_short.tsv, td_pred.pep, pep_short.pep_short )
        dia_rev   = DIAMOND_BLASTP_REV( rev.revised_pep, db.db )

        hc_prot = HC_SUPPORTED_BY_PROTEINS( dia_rev.tsv, rev.revised_pep, CH_PROTEINS )
        hc_intr = HC_SUPPORTED_BY_INTRINSIC_AND_TRAINING( hc_prot.qdict, asm.gtf, file("stringtie/transcripts.gff3"), td_fa.cdna, miniprot_gff_for_conf )
    }

    emit: all_hints.hints
}

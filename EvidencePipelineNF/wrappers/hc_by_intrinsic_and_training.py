#!/usr/bin/env python3
import json, argparse
from Bio.Seq import Seq
from EvidencePipeline import generate_hc_genes as hc
ap = argparse.ArgumentParser()
ap.add_argument("--qdict_json", required=True)
ap.add_argument("--stringtie_gtf", required=True)
ap.add_argument("--stringtie_gff3", required=True)
ap.add_argument("--transcripts_fasta", required=True)
ap.add_argument("--miniprot_gff", required=True)
ap.add_argument("--transdecoder_util", required=True)
ap.add_argument("--bedtools_path", required=True)
ap.add_argument("--outdir", required=True)
ap.add_argument("--training_gff", required=True)
ap.add_argument("--min_cds_len", type=int, default=300)
ap.add_argument("--score_threshold", type=float, default=50.0)
args = ap.parse_args()
with open(args.qdict_json) as f:
    q_serial = json.load(f)
q_dict = {k: (v[0], Seq(v[1])) for k, v in q_serial.items()}
hc_genes, lc_genes = hc.getting_hc_supported_by_intrinsic(
    q_dict, args.stringtie_gtf, args.stringtie_gff3, args.transcripts_fasta,
    args.miniprot_gff, args.transdecoder_util, args.bedtools_path,
    output_dir=args.outdir, min_cds_len=args.min_cds_len, score_threshold=args.score_threshold
)
hc_single_pep = hc.choose_one_isoform(hc_genes, f"{args.outdir}/hc_one_isoform.pep")
hc_single_gff = hc.from_pep_file_to_gff3(hc_single_pep, args.stringtie_gtf, f"{args.outdir}/hc_one_isoform.gff3")
hc.from_transcript_to_genome(hc_single_gff, args.stringtie_gff3, args.transcripts_fasta, args.training_gff, args.transdecoder_util)

#!/usr/bin/env python3
import json, argparse
from EvidencePipeline import prep_evidence as pe
ap = argparse.ArgumentParser()
ap.add_argument("--diamond_normal", required=True)
ap.add_argument("--diamond_short", required=True)
ap.add_argument("--transdecoder_pep", required=True)
ap.add_argument("--shortened_pep", required=True)
ap.add_argument("--revised_pep", required=True)
ap.add_argument("--classifications_json", required=True)
args = ap.parse_args()
classes = pe.get_cds_classification(args.diamond_normal, args.diamond_short)
with open(args.classifications_json, "w") as f:
    json.dump(classes, f)
pe.get_optimized_pep_file(args.transdecoder_pep, args.shortened_pep, classes, args.revised_pep)

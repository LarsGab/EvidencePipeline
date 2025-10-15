#!/usr/bin/env python3
import json, argparse
from EvidencePipeline import generate_hc_genes as hc
ap = argparse.ArgumentParser()
ap.add_argument("--diamond_revised", required=True)
ap.add_argument("--revised_pep", required=True)
ap.add_argument("--proteins", required=True)
ap.add_argument("--hc_pep_out", required=True)
ap.add_argument("--qdict_json", required=True)
args = ap.parse_args()
q_dict = hc.getting_hc_supported_by_proteins(args.diamond_revised, args.revised_pep, args.proteins, args.hc_pep_out)
q_serial = {k: [v[0], str(v[1])] for k, v in q_dict.items()}
with open(args.qdict_json, "w") as f:
    json.dump(q_serial, f)

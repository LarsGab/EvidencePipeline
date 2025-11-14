#!/usr/bin/env python3
import json, argparse
import pandas as pd 
import logging

def get_cds_classification(normal_tsv, shortened_tsv):
    """
    Compare DIAMOND search results to classify codingseq candidates as complete or incomplete.

    This function parses the DIAMOND results from normal and shortened ORFs, calculates 
    support scores, and determines whether each codingseq candidate should be classified as 
    "complete" or "incomplete" based on its best protein alignments.

    :param normal_tsv: Path to the DIAMOND results for the original TransDecoder output.
    :param shortened_tsv: Path to the DIAMOND results for the shortened ORFs.
    :return: Dictionary with codingseq IDs as keys and classifications ("complete" or "incomplete") as values.
    """
    logging.info("Analyzing incomplete ORF predictions...")

    # Define DIAMOND output column names
    header_list = ["cdsID", "proteinID", "percIdentMatches", "alignLength", "mismatches", "gapOpenings",
                   "queryStart", "queryEnd", "targetStart", "targetEnd", "eValue", "bitScore"]

    # Load DIAMOND results into pandas DataFrames
    df_normal = pd.read_csv(normal_tsv, delimiter='\t', header=None, names=header_list)
    df_shortened = pd.read_csv(shortened_tsv, delimiter='\t', header=None, names=header_list)

    # Merge results on "cdsID" and "proteinID" to compare normal and shortened ORFs
    merged_df = pd.merge(df_shortened, df_normal, on=["cdsID", "proteinID"], suffixes=('_short', '_normal'))

    # Drop unnecessary columns
    columns_to_drop = ["alignLength_short", "alignLength_normal", "mismatches_short", "mismatches_normal",
                       "gapOpenings_short", "gapOpenings_normal", "queryEnd_short", "queryEnd_normal",
                       "targetEnd_short", "targetEnd_normal", "eValue_short", "eValue_normal"]
    merged_df.drop(columns=columns_to_drop, inplace=True)

    # Add support score column
    merged_df["supportScore"] = None

    # Compute support scores
    for i, row in merged_df.iterrows():
        q_incomplete_start = row["queryStart_normal"]
        t_incomplete_start = row["targetStart_normal"]
        t_complete_start = row["targetStart_short"]
        aai_incomplete = row["percIdentMatches_normal"]
        aai_complete = row["percIdentMatches_short"]

        # Prevent division by zero by assigning a small value if AAI is zero
        if aai_complete == 0:
            aai_complete = 0.0001

        match_log = math.log(aai_incomplete / aai_complete)
        support_score = (t_complete_start - t_incomplete_start) - (q_incomplete_start - 1) + match_log**1000
        merged_df.at[i, "supportScore"] = support_score

    # Compute maximum bit score for each candidate
    merged_df["bitScore_max"] = merged_df[["bitScore_short", "bitScore_normal"]].max(axis=1)

    # Group by cdsID and analyze the best 25 alignments
    classifications = {}
    for cds_id, group in merged_df.groupby("cdsID"):
        sorted_group = group.sort_values(by="bitScore_max", ascending=False)

        # Determine if any of the top 25 alignments has a positive support score
        incomplete = any(sorted_group.head(25)["supportScore"] > 0)

        classifications[cds_id] = "incomplete" if incomplete else "complete"

    logging.info("CDS classification completed successfully.")
    return classifications


def get_optimized_pep_file(normal_pep, shortened_pep, classifications, output_path="revised_candidates.pep"):
    """
    Create a revised PEP file by incorporating classifications and shortened ORFs.

    This function reads a normal PEP file and a shortened PEP file, updates records based on
    classifications, and writes the revised records to a new file.

    :param normal_pep: Path to the original PEP file with all ORFs.
    :param shortened_pep: Path to the PEP file containing shortened ORFs.
    :param classifications: Dictionary of classifications for ORFs.
    :param output_path: Path to the output revised PEP file (default: "revised_candidates.pep").
    :return: Path to the output PEP file.
    """
    # Ensure the output directory exists
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    # Parse the shortened ORFs into a dictionary for quick access
    shortened_pep_dict = {
        record.id: (record.seq, record.description)
        for record in SeqIO.parse(shortened_pep, "fasta")
    }

    logging.info("Creating a revised PEP file using classifications...")

    with open(output_path, "w") as output:
        # Process each record in the normal PEP file
        for record in SeqIO.parse(normal_pep, "fasta"):
            orf_id = record.id

            # Determine the classification type from the record's description
            if "type:complete" in record.description:
                classification = "complete"
            elif "type:5prime_partial" in record.description:
                classification = "5prime_partial"
            elif "type:internal" in record.description:
                classification = "internal"
            elif "type:3prime_partial" in record.description:
                classification = "3prime_partial"
            else:
                classification = None

            # Handle incomplete ORFs (5prime_partial or internal)
            if classification in {"5prime_partial", "internal"}:
                if orf_id in classifications:
                    # Write the original record if classified as incomplete
                    if classifications[orf_id] == "incomplete":
                        SeqIO.write(record, output, "fasta")
                    else:
                        # Update the record with shortened sequence and adjusted description
                        seq, description = shortened_pep_dict[orf_id]
                        record.seq = seq
                        if "type:5prime_partial" in description:
                            description = description.replace("type:5prime_partial", "type:complete")
                        elif "type:internal" in description:
                            description = description.replace("type:internal", "type:3prime_partial")
                        record.description = description
                        SeqIO.write(record, output, "fasta")
                else:
                    # Write the original record if no classification is found
                    SeqIO.write(record, output, "fasta")
            else:
                # Write original records for complete or 3prime_partial ORFs
                SeqIO.write(record, output, "fasta")

    logging.info(f"Revised PEP file created successfully at: {output_path}")
    return output_path

ap = argparse.ArgumentParser()
ap.add_argument("--diamond_normal", required=True)
ap.add_argument("--diamond_short", required=True)
ap.add_argument("--transdecoder_pep", required=True)
ap.add_argument("--shortened_pep", required=True)
ap.add_argument("--revised_pep", required=True)
ap.add_argument("--classifications_json", required=True)
args = ap.parse_args()
classes = get_cds_classification(args.diamond_normal, args.diamond_short)
with open(args.classifications_json, "w") as f:
    json.dump(classes, f)
get_optimized_pep_file(args.transdecoder_pep, args.shortened_pep, classes, args.revised_pep)

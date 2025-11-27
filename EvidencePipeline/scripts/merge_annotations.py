#!/usr/bin/env python3
"""
Merge multiple GTF/GFF files into a single GFF3 file.

Two modes:

1) Full merge (mode=full)
   - Take the union of all transcripts across all inputs.
   - Define each transcript by its exon structure (seqid, strand, exon coordinates).
   - Cluster overlapping transcripts into "genes" (loci) per seqid/strand.
   - Within each locus, collapse identical transcript structures.
   - Assign new gene IDs: gene_000001, gene_000002, ...
   - Assign new transcript IDs: gene_000001.t1, gene_000001.t2, ...

2) Priority merge (mode=priority)
   - One input file is chosen as priority via --priority-file.
   - Build genes from priority file transcripts (same clustering as in full mode).
   - Keep all these priority genes.
   - For transcripts from non-priority files:
        * Discard any transcript that overlaps any priority gene (same seqid & strand, any bp overlap).
        * Keep non-overlapping transcripts, build additional genes from them.
   - Collapse identical transcripts within each locus as in full mode.

Input:
    One or more GTF/GFF files (positional arguments). Each may be GTF or GFF3.
    We assume transcripts are defined by 'exon' and/or 'CDS' features.

Output:
    A GFF3 written to stdout (or redirected by the user). Source (column 2)
    preserves the originating annotation source for each gene and child feature.

Usage examples:

    # Full union
    python merge_annotations.py --mode full in1.gff3 in2.gtf > merged_full.gff3

    # Priority merge with in1.gff3 as priority
    python merge_annotations.py --mode priority --priority-file in1.gff3 \
        in1.gff3 in2.gtf in3.gff3 > merged_priority.gff3
"""

import sys
import argparse
from collections import defaultdict
from dataclasses import dataclass, field
from typing import List, Tuple, Dict, Optional, Set


# ----------------- Attribute parsing ----------------- #

def parse_attributes_mixed(attr_str: str) -> Dict[str, str]:
    """
    Parse a GTF or GFF3 attributes column into a dict.

    Supports:
        - GFF3 style: key=value;key2=value2
        - GTF style:  key "value"; gene_id "X"; transcript_id "Y";
    """
    attrs: Dict[str, str] = {}
    if not attr_str or attr_str.strip() in {".", ""}:
        return attrs

    parts = [p for p in attr_str.strip().split(";") if p.strip()]
    for part in parts:
        part = part.strip()
        if not part:
            continue
        if "=" in part:  # GFF3-like
            key, value = part.split("=", 1)
            attrs[key.strip()] = value.strip()
        else:  # likely GTF-like: key "value"
            toks = part.split(None, 1)
            key = toks[0].strip()
            if len(toks) > 1:
                value = toks[1].strip().strip('"')
            else:
                value = ""
            attrs[key] = value
    return attrs


def attrs_to_str(attrs: Dict[str, str]) -> str:
    """Convert attributes dict back to a GFF3 attribute string."""
    if not attrs:
        return "."
    items = []
    # ID first, Parent second, others in key order
    if "ID" in attrs:
        items.append(f"ID={attrs['ID']}")
    if "Parent" in attrs:
        items.append(f"Parent={attrs['Parent']}")
    for k in sorted(attrs.keys()):
        if k in {"ID", "Parent"}:
            continue
        items.append(f"{k}={attrs[k]}")
    return ";".join(items) if items else "."


# ----------------- Data structures ----------------- #

@dataclass
class Transcript:
    """Transcript defined by its exon structure."""
    internal_id: str           # unique inside this script (includes source index)
    source_index: int          # which input file
    seqid: str
    strand: str
    source_labels: Set[str] = field(default_factory=set)
    exons: List[Tuple[int, int]] = field(default_factory=list)
    cds: List[Tuple[int, int, str]] = field(default_factory=list)

    @property
    def source(self) -> str:
        """Return a stable source string for output."""
        labels = sorted({s for s in self.source_labels if s})
        if not labels:
            return "merge"
        return ",".join(labels)

    @property
    def span(self) -> Tuple[int, int]:
        """Return (start, end) of the transcript span."""
        coords = self.exon_intervals
        starts = [s for s, _ in coords]
        ends = [e for _, e in coords]
        return min(starts), max(ends)

    @property
    def exon_intervals(self) -> List[Tuple[int, int]]:
        """
        Exon-like intervals used for structure/output. Falls back to CDS if no exons.
        """
        if self.exons:
            return self.exons
        return [(s, e) for s, e, _ in self.cds]

    @property
    def struct_key(self) -> Tuple[str, str, Tuple[Tuple[int, int], ...]]:
        """
        Structure key for deduplication:
        (seqid, strand, sorted_exons).
        """
        exons_sorted = tuple(sorted(self.exon_intervals))
        return self.seqid, self.strand, exons_sorted


@dataclass
class GeneCluster:
    """Group of transcripts representing one merged gene locus."""
    gene_id: str
    seqid: str
    strand: str
    start: int
    end: int
    transcripts: List[Transcript] = field(default_factory=list)


# ----------------- Parsing inputs ----------------- #

def parse_transcripts_from_file(path: str, source_index: int) -> List[Transcript]:
    """
    Parse one GTF/GFF file and return a list of Transcript objects.
    Transcripts are defined by exon/CDS features. We do not rely on gene features.
    """
    # internal key: (source_index, external_transcript_id)
    tx_map: Dict[str, Transcript] = {}

    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) != 9:
                continue
            seqid, source, ftype, start, end, score, strand, phase, attr_str = cols
            attrs = parse_attributes_mixed(attr_str)

            # We care about exons/CDS for structure
            if ftype not in {"exon", "CDS"}:
                continue

            # Get a transcript identifier
            tid: Optional[str] = (
                attrs.get("transcript_id")
                or attrs.get("transcriptId")
            )
            if not tid:
                parent = attrs.get("Parent")
                if parent:
                    # Parent may be comma-separated; take the first
                    tid = parent.split(",")[0].strip()

            if not tid:
                # No usable transcript ID; skip
                continue

            internal_id = f"{source_index}:{tid}"
            if internal_id not in tx_map:
                tx_map[internal_id] = Transcript(
                    internal_id=internal_id,
                    source_index=source_index,
                    seqid=seqid,
                    strand=strand,
                    source_labels={source},
                    exons=[],
                    cds=[],
                )

            tx = tx_map[internal_id]
            tx.source_labels.add(source)
            if ftype == "exon":
                tx.exons.append((int(start), int(end)))
            elif ftype == "CDS":
                tx.cds.append((int(start), int(end), phase))

    # Filter out transcripts with no exon/CDS signal
    transcripts = [t for t in tx_map.values() if t.exons or t.cds]
    return transcripts


def parse_all_inputs(input_paths: List[str]) -> List[Transcript]:
    """Parse all input files and collect all transcripts."""
    all_tx: List[Transcript] = []
    for idx, path in enumerate(input_paths):
        txs = parse_transcripts_from_file(path, source_index=idx)
        all_tx.extend(txs)
    return all_tx


# ----------------- Gene building / clustering ----------------- #

def build_gene_clusters_from_transcripts(transcripts: List[Transcript]) -> List[GeneCluster]:
    """
    Given a list of Transcript objects, cluster them into genes by overlap.

    Strategy:
        - Group transcripts by (seqid, strand).
        - Within each group, sort by start.
        - Greedy clustering: consecutive transcripts that overlap (by span) form one gene.
        - Within each gene, collapse identical transcript structures.
        - New gene IDs are assigned later by the caller.
    """
    clusters: List[GeneCluster] = []
    # group by seqid, strand
    by_chr_strand: Dict[Tuple[str, str], List[Transcript]] = defaultdict(list)
    for tx in transcripts:
        by_chr_strand[(tx.seqid, tx.strand)].append(tx)

    for (seqid, strand), tx_list in by_chr_strand.items():
        # sort by transcript start
        tx_list_sorted = sorted(tx_list, key=lambda t: t.span[0])

        current_group: List[Transcript] = []
        current_start = None
        current_end = None

        for tx in tx_list_sorted:
            t_start, t_end = tx.span
            if not current_group:
                current_group = [tx]
                current_start, current_end = t_start, t_end
            else:
                # overlap?
                if t_start <= current_end:
                    current_group.append(tx)
                    current_end = max(current_end, t_end)
                else:
                    # finalize previous group
                    clusters.append(
                        GeneCluster(
                            gene_id="",  # filled later
                            seqid=seqid,
                            strand=strand,
                            start=current_start,
                            end=current_end,
                            transcripts=current_group,
                        )
                    )
                    # start new group
                    current_group = [tx]
                    current_start, current_end = t_start, t_end

        if current_group:
            clusters.append(
                GeneCluster(
                    gene_id="",
                    seqid=seqid,
                    strand=strand,
                    start=current_start,
                    end=current_end,
                    transcripts=current_group,
                )
            )

    # Deduplicate identical transcripts within each cluster
    for cluster in clusters:
        seen_structs: Dict[Tuple[str, str, Tuple[Tuple[int, int], ...]], Transcript] = {}
        unique_txs: List[Transcript] = []
        for tx in cluster.transcripts:
            key = tx.struct_key
            if key not in seen_structs:
                seen_structs[key] = tx
                unique_txs.append(tx)
            # else: duplicate, skip
        cluster.transcripts = unique_txs

    return clusters


# ----------------- Priority overlap filtering ----------------- #

def build_priority_interval_index(priority_clusters: List[GeneCluster]):
    """
    Build an index of priority gene intervals for fast-ish overlap checks.

    Returns:
        dict[(seqid, strand)] -> list[(start, end)]
    """
    idx: Dict[Tuple[str, str], List[Tuple[int, int]]] = defaultdict(list)
    for gc in priority_clusters:
        idx[(gc.seqid, gc.strand)].append((gc.start, gc.end))

    # sort intervals by start for each key
    for key in idx:
        idx[key].sort(key=lambda x: x[0])

    return idx


def transcript_overlaps_priority(tx: Transcript,
                                 interval_index: Dict[Tuple[str, str], List[Tuple[int, int]]]) -> bool:
    """Return True if transcript overlaps any priority gene."""
    t_start, t_end = tx.span
    intervals = interval_index.get((tx.seqid, tx.strand))
    if not intervals:
        return False
    for g_start, g_end in intervals:
        if g_start > t_end:
            break
        if g_end >= t_start and g_start <= t_end:
            return True
    return False


# ----------------- Output GFF3 ----------------- #

def write_clusters_as_gff3(clusters: List[GeneCluster], out_handle):
    """
    Write gene clusters as a GFF3 to out_handle.

    Features:
        - gene
        - transcript
        - exon
        - CDS
    """
    out_handle.write("##gff-version 3\n")

    # assign stable new gene IDs
    # sort clusters for stable output: by seqid, start
    clusters_sorted = sorted(
        clusters,
        key=lambda gc: (gc.seqid, gc.start, gc.strand)
    )

    gene_counter = 1

    for gc in clusters_sorted:
        gc.gene_id = f"gene_{gene_counter:06d}"
        gene_counter += 1

        # update gene span from transcripts (in case dedup shrunk it)
        starts = [tx.span[0] for tx in gc.transcripts]
        ends = [tx.span[1] for tx in gc.transcripts]
        gc.start = min(starts)
        gc.end = max(ends)

        # gene feature
        gene_sources = sorted({
            src for tx in gc.transcripts for src in tx.source.split(",") if src
        })
        gene_source = ",".join(gene_sources) if gene_sources else "merge"

        gene_attrs = {"ID": gc.gene_id}
        gene_row = [
            gc.seqid,
            gene_source,
            "gene",
            str(gc.start),
            str(gc.end),
            ".",
            gc.strand,
            ".",
            attrs_to_str(gene_attrs),
        ]
        out_handle.write("\t".join(gene_row) + "\n")

        # sort transcripts by span
        txs_sorted = sorted(gc.transcripts, key=lambda t: t.span[0])

        tx_counter = 1
        for tx in txs_sorted:
            tx_id = f"{gc.gene_id}.t{tx_counter}"
            tx_counter += 1
            t_start, t_end = tx.span

            tx_attrs = {"ID": tx_id, "Parent": gc.gene_id}
            tx_row = [
                tx.seqid,
                tx.source,
                "transcript",
                str(t_start),
                str(t_end),
                ".",
                tx.strand,
                ".",
                attrs_to_str(tx_attrs),
            ]
            out_handle.write("\t".join(tx_row) + "\n")

            # exons: sorted by start
            exons_sorted = sorted(tx.exon_intervals)
            exon_num = 1
            for s, e in exons_sorted:
                exon_id = f"{tx_id}.exon{exon_num}"
                exon_num += 1
                exon_attrs = {"ID": exon_id, "Parent": tx_id}
                exon_row = [
                    tx.seqid,
                    tx.source,
                    "exon",
                    str(s),
                    str(e),
                    ".",
                    tx.strand,
                    ".",
                    attrs_to_str(exon_attrs),
                ]
                out_handle.write("\t".join(exon_row) + "\n")

            if tx.cds:
                cds_sorted = sorted(tx.cds, key=lambda item: item[0])
                cds_num = 1
                for s, e, phase in cds_sorted:
                    cds_id = f"{tx_id}.cds{cds_num}"
                    cds_num += 1
                    cds_attrs = {"ID": cds_id, "Parent": tx_id}
                    cds_row = [
                        tx.seqid,
                        tx.source,
                        "CDS",
                        str(s),
                        str(e),
                        ".",
                        tx.strand,
                        phase if phase else ".",
                        attrs_to_str(cds_attrs),
                    ]
                    out_handle.write("\t".join(cds_row) + "\n")


# ----------------- Main logic ----------------- #

def main():
    ap = argparse.ArgumentParser(
        description="Merge GTF/GFF annotations into a single GFF3 file."
    )
    ap.add_argument(
        "--mode",
        choices=["full", "priority"],
        required=True,
        help="Merge mode: 'full' for union; 'priority' to respect one priority file."
    )
    ap.add_argument(
        "--priority-file",
        help="Path to the priority annotation file (required for --mode priority)."
    )
    ap.add_argument(
        "inputs",
        nargs="+",
        help="Input GTF/GFF files."
    )

    args = ap.parse_args()

    if args.mode == "priority":
        if not args.priority_file:
            ap.error("--priority-file is required when --mode priority is used.")
        if args.priority_file not in args.inputs:
            ap.error("--priority-file must be one of the input files.")

    # Parse all transcripts
    transcripts = parse_all_inputs(args.inputs)

    if not transcripts:
        sys.stderr.write("No transcripts (exon/CDS features) found in inputs.\n")
        return 1

    if args.mode == "full":
        clusters = build_gene_clusters_from_transcripts(transcripts)
        write_clusters_as_gff3(clusters, sys.stdout)

    else:  # priority mode
        priority_index = args.inputs.index(args.priority_file)

        # split into priority vs non-priority transcripts
        priority_tx = [t for t in transcripts if t.source_index == priority_index]
        other_tx = [t for t in transcripts if t.source_index != priority_index]

        # build genes from priority transcripts
        priority_clusters = build_gene_clusters_from_transcripts(priority_tx)

        # index priority intervals
        interval_index = build_priority_interval_index(priority_clusters)

        # filter other transcripts by overlap with priority genes
        kept_other_tx = [
            t for t in other_tx
            if not transcript_overlaps_priority(t, interval_index)
        ]

        # build genes from kept non-priority transcripts
        other_clusters = build_gene_clusters_from_transcripts(kept_other_tx)

        # combined clusters: priority first, then others
        all_clusters = priority_clusters + other_clusters

        write_clusters_as_gff3(all_clusters, sys.stdout)

    return 0


if __name__ == "__main__":
    sys.exit(main())

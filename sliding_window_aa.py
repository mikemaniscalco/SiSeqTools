#!/usr/bin/env python3
"""
sliding_window_aa.py

Author: Michael Maniscalco
Date: 2025-09-18
License: MIT
---------------------

Description:
    This script scans protein FASTA sequences with a sliding window
    to measure amino acid composition (% Lysine, % Serine).
    Only windows meeting thresholds are reported:
        - ≥ 10% Lysine, or
        - ≥ 18% Serine

    To reduce redundant output, overlapping windows with
    identical length, identical percentages, and overlapping
    coordinates are collapsed into a single row (keeping the first
    occurrence). If overlapping windows differ slightly, only the
    "best" one (higher %K or %S) is retained.

    Default step size is 5 amino acids.

Usage:
    python sliding_window_aa.py -i input.fasta -w 100-2000 -o output.tsv -s 1 -t 2

Arguments:
    -i         Input FASTA file of amino acid sequences
    -w/--window Window size: either a single integer (e.g. 100),
                or a range MIN-MAX (e.g. 100-2000)
    -s/--step   Step size for sliding window (default: 5)
    -t/--threads Number of CPU threads (default: 1)
    -o/--output Output file (TSV format)
"""

import argparse
import multiprocessing as mp
from Bio import SeqIO

THRESH_K = 10
THRESH_S = 18
AA_1 = "K" # first of 2 amino acids to be checking the % of 
AA_2 = "S" # second of 2 amino acids to be checking the % of 


def compute_window(seq, start, window):
    """Compute %K and %S for a given window of a sequence."""
    subseq = seq[start:start + window]
    length = len(subseq)
    pct_k = subseq.count(AA_1) / length * 100 if length > 0 else 0
    pct_s = subseq.count(AA_2) / length * 100 if length > 0 else 0
    return pct_k, pct_s


def process_record(args):
    """
    Slide windows over a sequence for all window sizes in [min_w, max_w].
    Return only consolidated windows:
      - If overlapping windows have identical length, %K, and %S → keep only the first.
      - If multiple overlapping windows differ → keep the one with highest %K+%S.
    """
    record, min_w, max_w, step = args
    seq = str(record.seq).upper()
    rec_id = record.id
    raw_results = []

    # collect all passing windows first
    for window in range(min_w, max_w + 1, step if min_w != max_w else 1):
        for start in range(0, len(seq) - window + 1, step):
            subseq = seq[start:start + window]
            length = len(subseq)
            if length == 0:
                continue
            pct_k = subseq.count(AA_1) / length * 100
            pct_s = subseq.count(AA_2) / length * 100
            if pct_k >= THRESH_K or pct_s >= THRESH_S:
                raw_results.append((window, start, start + window, pct_k, pct_s))

    if not raw_results:
        return []

    # sort by start so overlapping regions are adjacent
    raw_results.sort(key=lambda x: (x[1], x[2]))

    consolidated = []
    prev = None

    for window, start, end, pct_k, pct_s in raw_results:
        if prev is None:
            prev = (window, start, end, pct_k, pct_s)
            continue

        pw, ps, pe, pk, ps_aa = prev

        # check if current overlaps with previous
        if start < pe:
            # identical percentages & length → keep only the first (prev)
            if window == pw and abs(pk - pct_k) < 1e-9 and abs(ps_aa - pct_s) < 1e-9:
                continue
            else:
                # choose the "better" window by higher sum %K+%S
                if (pct_k + pct_s) > (pk + ps_aa):
                    prev = (window, start, end, pct_k, pct_s)
                # else keep prev
        else:
            # no overlap, flush prev
            consolidated.append(prev)
            prev = (window, start, end, pct_k, pct_s)

    # don’t forget to append the last one
    if prev is not None:
        consolidated.append(prev)

    # attach record id to results
    return [(rec_id, w, s, e, pk, ps) for (w, s, e, pk, ps) in consolidated]



def filter_windows(windows):
    """
    Collapse redundant overlapping windows:
      - If length, %K, and %S are identical, keep only the first.
      - If overlapping windows differ, keep the one with
        higher %K or %S.
    Works per sequence.
    """
    filtered = []
    windows.sort(key=lambda x: (x[0], x[2]))  # sort by seq_id, start

    current = None
    for w in windows:
        if current is None:
            current = w
            continue

        same_seq = w[0] == current[0]
        same_len = w[1] == current[1]
        same_pct = abs(w[4] - current[4]) < 1e-6 and abs(w[5] - current[5]) < 1e-6
        overlap = w[2] < current[3]

        if same_seq and same_len and same_pct and overlap:
            # redundant identical window → skip
            continue
        elif same_seq and overlap:
            # overlapping but not identical → keep best (%K or %S higher)
            if (w[4] > current[4]) or (w[5] > current[5]):
                current = w
            continue
        else:
            filtered.append(current)
            current = w

    if current is not None:
        filtered.append(current)

    return filtered


def main():
    parser = argparse.ArgumentParser(description="Sliding window amino acid analysis.")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file")
    parser.add_argument("-o", "--output", required=True, help="Output TSV file")
    parser.add_argument("-w", "--window", type=str, required=True, help="Window size as N or MIN-MAX (e.g. 100-400)")
    parser.add_argument("-s", "--step", type=int, default=5, help="Step size (default=5)")
    parser.add_argument("-t", "--threads", type=int, default=1, help="Number of threads")
    args = parser.parse_args()

    # Parse window argument
    if "-" in args.window:
        min_w, max_w = map(int, args.window.split("-"))
    else:
        min_w = max_w = int(args.window)

    records = list(SeqIO.parse(args.input, "fasta"))
    tasks = [(record, min_w, max_w, args.step) for record in records]

    with open(args.output, "w", buffering=1) as f:
        f.write("seq_id\twindow\tstart\tend\tpercent_K\tpercent_S\n")

        # Use imap_unordered so results stream in as they're ready
        with mp.Pool(processes=args.threads) as pool:
            for result_list in pool.imap_unordered(process_record, tasks, chunksize=1):
                for rec_id, window, start, end, pct_k, pct_s in result_list:
                    f.write(
                        f"{rec_id}\t{window}\t{start}\t{end}\t{pct_k:.4f}\t{pct_s:.4f}\n"
                    )


if __name__ == "__main__":
    main()


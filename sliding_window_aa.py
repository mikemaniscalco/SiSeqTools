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
    Return list of windows that pass thresholds (≥10% K or ≥18% S).
    """
    record, min_w, max_w, step = args
    seq = str(record.seq).upper()
    rec_id = record.id
    results = []

    for window in range(min_w, max_w + 1, step if min_w != max_w else 1):
        for start in range(0, len(seq) - window + 1, step):
            pct_k, pct_s = compute_window(seq, start, window)
            if pct_k >= THRESH_K or pct_s >= THRESH_S:
                results.append((rec_id, window, start, start + window,
                                pct_k, pct_s))
    return results


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
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", dest="fasta", help="Input FASTA file")
    parser.add_argument("-w", "--window", type=str, required=True,
                        help="Window size as N or MIN-MAX (e.g. 200-400)")
    parser.add_argument("-s", "--step", type=int, default=5,
                        help="Step size for sliding window (default: 5)")
    parser.add_argument("-t", "--threads", type=int, default=1,
                        help="Number of threads (default: 1)")
    parser.add_argument("-o", "--output", type=str, required=True,
                        help="Output file (TSV)")
    args = parser.parse_args()

    # Parse window argument
    if "-" in args.window:
        min_w, max_w = map(int, args.window.split("-"))
    else:
        min_w = max_w = int(args.window)

    records = list(SeqIO.parse(args.fasta, "fasta"))

    # Prepare args for multiprocessing
    tasks = [(rec, min_w, max_w, args.step) for rec in records]

    if args.threads > 1:
        with mp.Pool(processes=args.threads) as pool:
            all_results = pool.map(process_record, tasks)
    else:
        all_results = [process_record(task) for task in tasks]

    # flatten results
    all_results = [w for rec_res in all_results for w in rec_res]

    # filter redundant results
    filtered_results = filter_windows(all_results)

    # write output
    with open(args.output, "w") as f:
        f.write("seq_id\twindow\tstart\tend\tpercent_K\tpercent_S\n")
        for rec_id, window, start, end, pct_k, pct_s in filtered_results:
            f.write(f"{rec_id}\t{window}\t{start}\t{end}\t{pct_k:.1f}\t{pct_s:.1f}\n")


if __name__ == "__main__":
    main()


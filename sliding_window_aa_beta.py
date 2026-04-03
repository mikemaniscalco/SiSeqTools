#!/usr/bin/env python3
"""
sliding_window_aa.py

Author: Michael Maniscalco
Date: 2026-03-23
License: MIT
---------------------

Description:
    This script scans protein FASTA sequences with a variable-length 
    sliding window to measure specific amino acid compositions.
    
    ALGORITHM OVERVIEW:
    To achieve maximum performance, this script utilizes NumPy Prefix Sum 
    (cumulative sum) arrays. Instead of repeatedly slicing and counting 
    strings, the script precomputes the sequence mathematically once. 
    This reduces the computational complexity of evaluating any window 
    from O(N) to O(1), allowing rapid evaluation of massive datasets 
    even at a step size of 1.

    THRESHOLDS:
    Target amino acids and their minimum percentage thresholds are provided 
    dynamically via the command line. Only windows meeting AT LEAST ONE of 
    the provided thresholds are reported.

    CONSOLIDATION LOGIC:
    To reduce redundant output, overlapping windows (defined as > 25% 
    overlap fraction) are collapsed into a single row. The script 
    prioritizes the "best" window based on the following tie-breakers:
        1. Windows meeting the highest number of amino acid thresholds.
        2. The longer of the two windows.

Usage:
    python sliding_window_aa.py -i input.fasta -w 100-2000 -a K=10 S=18 -o output.tsv -s 1 -t 2

Arguments:
    -i/--input       Input FASTA file of amino acid sequences
    -w/--window      Window size: either a single integer (e.g. 100),
                     or a range MIN-MAX (e.g. 100-2000)
    -a/--amino-acids Target amino acids and thresholds as AA=PCT (e.g., K=10 S=18)
    -s/--step        Step size for sliding window (default: 1)
    -t/--threads     Number of CPU threads (default: 1)
    -o/--output      Output file (TSV format)
"""

import argparse
import multiprocessing as mp
import numpy as np
from Bio import SeqIO

def process_record(args):
    """
    Slide windows over a sequence using Prefix Sums for O(1) evaluation.
    Consolidate overlapping windows based on dynamic priority rules.
    """
    record, min_w, max_w, step, criteria = args
    seq = str(record.seq).upper()
    rec_id = record.id
    seq_len = len(seq)
    raw_results = []

    if seq_len < min_w:
        return []

    # --- 1. PRECOMPUTE PREFIX SUMS ---
    seq_bytes = np.frombuffer(seq.encode('ascii'), dtype=np.uint8)
    prefix_sums = {}
    
    for aa in criteria:
        is_aa = seq_bytes == ord(aa)
        p_sums = np.zeros(seq_len + 1, dtype=np.int32)
        p_sums[1:] = np.cumsum(is_aa)
        prefix_sums[aa] = p_sums

    # --- 2. FAST WINDOW SCANNING ---
    for window in range(min_w, max_w + 1, step if min_w != max_w else 1):
        for start in range(0, seq_len - window + 1, step):
            end = start + window
            
            pcts = {}
            met_count = 0
            
            # O(1) calculation for every targeted amino acid
            for aa, threshold in criteria.items():
                count = prefix_sums[aa][end] - prefix_sums[aa][start]
                pct = (count / window) * 100
                pcts[aa] = pct
                
                if pct >= threshold:
                    met_count += 1
            
            # If at least one threshold is met (OR logic)
            if met_count > 0:
                raw_results.append((window, start, end, pcts, met_count))

    if not raw_results:
        return []

    # --- 3. CONSOLIDATION LOGIC ---
    raw_results.sort(key=lambda x: (x[1], x[2]))

    consolidated = []
    prev = None

    for window, start, end, pcts, met_count in raw_results:
        if prev is None:
            prev = (window, start, end, pcts, met_count)
            continue

        pw, ps, pe, ppcts, pmet = prev

        overlap_len = max(0, min(end, pe) - max(start, ps))
        overlap_fraction = overlap_len / min(end - start, pe - ps)

        if overlap_fraction > 0.25:
            # Priority 1: Meets more thresholds simultaneously
            if met_count > pmet:
                prev = (window, start, end, pcts, met_count)
            # Priority 2: Tie-breaker - longer window wins
            elif met_count == pmet:
                if window > pw:
                    prev = (window, start, end, pcts, met_count)
        else:
            consolidated.append(prev)
            prev = (window, start, end, pcts, met_count)

    if prev is not None:
        consolidated.append(prev)

    # Return tuples of (rec_id, window, start, end, dictionary_of_percentages)
    return [(rec_id, w, s, e, pcts) for (w, s, e, pcts, _) in consolidated]


def main():
    parser = argparse.ArgumentParser(description="Sliding window amino acid analysis.")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file")
    parser.add_argument("-o", "--output", required=True, help="Output TSV file")
    parser.add_argument("-w", "--window", type=str, required=True, help="Window size as N or MIN-MAX (e.g. 100-400)")
    parser.add_argument("-a", "--amino-acids", nargs='+', required=True, help="Target amino acids and thresholds as AA=PCT (e.g., K=10 S=18)")
    parser.add_argument("-s", "--step", type=int, default=1, help="Step size (default=1)")
    parser.add_argument("-t", "--threads", type=int, default=1, help="Number of threads")
    args = parser.parse_args()

    # Parse window argument
    if "-" in args.window:
        min_w, max_w = map(int, args.window.split("-"))
    else:
        min_w = max_w = int(args.window)

    # Parse dynamic amino acid criteria into a dictionary
    criteria = {}
    for item in args.amino_acids:
        aa, pct = item.split('=')
        criteria[aa.upper()] = float(pct)

    # Note: Loading all records into a list uses significant RAM for massive files.
    records = list(SeqIO.parse(args.input, "fasta"))
    
    # Pass the criteria dictionary into the mp.Pool tasks
    tasks = [(record, min_w, max_w, args.step, criteria) for record in records]
    
    # Keep track of the order of amino acids for the TSV header
    aa_keys = list(criteria.keys())

    with open(args.output, "w", buffering=1) as f:
        # Dynamically build the TSV header
        header_cols = ["seq_id", "window", "start", "end"] + [f"percent_{aa}" for aa in aa_keys]
        f.write("\t".join(header_cols) + "\n")

        with mp.Pool(processes=args.threads) as pool:
            for result_list in pool.imap_unordered(process_record, tasks, chunksize=1):
                for rec_id, window, start, end, pcts in result_list:
                    # Dynamically build the row values in the correct order
                    row_vals = [rec_id, str(window), str(start), str(end)]
                    row_vals.extend([f"{pcts[aa]:.4f}" for aa in aa_keys])
                    
                    f.write("\t".join(row_vals) + "\n")


if __name__ == "__main__":
    main()

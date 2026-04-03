#!/usr/bin/env python3

import os
import subprocess

# download and pip install https://services.healthtech.dtu.dk/services/SignalP-6.0/

def run_signalp(fasta_file, output_dir, organism_type='other', run_mode='fast'):
    """
    Runs the SignalP 6.0 command-line tool using subprocess.

    Args:
        fasta_file (str): Path to the input FASTA file with protein sequences.
        output_dir (str): Directory where output files will be saved.
        organism_type (str): Organism type, e.g., 'eukarya' or 'other' (default).
        run_mode (str): Prediction mode, 'fast' (default), 'slow', or 'slow-sequential'.
    """
    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Construct the command
    command = [
        'signalp6',
        '--fastafile', fasta_file,
        '--organism', organism_type,
        '--output_dir', output_dir,
        '--format', 'txt', # 'txt' produces a tabular .gff file
        '--mode', run_mode
    ]
    
    print(f"Running command: {' '.join(command)}")

    try:
        # Execute the command
        result = subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        print("SignalP search completed successfully.")
        return result.stdout
    except subprocess.CalledProcessError as e:
        print(f"Error during SignalP run: {e.stderr}")
        return None

# Example usage:
input_fasta = 'my_sequences.fasta'
output_folder = 'signalp_results'

# Example ti parse the generated output files (.gff or summary files) in Python
run_signalp(input_fasta, output_folder, organism_type='eukarya', run_mode='fast')

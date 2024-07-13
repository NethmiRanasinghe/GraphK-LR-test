import argparse
import sys
import numpy as np
import subprocess
from Bio import SeqIO
from collections import defaultdict

def read_tsv_to_clusters(tsv_file, sequence_file):
    # Determine file format
    file_format = "fasta" if sequence_file.endswith(('.fasta', '.fa')) else "fastq"

    # Count total reads
    if file_format == "fastq":
        count_command = f"grep -c '^@' {sequence_file}"
    else:  # fasta
        count_command = f"grep -c '^>' {sequence_file}"

    total_sequences = int(subprocess.run(
        count_command,
        shell=True,
        capture_output=True,
        text=True,
        check=True
    ).stdout.strip())

    # Initialize the clusters array
    clusters = np.full(total_sequences, -1, dtype=int)

    # Read TSV file
    readId_bin_dict = {}
    with open(tsv_file, 'r') as file:
        next(file)  # Skip header
        for line in file:
            read_id, bin_str = line.strip().split('\t')
            bin_index = int(bin_str)
            readId_bin_dict[read_id] = bin_index


    # Process sequence file
    for i, record in enumerate(SeqIO.parse(sequence_file, file_format)):
        if record.id in readId_bin_dict:
            clusters[i] = readId_bin_dict[record.id]

    return clusters

def bins_to_npy(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()

    # Extract numerical values from each line
    values = [int(line.strip().split('-')[1]) for line in lines]

    # Convert the list to a NumPy array
    array = np.array(values)

    return array
    
    

def other_function():
    # Other function implementation
    print("""The name of the binning tool is unrecongnized. Please enter one of these - oblr, lrbinner, metabcclr, semibin2""")
    

    exit()
    

def main():
    parser = argparse.ArgumentParser(description="Binning tool processing script")
    parser.add_argument("-i", "--initial_tool_results", help="Path to the initial binning results file")
    parser.add_argument("-r", "--reads_file", help="Path to the reads file (fastq, fasta)")
    parser.add_argument("-o", "--out_folder", help="Path to the output folder")
    parser.add_argument("-t", "--binning_tool", help="Name of the binning tool used (oblr, lrbinner, metabcclr, semibin2)")
    args = parser.parse_args()

            
    tool = args.binning_tool

    if tool == "oblr":
        clusters = np.load(args.initial_tool_results)['classes']
        np.save(args.out_folder + '/initial_bins.npy', clusters)
        
    elif tool == "lrbinner":
        clusters = np.loadtxt(args.initial_tool_results, dtype=int) 
        np.save(args.out_folder + '/initial_bins.npy', clusters)

    elif tool == "metabcclr":
        clusters = bins_to_npy(args.initial_tool_results)
        np.save(args.out_folder + '/initial_bins.npy', clusters)

    elif tool == "semibin2":
        clusters = read_tsv_to_clusters(args.initial_tool_results, args.fastq_file)
        np.save(args.out_folder + '/initial_bins.npy', clusters)

    else:
        other_function()
        
    print(f'Results parsed succesfully! You can find the results at {args.out_folder}/initial_bins.npy')

if __name__ == "__main__":
    main()
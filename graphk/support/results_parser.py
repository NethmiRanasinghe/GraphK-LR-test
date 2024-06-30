import argparse
import sys
import numpy as np
import subprocess
from Bio import SeqIO

def read_tsv_to_clusters(tsv_file, fastq_file):

    total_reads = subprocess.run(
            f"cat {fastq_file} | wc -l | awk '{{print $1/4}}'",
            shell=True,
            capture_output=True,
            text=True,
            check=True
        )
    
    # fastq_handle = open(fastq_file, 'r')
    # total_reads = sum(1 for _ in SeqIO.parse(fastq_handle, "fastq"))
    # fastq_handle.seek(0)  # Reset file pointer to the beginning
    
    # Initialize the clusters array with -1
    clusters = np.full(total_reads, -1, dtype=int)

    readId_bin_dict = {}

    with open(tsv_file, 'r') as file:
            file.readline()
            for line in file:
                read_id, bin = line.strip().split('\t')
                readId_bin_dict[read_id] = bin

    i = 0
    for record in SeqIO.parse(fastq_file, "fastq"):
        if record.id in readId_bin_dict:
            clusters[i] = bin
            i += 1   

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
    print("""This parser does not support for other tools.
          Please convert the results file into a .npy file with respective bin of each seq in order of the original reads file.
          Bin number should start from 0. """)
    

    exit()
    

def main():
    parser = argparse.ArgumentParser(description="Binning tool processing script")
    parser.add_argument("initial_tool_results", help="Path to the input txt file")
    parser.add_argument("fastq_file", help="Path to the input fastq file")
    parser.add_argument("out_folder", help="Path to the output folder")
    args = parser.parse_args()

    print("Which tool did you use for initial binning?")
    print("1) OBLR")
    print("2) LRBinner")
    print("3) MetaBCC")
    print("4) semibin2")
    print("5) other")

    while True:
        try:
            choice = int(input("Enter your choice (1-5): "))
            if 1 <= choice <= 5:
                break
            else:
                print("Please enter a number between 1 and 5.")
        except ValueError:
            print("Please enter a valid number.")

    if choice == 1:
        clusters = np.load(args.initial_tool_results)['classes']
        np.save(args.out_folder + '/initial_bins.npy', clusters)
        
    elif choice == 2:
        clusters = np.loadtxt(args.initial_tool_results, dtype=int) 
        np.save(args.out_folder + 'initial_bins.npy', clusters)

    elif choice == 3:
        clusters = bins_to_npy(args.initial_tool_results)
        np.save(args.out_folder + 'initial_bins.npy', clusters)

    elif choice == 4:
        clusters = read_tsv_to_clusters(args.initial_tool_results, args.fastq_file)
        np.save(args.out_folder + 'initial_bins.npy', clusters)

    else:
        other_function()

if __name__ == "__main__":
    main()
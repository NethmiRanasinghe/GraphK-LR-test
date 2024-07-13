import argparse
import scipy.special
import csv
from Bio import SeqIO

from tabulate import tabulate

def run(out_dir, truthfile, reads_file):


    # Get paths to binning result and ground truth
    binned_file = out_dir+'/final_refined_bins.tsv'
    ground_truth_file = truthfile

    print("\nStarting evaluate.py...")

    # Get the number of bins from the ground truth
    #---------------------------------------------------------
    ground_truth_n_bins = 0

    all_ground_truth_bins_list = []

    with open(ground_truth_file) as csvfile:
        readCSV = csv.reader(csvfile, delimiter='\t')
        for row in readCSV:
            all_ground_truth_bins_list.append(row[1])
            
    ground_truth_bins_list = list(set(all_ground_truth_bins_list))
    ground_truth_n_bins = len(ground_truth_bins_list)

    print("\nNumber of bins available in the ground truth: ", ground_truth_n_bins)


    # Get the ground truth
    #----------------------------
    ground_truth_bins = [[] for x in range(ground_truth_n_bins)]

    ground_truth_count = 0
    ground_truth_bins_1 = {}

    with open(ground_truth_file) as contig_bins:
        if ground_truth_file.endswith('.tsv'):
            readCSV = csv.reader(contig_bins, delimiter='\t')
        elif ground_truth_file.endswith('.csv'):
            readCSV = csv.reader(contig_bins, delimiter=',')
        
        for row in readCSV:
            ground_truth_count += 1
            contig = row[0]
            bin_num = ground_truth_bins_list.index(row[1])
            ground_truth_bins[bin_num].append(contig)
            ground_truth_bins_1[contig] = bin_num

    print("Number of contigs available in the ground truth: ", ground_truth_count)

    # Get the number of bins from the initial binning result
    #---------------------------------------------------------
    n_bins = 0
    binned_count = 0

    all_bins_list = []

    with open(binned_file, 'r') as binned_txt:
        lines = binned_txt.readlines()
        
        for row in lines:
            all_bins_list.append(int(row.strip()))
            binned_count += 1 
            
    bins_list = list(set(all_bins_list))
    n_bins = len(bins_list)

    print("Number of bins available in the binning result: ", n_bins)


    # Get initial binning result
    #----------------------------
    bins = [[] for x in range(n_bins)]

    bins_1 = {}
    binned_contigs = []

    file_extension = reads_file.split('.')[-1].lower()
    
    if file_extension in ['fasta', 'fa']:
        file_format = 'fasta'
    elif file_extension in ['fastq', 'fq']:
        file_format = 'fastq'
    else:
        raise ValueError("Unsupported file format. Please provide a FASTA or FASTQ file.")
    
    with open(reads_file, 'r') as file:
        i = 0
        for record in SeqIO.parse(file, file_format):
            contig = record.id

            bin_num = all_bins_list[i]
            bins[bin_num].append(contig)
            bins_1[contig] = bin_num
            binned_contigs.append(contig)
            i += 1
            

    print("Number of contigs available in the binning result: ", binned_count)


    # Functions to determine precision, recall, F1-score and ARI
    #------------------------------------------------------------

    # Get precicion
    def getPrecision(mat, k, s, total):
        sum_k = 0
        for i in range(k):
            max_s = 0
            for j in range(s):
                if mat[i][j] > max_s:
                    max_s = mat[i][j]
            sum_k += max_s
        return sum_k/total

    # Get recall
    def getRecall(mat, k, s, total, unclassified):
        sum_s = 0
        for i in range(s):
            max_k = 0
            for j in range(k):
                if mat[j][i] > max_k:
                    max_k = mat[j][i]
            sum_s += max_k
        return sum_s/(total+unclassified)

    # Get ARI
    def getARI(mat, k, s, N):
        t1 = 0    
        for i in range(k):
            sum_k = 0
            for j in range(s):
                sum_k += mat[i][j]
            t1 += scipy.special.binom(sum_k, 2)
            
        t2 = 0
        for i in range(s):
            sum_s = 0
            for j in range(k):
                sum_s += mat[j][i]
            t2 += scipy.special.binom(sum_s, 2)
            
        t3 = t1*t2/scipy.special.binom(N, 2)
        
        t = 0
        for i in range(k):
            for j in range(s):
                t += scipy.special.binom(mat[i][j], 2)
            
        ari = (t-t3)/((t1+t2)/2-t3)
        return ari

    # Get F1-score
    def getF1(prec, recall):
        return 2*prec*recall/(prec+recall)


    # Determine precision, recall, F1-score and ARI for binning result
    #------------------------------------------------------------------

    total_binned = 0

    bins_species = [[0 for x in range(ground_truth_n_bins)] for y in range(n_bins)]

    for i in bins_1:
        if i in ground_truth_bins_1:
            total_binned += 1
            bins_species[bins_1[i]][ground_truth_bins_1[i]] += 1


    print("Number of contigs available in the binning result that are present in the ground truth:", total_binned)
    print(ground_truth_count, (ground_truth_count-total_binned))

    #print()
    #print("Ground truth bin labels:")

    #for i in range(len(ground_truth_bins_list)):
    #    print(i+1, ground_truth_bins_list[i])

    #print()
    #print("KxS Matrix:")
    #print(tabulate(bins_species))
    # print()


    my_precision = getPrecision(bins_species, n_bins, ground_truth_n_bins, total_binned)
    my_recall = getRecall(bins_species, n_bins, ground_truth_n_bins, total_binned, (ground_truth_count-total_binned))
    my_ari = getARI(bins_species, n_bins, ground_truth_n_bins, total_binned)
    my_f1 = getF1(my_precision, my_recall)


    print("Evaluation Results:")
    print("Precision =", my_precision*100)
    print("Recall =", my_recall*100)
    print("F1-score =", my_f1*100)
    print("ARI =", my_ari*100)

    print()


if __name__ == "__main__":
    main()
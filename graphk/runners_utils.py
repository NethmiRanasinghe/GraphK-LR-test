import os
import pickle
import logging
import subprocess
from Bio import SeqIO
import numpy as np
from .steps import step1, step2, step3, step4
from .support import evaluate

logger = logging.getLogger('GraphKLR')

class Checkpointer():
    def __init__(self, checkpoint_path, _load_to_resume=False):
        self.cpath = checkpoint_path
        if _load_to_resume and os.path.isfile(self.cpath):
            self.completed = pickle.load(open(self.cpath, "rb"))
        else:
            self.completed = {}

    def should_run_step(self, stage, params):
        if stage not in self.completed:
            return True
        if self.completed[stage] != params:
            return True
        return False


    def log(self, stage, params):
        self.completed[stage] = params
        ps, pc = [int(x) for x in stage.split("_")][:2]

        # remove all child stages
        for s in list(self.completed.keys()):
            p, c = [int(x) for x in s.split("_")][:2]
            # is stage is seen after current checkpoint
            if p > ps:
                del self.completed[s]

        self._save()


    def _save(self):
        pickle.dump(self.completed, open(self.cpath, "wb+"))


    def __str__(self):
        return str(self.completed)

# Define the directory where binaries are installed
# BIN_DIR = os.path.join('bin')
BIN_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'bin')
OBLR_UTILS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'oblr_utils')

def rename_reads(exp_dir, fastq_file):
    seqtk_path = os.path.join(BIN_DIR, 'seqtk')
    subprocess.run(f"{seqtk_path} rename {fastq_file} read_ | {seqtk_path} seq -A > {exp_dir}/reads.fasta", shell=True, check=True)

def obtain_read_ids(exp_dir):
    subprocess.run(f"grep '>' {exp_dir}/reads.fasta > {exp_dir}/read_ids", shell=True, check=True)
    
    
def run_seq2vec(exp_dir):
    #seq2vec_path = os.path.join(BIN_DIR, 'seq2vec')
    #subprocess.run(f"{seq2vec_path} -k 4 -o {exp_dir}/4mers -f {exp_dir}/reads.fasta", shell=True, check=True)
    subprocess.run(f"kmertools comp oligo -p csv -i {exp_dir}/reads.fasta -o {exp_dir}/4mers -k 4", shell=True, check=True)
    

def run_seq2covvec(exp_dir, fastq_file):
    #subprocess.run(f"python {BIN_DIR}/seq2covvec/seq2covvec.py -k 16 -o {exp_dir}/16mers -r {fastq_file}", shell=True, check=True)
    subprocess.run(f"kmertools cov -s 10 -c 32 -t 8 -i {fastq_file} -o {exp_dir}/16mers -k 16", shell=True, check=True)
    

def create_overlaps(exp_dir):
    subprocess.run(f"bash {OBLR_UTILS_DIR}/buildgraph_with_chunks.sh -r {exp_dir}/reads.fasta -c 250000 -o {exp_dir}/", shell=True, check=True)
    

def run_step1(in_file, exp_dir, out_dir):
    step1.run(in_file, exp_dir, out_dir)
    
    
def run_step2(exp_dir, out_dir):
    step2.run(exp_dir, out_dir)

def run_step3(exp_dir, out_dir):
    step3.run(exp_dir, out_dir)

def run_step4(exp_dir, out_dir, epochs, fastq_file):
    # subprocess.run(f"awk 'NR%4==1 {{print substr($1,2)}}' {fastq_file} > {exp_dir}/reads_original_ids", shell=True, check=True)
    step4.run(exp_dir, out_dir, epochs)
    
def run_eval(out_dir, groundtruth,):
    evaluate.run(out_dir, groundtruth)

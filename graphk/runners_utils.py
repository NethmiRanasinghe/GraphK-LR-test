import os
import pickle
import logging
import subprocess
from Bio import SeqIO
import numpy as np
from .steps import step1, step2, step3, step4, create_graph
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
    subprocess.run(f"{seqtk_path} seq -A {fastq_file} | {seqtk_path} rename - read_ > {exp_dir}/reads.fasta", shell=True, check=True)

def obtain_read_ids(exp_dir):
    subprocess.run(f"grep '>' {exp_dir}/reads.fasta > {exp_dir}/read_ids", shell=True, check=True)
    
def run_seq2vec(exp_dir):
    subprocess.run(f"kmertools comp oligo -p csv -i {exp_dir}/reads.fasta -o {exp_dir}/4mers -k 4", shell=True, check=True)

def run_seq2covvec(exp_dir, fastq_file):
    subprocess.run(f"kmertools cov -s 10 -c 32 -t 8 -i {fastq_file} -o {exp_dir}/16mers -k 16", shell=True, check=True)

def create_overlaps(exp_dir):
    subprocess.run(f"bash {OBLR_UTILS_DIR}/buildgraph_with_chunks.sh -r {exp_dir}/reads.fasta -c 250000 -o {exp_dir}/", shell=True, check=True)
    
def run_create_graph(exp_dir):
    create_graph.run(exp_dir)

def run_step1(exp_dir, in_file, out_dir):
    step1.run(exp_dir, in_file, out_dir)
    
def run_step2(exp_dir, out_dir):
    step2.run(exp_dir, out_dir)

def run_step3(exp_dir, out_dir):
    step3.run(exp_dir, out_dir)

def run_step4(exp_dir, out_dir, epochs):
    step4.run(exp_dir, out_dir, epochs)
    
def run_eval(out_dir, groundtruth, fastq_file):
    evaluate.run(out_dir, groundtruth, fastq_file)

import os
import logging
from .runners_utils import *

logger = logging.getLogger('GraphKLR')

def run_pipeline(args):

    exp_dir = args.exp_dir
    in_file = args.in_file
    fastq_file = args.fastq_file
    epochs = args.epochs
    resume = args.resume
    groundtruth = args.groundtruth

    # Ensure the experiment directory exists
    os.makedirs(exp_dir, exist_ok=True)
    
    out_dir = os.path.join(exp_dir, "refined_output")
    os.makedirs(out_dir, exist_ok=True)

    checkpoints_path = f"{exp_dir}/checkpoints"

    if not resume:
        checkpoint = Checkpointer(checkpoints_path)
    else:
        logger.info("Resuming the program from previous checkpoints")
        checkpoint = Checkpointer(checkpoints_path, True)

    # rename reads
    stage = "1_1"
    stage_params = [exp_dir, fastq_file]

    if checkpoint.should_run_step(stage, stage_params):
        logger.info("Renaming reads")
        rename_reads(exp_dir, fastq_file)
        
        checkpoint.log(stage, stage_params)
        logger.info("Renaming reads complete")
    else:
        logger.info("Reads already renamed")

    # obtain_read_ids
    stage = "1_2"
    stage_params = [exp_dir]

    if checkpoint.should_run_step(stage, stage_params):
        logger.info("Obtaining read ids")
        obtain_read_ids(exp_dir)
        
        checkpoint.log(stage, stage_params)
        logger.info("Obtaining read ids complete")
    else:
        logger.info("Read ids already obtained")

    # run_seq2vec
    stage = "1_3"
    stage_params = [exp_dir]

    if checkpoint.should_run_step(stage, stage_params):
        logger.info("Running seq2vec")
        run_seq2vec(exp_dir)
        
        checkpoint.log(stage, stage_params)
        logger.info("Running seq2vec complete")
    else:
        logger.info("seq2vec already ran")

    # # compile_splitreads
    # stage = "1_4"
    # stage_params = [exp_dir]

    # if checkpoint.should_run_step(stage, stage_params):
    #     logger.info("Compiling split reads")
    #     compile_splitreads(exp_dir)
        
    #     checkpoint.log(stage, stage_params)
    #     logger.info("Compiling split reads complete")
    # else:
    #     logger.info("Split reads already compiled")

    # create_overlaps
    stage = "1_4"
    stage_params = [exp_dir]

    if checkpoint.should_run_step(stage, stage_params):
        logger.info("Creating overlaps")
        create_overlaps(exp_dir)
        
        checkpoint.log(stage, stage_params)
        logger.info("Creating overlaps complete")
    else:
        logger.info("Overlaps already created")
        
        
    # run_graph_create
    stage = "2_1"
    stage_params = [exp_dir]

    if checkpoint.should_run_step(stage, stage_params):
        logger.info("Creating Graph ")
        run_create_graph(exp_dir)
        
        checkpoint.log(stage, stage_params)
        logger.info("Creating Graph complete")
    else:
        logger.info("Creating Graph executed")

    # run_step1
    stage = "2_1"
    stage_params = [exp_dir, in_file, out_dir]

    if checkpoint.should_run_step(stage, stage_params):
        logger.info("Running step 1 ")
        run_step1(exp_dir, in_file, out_dir)
        
        checkpoint.log(stage, stage_params)
        logger.info("Running step 1 complete")
    else:
        logger.info("Step1 already executed")        

    # run_step2
    stage = "2_2"
    stage_params = [exp_dir, out_dir]

    if checkpoint.should_run_step(stage, stage_params):
        logger.info("Running step 2")
        run_step2(exp_dir, out_dir)
        
        checkpoint.log(stage, stage_params)
        logger.info("Running step 2 complete")
    else:
        logger.info("Step2 already executed")

    # run_step3
    stage = "2_1"
    stage_params = [exp_dir, out_dir]

    if checkpoint.should_run_step(stage, stage_params):
        logger.info("Running step 3 ")
        run_step3(exp_dir, out_dir)
        
        checkpoint.log(stage, stage_params)
        logger.info("Running step 3 complete")
    else:
        logger.info("Step3 already executed")


    # run_seq2covvec
    stage = "4_1"
    stage_params = [exp_dir, fastq_file]

    if checkpoint.should_run_step(stage, stage_params):
        logger.info("Running seq2covvec")
        run_seq2covvec(exp_dir, fastq_file)
        
        checkpoint.log(stage, stage_params)
        logger.info("Running seq2covvec complete")
    else:
        logger.info("seq2covvec already ran")

    # run_step4
    stage = "4_2"
    stage_params = [exp_dir, out_dir, epochs]

    if checkpoint.should_run_step(stage, stage_params):
        logger.info("Running step 4")
        run_step4(exp_dir, out_dir, epochs)
        
        checkpoint.log(stage, stage_params)
        logger.info("Running step 4 complete")
    else:
        logger.info("step 4 already executed")
        
    # run_step4
    if (groundtruth != None):
      stage = "5_1"
      stage_params = [out_dir, groundtruth, fastq_file]

      if checkpoint.should_run_step(stage, stage_params):
          logger.info("Evaluating Results ... ")
          run_eval(out_dir, groundtruth, fastq_file)
        
          checkpoint.log(stage, stage_params)
          logger.info("Evaluation complete")
      else:
          logger.info("Evaluation already executed")
    else:
          logger.info("Groundtruth file not provided")

    # Todo - Check output files

    # Final message
    logger.info("Pipeline completed successfully!")
    

import os
import argparse
import logging
import time
from graphk import pipeline

def banner():
  """
  This function defines the banner text for your tool.
  """
  return """

 _____                 _     _   __      _     ______ 
|  __ \               | |   | | / /     | |    | ___ \\
| |  \/_ __ __ _ _ __ | |__ | |/ /______| |    | |_/ /
| | __| '__/ _` | '_ \| '_ \|    \______| |    |    / 
| |_\ \ | | (_| | |_) | | | | |\  \     | |____| |\ \ 
 \____/_|  \__,_| .__/|_| |_\_| \_/     \_____/\_| \_|
                | |                                   
                |_|                                   


  GraphK-LR - Long Read Binning Refiner Tool
  """

def main():
    print(banner())
    parser = argparse.ArgumentParser(description='Run the sequencing pipeline.')

    parser.add_argument('-o', '--exp_dir', type=str, help='Path to the output directory')
    parser.add_argument('-i', '--in_file', type=str, help='Path to the input file of inial binning results')
    parser.add_argument('-r', '--fastq_file', type=str, help='Path to the original reads file (fastq)')
    parser.add_argument('-e', '--epochs', type=int, default=100, help='Number of epochs (default: 100)')
    parser.add_argument('--resume', action='store_true', help='Resume the program from previous checkpoints')
    parser.add_argument('-g', '--groundtruth', type=str, help='Groundtruth file if available')

    args = parser.parse_args()

    exp_dir = args.exp_dir
    in_file = args.in_file

    # init logger
    logger = logging.getLogger('GraphKLR')
    logger.setLevel(logging.DEBUG)

    # init logger console handler - To log messages to console
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    consoleHeader = logging.StreamHandler()
    consoleHeader.setFormatter(formatter)
    consoleHeader.setLevel(logging.INFO)
    logger.addHandler(consoleHeader)

    # Ensure the experiment directory exists
    os.makedirs(exp_dir, exist_ok=True)
    
    out_dir = os.path.join(exp_dir, "refined_output")
    os.makedirs(out_dir, exist_ok=True)

    # init logger file handler - To log messages to a file
    fileHandler = logging.FileHandler(f"{exp_dir}/GraphKLR.log", mode='w')
    fileHandler.setLevel(logging.DEBUG)
    fileHandler.setFormatter(formatter)
    logger.addHandler(fileHandler)

    consoleHandler = logging.StreamHandler()

    # starting the pipeline
    start_time = time.time()

    try:
        pipeline.run_pipeline(args)
    except Exception as e:
        logger.error(f"An error occurred: {e}", exc_info=True)
    finally:
        # pipeline finished
        end_time = time.time()
        time_taken = end_time - start_time

        logger.info(f"Program Finished! Please find the output.")
        logger.info(f"Total time consumed = {time_taken:.2f} seconds")
        logger.info(f"Thank you for using GraphKLR.")

        # Remove handlers
        logger.removeHandler(fileHandler)
        logger.removeHandler(consoleHandler)
        fileHandler.close()
        consoleHandler.close()   

if __name__ == "__main__":
    main()

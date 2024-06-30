# GraphK-LR : Long-reads metagenomic bin refiner

![GitHub](https://github.com/NethmiRanasinghe/GraphK-LR)

## Dependencies
<!-- GraphK-LR is coded using C++ (v9+) and Python 3.9. To run GraphK-LR, you will need to install the following python and C++ modules. -->

A possible conda environment to work

For a system with CUDA:
```sh
conda env create -f environment.yml 
```
For a system with only CPU:
```sh
conda env create -f env_cpu.yml 
```

<!-- ### Python dependencies -->
<!-- Essential libraries

* numpy 1.16.4 
* scipy 1.3.0 
* seaborn 0.9.0
* h5py 2.9.0
* tabulate 0.8.7
* pytorch 1.4.0 -->

<!-- Essential for contig binning -->
<!-- * umap-learn -->
<!-- * fraggenescan 1.31
* hmmer 3.3.2
* HDBSCAN -->


<!-- ### C++ requirements -->
<!-- * GCC version 9.1.0 or later
* OpenMP 4.5 for multi processing -->

## Downloading the tool
To download GraphK-LR, you have to clone the GraphK-LR repository to your machine.

```
git clone https://github.com/NethmiRanasinghe/GraphK-LR
```

## Compiling the source code
* Build the binaries
```
cd GraphK-LR
```
```
sh build.sh
```    

## Test run data 
Download dataset from [here](https://anu365-my.sharepoint.com/:u:/g/personal/u6776114_anu_edu_au/EZfR_olp6fxOpoKEexcm1ZABk1fTQTnW0-3ja772k22WbA?e=aoXr7N);

```
python main.py -r <path to fastq file> -i <path to initial tool results .txt file> -o <path to output folder> --resume
```

### Parameters


### Available Commands 

Use the `-h` argument to list all the available commands.

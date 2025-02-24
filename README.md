# Tardigrades in concrete holes - data analysis
This repository contains the scripts for the data analysis associated with the paper: **PAPER NAME TO ADD**

## 1 - Create and activate the conda environment

Package management is done with conda (it needs to be installed in your system: https://docs.anaconda.com/anaconda/install/).  

You can create and install the conda environment from the yml file by typing the followinc commands in our terminal:
```
conda env create -f ./conda_envs/environment.yml
conda activate metabarcoding_env
```

If you have issues with the environment.yml file, try with the environment_detailed.yml file where every dependency version is specified:
```
conda env create -f ./conda_envs/environment_detailed.yml
conda activate metabarcoding_env
```

## 2 - Clone the repository
Clone this repository with
```
git clone https://github.com/tardipede/tardigrades_concrete_holes.git
cd microcosmos_rock_pools
```
## 3 - Download reads
Create in the working directory a folder named *data* and download there the files containing the raw reads

# SCARAB

[![Codacy Badge](https://api.codacy.com/project/badge/Grade/1a2954edef114b81a583bb23ffba2ace)](https://app.codacy.com/gh/hallamlab/scarab?utm_source=github.com&utm_medium=referral&utm_content=hallamlab/scarab&utm_campaign=Badge_Grade_Dashboard)

SAG Anchored Binner for recruiting metagenomic reads using single-cell amplified genomes as references

Check out the [wiki](https://github.com/hallamlab/scarab/wiki) for tutorials and more information on SCARAB!!

### Install SCARAB and Dependencies
Currently the easiest way to install SCARAB is to use a conda virtual environment.  
This will require the installation of conda. It can be [Anaconda](https://www.anaconda.com/download), [MambaForge](https://mamba.readthedocs.io/en/latest/installation.html), or [Miniconda](https://docs.conda.io/en/latest/miniconda.html).\
Note: that SCARAB is written in Python 3, so the conda has to support Python 3.\
**Warning!!** SCARAB is developed and tested on Linux only.

Once one of the "conda"s is installed, you can follow the directions below to install all dependencies and SCARAB within a conda environment.
```sh
git clone https://github.com/hallamlab/scarab.git
cd scarab
```
 Now use `make` to create the conda env, activate it, and install SCARAB via pip.
```sh
make install-scarabenv
conda activate scarab_cenv
make install-scarab
```

### Test SCARAB Install
Here is a small [demo dataset](https://drive.google.com/file/d/1yUoPpoNRl6-CZHkRoUYDbikBJk4yC-3V/view?usp=sharing) to make sure your SCARAB install was successful.
Just download and follow along below to run SCARAB. (make sure you've activated the SCARAB conda env)
```sh
unzip demo.zip
cd demo
scarab recruit -m k12.gold_assembly.fasta -l read_list.txt -o SCARAB_out -s SAG
```
The result of the above commands is a new directory named `SCARAB_out` that contains all the intermediate and final outputs for the SCARAB analysis. 

### Docker and Singularity containers
If you would like to use a [Docker](https://docs.docker.com/engine/install/) or [Singularity](https://docs.sylabs.io/guides/3.0/user-guide/installation.html) container of SCARAB they are available.\
In either case, they can be pulled from [Quay.IO](https://quay.io/repository/hallamlab/scarab) with one of the following commands:
```sh
# For Docker (need sudo access)
sudo docker pull quay.io/hallamlab/scarab
# Run above demo with Docker
sudo docker run -it --network=host --rm -v ./:/cwd quay.io/hallamlab/scarab:latest scarab recruit -m cwd/k12.gold_assembly.fasta -l cwd/docker_read_list.txt -o cwd/SCARAB_out -s cwd/SAG

#For Singularity
singularity pull docker://quay.io/hallamlab/scarab
# Run the above demo with Singularity
singularity exec scarab_latest.sif scarab recruit -m k12.gold_assembly.fasta -l read_list.txt -o SCARAB_out -s SAG
```
Note: Make sure you are in the `demo` directory when running the above commands as the examples assume this.

They can also be build from scratch using the following commands:\
Docker (need sudo access):
```sh
make docker-build
```
Singularity (assumes you built the Docker locally first):
```sh
make singularity-local-build
``` 

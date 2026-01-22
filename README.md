# SCARAB

SCARAB recruits metagenomic reads using single-cell amplified genomes as references.

Documentation and tutorials: [Wiki](https://github.com/RyloByte/scarab/wiki)

## Install

SCARAB is developed and tested on Linux only. A conda distribution is required (Mambaforge, Miniconda, or Anaconda).

Clone the repo:

```sh
git clone https://github.com/RyloByte/scarab.git
cd scarab
```

Create the conda environment and install SCARAB:

```sh
make install-scarabenv
conda activate scarab_cenv
make install-scarab
```

Or without make:

```sh
conda env create -f environment.yml
conda activate scarab_cenv
pip install .
```

## Test

Download the demo dataset: https://drive.google.com/file/d/1yUoPpoNRl6-CZHkRoUYDbikBJk4yC-3V/view?usp=sharing

```sh
unzip demo.zip
cd demo
scarab recruit -m k12.gold_assembly.fasta -l read_list.txt -o SCARAB_out -s SAG
```

The result is a new directory named `SCARAB_out` that contains all intermediate and final outputs.

## Containers

Container images are published to Quay.io:

```sh
# Docker (needs sudo access)
sudo docker pull quay.io/hallamlab/scarab
sudo docker run -it --network=host --rm -v ./:/cwd quay.io/hallamlab/scarab:latest \
  scarab recruit -m cwd/k12.gold_assembly.fasta -l cwd/docker_read_list.txt -o cwd/SCARAB_out -s cwd/SAG

# Singularity
singularity pull docker://quay.io/hallamlab/scarab
singularity exec scarab_latest.sif scarab recruit -m k12.gold_assembly.fasta -l read_list.txt -o SCARAB_out -s SAG
```

Build locally:

```sh
make docker-build
make singularity-local-build
```


## Developer workflow

Local conda build from the working tree (no version bump required):

```sh
make conda-build-local
make conda-test-env
conda activate scarab_test
scarab info
```

Release checklist (tagged source build):

```sh
# 1) Update version, tag, and push
# 2) Update conda-recipe/meta.yaml source URL + sha256
make conda-build-release
```

Upload to Anaconda Cloud:

```sh
ANACONDA_USER=yourname make conda-upload
```

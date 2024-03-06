#!/bin/bash

## Prepare virtual env ##
conda create -y -n flic_env python=3.10
conda activate flic_env

## Install required python and conda packages ##
pip3 install cutadapt

conda install -y -c conda-forge libstdcxx-ng  # need for install correct version of porechop
conda install -y -c conda-forge python_abi=3.10  # need for install correct version of porechop
conda install -y -c bioconda porechop

conda install -y libgcc-ng  # need for install correct version of minimap2
conda install -y -c conda-forge libzlib  # need for install correct version of minimap2
conda install -y -c bioconda minimap2=2.24

conda install -y -c bioconda samtools
conda install -y -c bioconda bedtools=2.30.0
conda install -y -c bioconda ucsc-bedgraphtobigwig

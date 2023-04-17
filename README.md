# TERA-Seq analysis and methods
This repository contains all the main tools and scripts needed to reproduce analysis used in publication **TERA-Seq: True end-to-end sequencing of native RNA molecules for transcriptome characterization** (Ibrahim, F.\*, Oppelt, J.\*, Maragkakis, M.\* and Mourelatos, Z., 2021, Nucleic Acids Research, DOI: [10.1093/nar/gkab713](https://doi.org/10.1093/nar/gkab713), PMID: [34428294](https://pubmed.ncbi.nlm.nih.gov/34428294/). \* The authors wish it to be known that, in their opinion, the first three authors should be regarded as Joint First Authors. If you use any part of the analysis, code or samples we kindly ask you to cite this publication.

This repository will be archived as soon as the publication is accepted to keep the code in the same *version* as in the time of publication.

<font size="10">For **primary TERA-Seq data processing** we strongly recommend to use [**TERA-Seq Snakemake workflow**](https://github.com/opplatek/TERA-Seq_snakemake).</font>

## Samples
**fast5** (basecalled, untrimmed) and **fastq** files can be found in the [`samples`](samples) directory.

## General information

The analysis was tested on Ubuntu 16.04 LTS, 16 threads, 32 GB RAM, 2 TB HDD. The only exception is `STAR` indexing and mapping which needs ~32GB RAM for human-sized genome. If you need to run `STAR` on a machine with lower amount of RAM, you can add `--genomeSAsparseD 2` (or higher value) to `STAR` genome indexings to decrease RAM requirements for at a cost of speed reduction.

## Requirements

### Docker
Using Docker is strongly recommended. The Docker container that could be used to reanalyze the data is available here:

    docker pull joppelt/teraseq:latest

If you don't have access do Docker, please follow the instructions for [manual installation below](#manual-installation).

### Manual installation
There are seven system requirements: **`git`**, **`gcc`**, **`make`**, **`wget`**, **`g++`**, **`zlib-dev`**, **`bzip2`**. `g++`, `zlib-dev`, and `bzip2` are Nanopolish requirements. If you don't plan to run Nanopolish, you can skip these three. You can install them (on Ubuntu) with:

    apt-get update
    apt-get install -y git gcc make wget g++ zlib1g-dev bzip2
    # rm -rf /var/lib/apt/lists/*

Make a local copy of this repository anywhere you want.

    git clone https://github.com/mourelatos-lab/TERA-Seq_manuscript

Then, we would like to specificy some system-wide variables we will use in this analysis. The variables used throughout the analysis can be found in `PARAMS.sh` in the main `TERA-Seq_manuscript` directory you just downloaded. You can open in it and edit the variables so it suits to your needs. We will assume you will work on the analysis in the `TERA-Seq_manuscript` directory.

The individual analysis scripts **will assume the structure** from `PARAMS.sh`.

#### Conda environment
Most of the software used for the analyses can be installed using the provided Conda environment yml (or text) file. Some additional software had to be installed outside  Conda environment for various reason.

To get the Conda environment, you need to get Conda environment manager ([Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)) first. We recommend Miniconda but Conda installation should work as well. The tested version is [Miniconda3-py37_23.1.0-1-Linux-x86_64](https://repo.anaconda.com/miniconda/Miniconda3-py37_23.1.0-1-Linux-x86_64.sh).

Once you installed Conda manager (tested with `conda 4.9.0`), we recommend to install [Mamba](https://mamba.readthedocs.io/en/latest/index.html) to significantly speed up the installation (90% of the time, it works every time):

    conda install mamba -n base -c conda-forge 


##### From provided YAML file (recommended)
Either w/ or w/o you can processed to install the main environment. We recommend to use the included yaml file `teraseq.yml`. Choose `mamba` or default `conda` installer, not both:

    cd TERA-Seq_manuscript/
    # mamba env create -f teraseq-env.yml # If you installed Mamba
    conda env create -f teraseq-env.yml # If you use base Conda installer

This will install the `teraseq` environment. In case yaml installation doesn't work, you can try to install the environment using the provided `teraseq-env.txt` using the same command. If you need to install the dependencies manually, please follow instructions in the [manual Conda environment install](#manual-conda-environment-install).

##### Manual Conda environment install
In case you need a manual installation, you can also use (substitute `conda` for `mamba` if you installed it, it is going to be much faster):

    # Intiate environment
    conda create --name teraseq python=3.7
    conda activate teraseq
    # Set channels priority
    conda config --add channels 'conda-forge'
    conda config --add channels 'r'
    conda config --add channels 'bioconda'
    conda config --add channels 'defaults'
    # Install software packages
    conda install -c bioconda \
    r=3.5.1 r-base \
    r-optparse=1.6.2 r-dplyr r-rio r-ggplot2 r-bit64 r-reshape2 r-devtools r-ggpubr \
    r-ggridges r-corrplot r-rsqlite r-futile.logger r-plyr r-readr r-scales r-viridis \
    r-data.table r-yaml r-treemap r-venndiagram r-extrafont r-hdf5r r-patchwork r-r.utils \
    bioconductor-ensembldb bioconductor-tximport bioconductor-dupradar bioconductor-gviz bioconductor-deseq2 \
    fastqc STAR=2.7.2b samtools">=1.10" bedtools=2.29.0 samblaster trimmomatic=0.39 \
    perl-cpan-shell perl-app-cpanminus perl-local-lib \
    umi_tools fastx_toolkit bbmap=38.67 minimap2=2.17 subread=2.0.0 pigz \
    seqkit seqtk picard=2.20.8 sambamba=0.7.1 gmap=2019.09.12 gffread \
    ucsc-bedgraphtobigwig ucsc-bedtobigbed ucsc-genePredToBed ucsc-gtfToGenePred ucsc-liftover ucsc-bigbedtobed \
    qualimap rseqc cd-hit bedops parallel

Note: If you get a lot of `ClobberError: This transaction has incompatible packages due to a shared path.` errors try to restart your Conda environment and/or exectue `conda clean --all` to [clean the Conda cache](https://github.com/conda/conda/issues/7038).

#### Additional software
`tools` directory contains information additional software required for some of the analyses. Installation of the additioonal software **is included in the main** `run.sh`. If you wish to install it manually you can execute `run.sh` in the `tools` directory. We use virtual environments wherever possible.

**Important: If your Conda main path is not exported by default, you have to specify `$CONDA_PREFIX` variable in `PARAMS.sh` file.**. If you **cannot** use Conda at this stage, you will need system-wide install **curl**, **Perl 5.22**, Perl **local::lib**, and **R 3.5.1**.

    cd tools/
    ./run.sh
    cd ../

#### Guppy basecaller
Guppy basecaller must be obtained from [Nanopore community](https://community.nanoporetech.com/downloads) section and must be installed in case you want to re-basecall the data. Please use **Guppy 3.2.2 GPU** version to replicate the basecall. The CPU version should provide the same results. Installation manual is provided in the same link.

### Miscellaneous scripts
`misc` directory contains various miscellaneous scripts and sources.

## Sample data
To re-analyze the data **please populate** the `samples` directory **first**. For more information see [`samples`](samples) directory [`README`](samples/README.md).

Note:  5TERA 5' adapter might and will be called *REL5/rel5*, and TERA3 3' end adapter might and will be called *REL3/rel3* throughout the analyses.

## Data and Analysis
The main `run.sh` in the cloned `TERA-Seq_manuscript` directory should **execute all the analyses** used in the manuscript including the installation of additional software and preparing the references.

**Reminder**: The only two exceptions are the **installation of `teraseq` Conda** environment and **populating the samples directory** which you have to do **before** you can run the analyses.

If you wish to run only selected analyses you can go to the individual directories and execute the *local* `run.sh`. Simply go to the selected analysis and run the `run.sh` script. For example:

    cd metacoord_correction/
    ./run.sh

Each analysis directory contains main `run.sh` script which shows how to reproduce the analysis and visualization if applicable. Visuzalization and other additional scripts (for example R visualization) can be found in `src` directory.  In most cases, executing the `run.sh` should be enough to reproduce the analysis. However, if this fails, the content of `run.sh` script should be comment enough for you to see where the error comes from. The `run.sh` also contains information about how to to execute the individual scripts in the `src` directory. `dev` directory contains additional used tools.

Most of the analyses has been created by [Jan Oppelt](mailto:jan.oppelt@pennmedicine.upenn.edu) ([Mourelatos lab](http://mourelatos.med.upenn.edu/), Department of Pathology and Laboratory Medicine,  Perelman School of Medicine, University of Pennsylvania, USA). Most of the `dev` tools have been developed by [Emmanouil "Manolis" Maragkakis](mailto:emmanouil.maragkakis@nih.gov) (currently [Laboratory of Genetics and Genomics, National Institute on Aging](https://www.nia.nih.gov/research/labs/lgg/computational-genomics-unit), NIH, USA) or [Panagiotis Alexiou](mailto:panagiotis.alexiou@ceitec.muni.cz) (currently [RBP Bioinformatics](https://www.ceitec.eu/rbp-bioinformatics-panagiotis-alexiou/rg281) & [Bioinformatics Core Facility](https://www.ceitec.eu/bioinformatics-core-facility/cf284), CEITEC, CZE).

### Data and preprocessing
#### Data preparation
The `data` directory `run.sh` scripts will download and format all the references and annotations required for the following analyses. It also downloads the additional external data.

Note: Pre-compiled references are temporarily also available to download [here](TODO).
#### Preprocessing
The `samples` directory contains several run scripts additonal to the main `run.sh`. This is because different library types are processed differently. However, it should be possible to run all the processing simply by executing the `run.sh`. The Guppy basecalling example is done only for one of the libraries. Reminder: the included data are only **examples** and not the complete libraries.
`run_guppy.sh` script contains commands used to basecall the libraries (`Guppy 3.3.2-1`). `run_5TERA.sh`, `run_TERA3`, `run_5TERA3.sh`, `run_5TERA-short.sh`, `run_5TERA-SIRV.sh`, `run_Akron5Seq.sh`, and `run_RiboSeq.sh` scripts contain steps to preprocess and align the reads and postproces the resulting bam files. Script `run_5TERA3-merge.sh` contains instructions to merge the 5TERA3 replicates to a single sample. It also shows how to create the `sqlite3` database ([CLIPSeqTools](http://mourelatos.med.upenn.edu/clipseqtools/)) necessary for some of the analyses.

### Analysis
#### Alignment statistics
The `align-stats` directory contains code to reproduce analysis summarized in: **Supplementary Table 2**; **Supplementary Figure 1d**.
#### Adapter length
The `adapter` directory contains code and scripts to reproduce analysis summarized in: **Supplementary Figure 1b**.
#### Meta-coordinates, re-annotation and heatmap
The `metacoord_correction` directory contains code and scripts to reproduce analysis summarized in: **Figure 1b,c,d** (left part of Figures 1b,c); **Figure 3b**; **Figure 5c**; **Supplementary Figure 2a,b,c,d**.
#### Changes after re-annotation
The `reannot-change` directory contains code and scripts to reproduce analysis summarized in: **Figure 1b,c** (right part of the figures).
#### Transcript coverage (transcriptome and genome)
The `trans-coverage` directory contains code and scripts to reproduce analysis summarized in: **Figure 2b,c**; **Figure 3c**; **Figure 4b**; **Figure 5d**; **Supplementary Figure 5**.
#### Poly(A) tail length
The `polya` directory contains code and scripts to reproduce analysis summarized in: **Figure 4a**; **Figure 5b**; **Supplementary Figure 4b,c,d,e**.
#### CAGE and APA
The `cage_apa` directory contains code and scripts to reproduce analysis summarized in: **Figure 2a (CAGE)**; **Figure 3d (APA)**; **Supplementary Figure 3a,b,c (CAGE)**, **Supplementary Figure 4a (APA)**.
#### Promoter signal
The `promoter-heatmap` directory contains code and scripts to reproduce analysis summazired in: **Figure 2d**.
#### Relative position distribution
The `relative-pos-distro` directory contains code and scripts to reproduce analysis summarized in: **Figure 6b,c,d**.
#### Conservation around cuts
The `conservation` directory contains code and scripts to reproduce analysis summarized in: **Figure 6e**.
#### SIRV expression distribution and alignment
The `sirv` directory contains code and scripts to reproduce analysis summarized in: **Supplementary Figure 1c**.
#### Meta-coordinates and poly(A) tail length correlation
The `metacoord-vs-polya` directory contains code and scripts to reproduce analysis summarized in: **Figure 4c,d**.
#### Expression correlation
The `expression` directory contains code and scripts to reproduce analysis summarized in: **Supplementary Figure 6**.

> Written with [StackEdit](https://stackedit.io/).

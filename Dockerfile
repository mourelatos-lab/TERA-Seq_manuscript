# Download base image ubuntu 16.04
FROM ubuntu:16.04

# LABEL about the custom image
LABEL maintainer="jan.oppelt@pennmedicine.upenn.edu"
LABEL version="0.1"
LABEL description="This is custom Docker Image for \
analysis of TERA-Seq publication (DOI: https://doi.org/10.1093/nar/gkab713)."

# Disable Prompt During Packages Installation
ARG DEBIAN_FRONTEND=noninteractive

# Set default shell
SHELL ["/bin/bash", "-c"]

### System-wide requirements; cpanminus is not required if Perl uses virtual environment method; g++ and zlib1g-dev are required only for Nanopolish
RUN apt-get update && apt-get install -y git gcc make wget g++ zlib1g-dev cpanminus && rm -rf /var/lib/apt/lists/*

### Main GitHub repo
WORKDIR /root
RUN git clone https://github.com/mourelatos-lab/TERA-Seq_manuscript.git

### Install Miniconda3
ENV PATH "/root/miniconda3/bin:${PATH}"
ARG PATH="/root/miniconda3/bin:${PATH}"
RUN echo ${PATH}

RUN wget \
    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && mkdir /root/.conda \
    && bash Miniconda3-latest-Linux-x86_64.sh -b \
    && rm -f Miniconda3-latest-Linux-x86_64.sh
RUN conda --version

# Install Mamba for faster installation
RUN conda install -c conda-forge mamba

# Get Conda yml and install environment
RUN mamba env create -f /root/TERA-Seq_manuscript/teraseq-env.yml

# Increase default FastQC RAM
RUN sed -i 's/-Xmx250m/-Xmx5g/g' /root/miniconda3/envs/teraseq/opt/fastqc-*/fastqc

#ENV PATH="${PATH}:/root/miniconda3/envs/teraseq/bin"

RUN ln -s /root/miniconda3/envs/teraseq/bin/R /bin/R \
    && ln -s /root/miniconda3/envs/teraseq/bin/curl /bin/curl

### Save default Conda path
RUN sed -i '/CONDA_PREFIX/d' /root/TERA-Seq_manuscript/PARAMS.sh \
    && echo -e "CONDA_PREFIX=\"/root/miniconda3\"" >> /root/TERA-Seq_manuscript/PARAMS.sh

### Perl
WORKDIR /root/TERA-Seq_manuscript/tools

## System-wide install (option 1)
# RUN cpanm inc::Module::Install \
#     && cpanm autodie \
#     && cpanm DBI \
#     && cpanm Devel::Size \
#     && cpanm Getopt::Long::Descriptive \
#     && cpanm IO::File \
#     && cpanm IO::Interactive \
#     && cpanm IO::Uncompress::Gunzip \
#     && cpanm Params::Validate \
#     && cpanm Params::Util \
#     && cpanm Sub::Install \
#     && cpanm Modern::Perl \
#     && cpanm --force MooseX::App::Simple \
#     && cpanm --force MooseX::App::Command \
#     && cpanm --force MooseX::Getopt::Meta::Attribute::Trait::NoGetopt
#
# WORKDIR /usr/local/share/perl/5.22.1
#
# RUN git clone --recursive https://github.com/genoo/GenOO.git GenOO_git
# WORKDIR /usr/local/share/perl/5.22.1/GenOO_git
# RUN git reset 6527029 --hard
# WORKDIR /usr/local/share/perl/5.22.1
# RUN mkdir GenOO \
#     && cp -r GenOO_git/lib/GenOO/* GenOO/
#
# RUN cpanm CLIPSeqTools
#
# RUN wget https://raw.githubusercontent.com/mourelatos-lab/TERA-Seq_manuscript/main/misc/GenOOx/Data/File/SAMminimap2.pm -O /usr/local/share/perl/5.22.1/GenOOx/Data/File/SAMminimap2.pm
# RUN mkdir /usr/local/share/perl/5.22.1/GenOOx/Data/File/SAMminimap2 \
#     && wget https://raw.githubusercontent.com/mourelatos-lab/TERA-Seq_manuscript/main/misc/GenOOx/Data/File/SAMminimap2/Record.pm -O /usr/local/share/perl/5.22.1/GenOOx/Data/File/SAMminimap2/Record.pm

## Virtual environment install (option 2)
RUN git clone https://github.com/jizhang/perl-virtualenv.git \
    && cd perl-virtualenv/ \
    && git reset f931774 --hard \
    && chmod u+x virtualenv.pl \
    && ./virtualenv.pl teraseq \
    && . teraseq/bin/activate \
    && curl -L https://cpanmin.us/ -o teraseq/bin/cpanm \
    && chmod +x teraseq/bin/cpanm

RUN . perl-virtualenv/teraseq/bin/activate \
    && cpanm inc::Module::Install \
    && cpanm autodie \
    && cpanm DBI \
    && cpanm Devel::Size \
    && cpanm Getopt::Long::Descriptive \
    && cpanm IO::File \
    && cpanm IO::Interactive \
    && cpanm IO::Uncompress::Gunzip \
    && cpanm Params::Validate \
    && cpanm Params::Util \
    && cpanm Sub::Install \
    && cpanm Modern::Perl \
    && cpanm --force MooseX::App::Simple \
    && cpanm --force MooseX::App::Command \
    && cpanm --force MooseX::Getopt::Meta::Attribute::Trait::NoGetopt

RUN git clone --recursive https://github.com/genoo/GenOO.git perl-virtualenv/teraseq/lib/perl5/GenOO_git \
    && cd perl-virtualenv/teraseq/lib/perl5/GenOO_git/ \
    && git reset 6527029 --hard \
    && cd ../ \
    && mkdir GenOO \
    && cp -r GenOO_git/lib/GenOO/* GenOO/

RUN . perl-virtualenv/teraseq/bin/activate \
    && cpanm CLIPSeqTools \
    && cp -r /root/TERA-Seq_manuscript/misc/GenOOx/* perl-virtualenv/teraseq/lib/perl5/GenOOx/

### Nanopolish
# Default version
RUN git clone --recursive https://github.com/jts/nanopolish.git \
    && mv nanopolish nanopolish-480fc85 \
    && cd nanopolish-480fc85/ \
    && git reset 480fc85 --hard \
    && sed -i 's#http://bitbucket.org/eigen/eigen/get/$(EIGEN_VERSION).tar.bz2#https://gitlab.com/libeigen/eigen/-/archive/$(EIGEN_VERSION)/eigen-$(EIGEN_VERSION).tar.bz2#' Makefile \
    && sed -i 's/tar -xjf $(EIGEN_VERSION).tar.bz2/tar -xjf eigen-$(EIGEN_VERSION).tar.bz2/' Makefile \
    && sed -i 's/eigen-eigen-\*/eigen-$(EIGEN_VERSION)/' Makefile \
    && rm -rf fast5 \
    && git clone https://github.com/mateidavid/fast5.git \
    && cd fast5/ \
    && git reset 18d6e34 --hard \
    && cd ../ \
    && rm -rf htslib \
    && git clone --recursive https://github.com/samtools/htslib.git \
    && cd htslib/ \
    && git reset 3dc96c5 --hard \
    && cd ../ \
    && make \
    && ln -s $(pwd)/nanopolish /root/miniconda3/envs/teraseq/bin/nanopolish

# New version with polya hmm scripts
RUN git clone --recursive https://github.com/jts/nanopolish.git \
    && mv nanopolish nanopolish-ab9722b \
    && cd nanopolish-ab9722b/ \
    && git reset ab9722b --hard

### Other dependencies
# Make sure to activate Conda
SHELL ["conda", "run", "-n", "teraseq", "/bin/bash", "-c"]

## GeneCycle
RUN Rscript -e 'install.packages("GeneCycle", repos="https://cloud.r-project.org")'

## Cutadapt
RUN mkdir cutadapt-2.5 \
    && cd cutadapt-2.5/ \
    && python3 -m venv venv \
    && source venv/bin/activate \
    && python3 -m pip install --upgrade pip \
    && pip3 install cutadapt==2.5 pysam numpy pandas matplotlib seaborn \
    && which cutadapt

## DeepTools
RUN mkdir deepTools-3.5.0 \
    && cd deepTools-3.5.0/ \
    && python3 -m venv venv \
    && source venv/bin/activate \
    && python3 -m pip install --upgrade pip \
    && pip3 install wheel \
    && pip3 install deeptools==3.5.0 \
    && deeptools --version

## ONT-Fast5-API
RUN mkdir ont-fast5-api \
    && cd ont-fast5-api/ \
    && python3 -m venv venv \
    && source venv/bin/activate \
    && pip install ont-fast5-api==3.3.0 h5py seaborn

## Jvarkit
RUN git clone "https://github.com/lindenb/jvarkit.git" \
    && cd jvarkit/ \
    && ./gradlew biostar84452 \
    && mkdir $CONDA_PREFIX/share/jvarkit \
    && ln -s /root/jvarkit/dist/biostar84452.jar /root/miniconda3/envs/teraseq/share/jvarkit/remove-softlip.jar

WORKDIR /root/TERA-Seq_manuscript
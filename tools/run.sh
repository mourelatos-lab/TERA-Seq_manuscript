#!/bin/bash
#
# Install additional software not included in the Conda environment
#

conda activate teraseq

source ../PARAMS.sh

echo ">>> INSTALL GeneCycle R PACKAGE <<<"
# This R-package is not available from Conda so we have to install it manually
# Installing packages manually in Conda environment is NOT recommended
Rscript -e 'install.packages("GeneCycle", repos="https://cloud.r-project.org")'

echo ">>> INSTALL CUTADAPT <<<"

cd $INSTALL/
mkdir cutadapt-2.5/
cd cutadapt-2.5/
python3 -m venv venv # Make Python virtual environment
source venv/bin/activate
pip3 install cutadapt==2.5
which cutadapt # Check installation
cutadapt --version
deactivate

echo ">>> INSTALL DEEPTOOLS <<<"
cd $INSTALL/
mkdir deepTools-3.5.0
cd deepTools-3.5.0/
python3 -m venv venv
source venv/bin/activate
pip3 install deeptools==3.5.0
# pip3 install -Iv deeptools==3.5.0
deeptools --version
deactivate

echo ">>> INSTALL ONT-FAST5-API <<<"
cd $INSTALL/
mkdir ont-fast5-api
cd ont-fast5-api/
python3 -m venv venv
source venv/bin/activate
pip install ont-fast5-api==3.3.0 h5py seaborn
deactivate

echo ">>> INSTALL JVARKIT <<<"

cd $INSTALL/
git clone "https://github.com/lindenb/jvarkit.git" # commit "ebcbaba"
cd jvarkit/
# Note: there has been substantial changes in the compilation recipe so we cannot reuse the 'ebcbaba' commit version; we provide the compiled jar in tools/utils but we recommend to compile your own
#git reset ebcbaba --hard
./gradlew biostar84452
mkdir $CONDA_PREFIX/share/jvarkit # the git commit version
ln -s $INSTALL/jvarkit/dist/biostar84452.jar $CONDA_PREFIX/share/jvarkit/remove-softlip.jar
#ln -s $INSTALL/tools/utils/biostar84452.jar $CONDA_PREFIX/share/jvarkit/remove-softlip.jar # Use this to link the used commit instead of the recent one

CONDA_BIN=$CONDA_PREFIX/bin
conda deactivate

echo ">>> INSTALL NANOPOLISH <<<"
# version used for polya tail estimates
cd $INSTALL/
git clone --recursive https://github.com/jts/nanopolish.git
mv nanopolish nanopolish-480fc85
cd nanopolish-480fc85/
git reset 480fc85 --hard
# fix outdated link to eigen and some code
sed -i 's#http://bitbucket.org/eigen/eigen/get/$(EIGEN_VERSION).tar.bz2#https://gitlab.com/libeigen/eigen/-/archive/$(EIGEN_VERSION)/eigen-$(EIGEN_VERSION).tar.bz2#' Makefile
sed -i 's/tar -xjf $(EIGEN_VERSION).tar.bz2/tar -xjf eigen-$(EIGEN_VERSION).tar.bz2/' Makefile
sed -i 's/eigen-eigen-\*/eigen-$(EIGEN_VERSION)/' Makefile
make
ln -s $(pwd)/nanopolish $CONDA_BIN/nanopolish

# new version with polya hmm scripts
cd $INSTALL/
git clone --recursive https://github.com/jts/nanopolish.git
mv nanopolish nanopolish-ab9722b
cd nanopolish-ab9722b/
git reset ab9722b --hard

echo ">>> INSTALL PERL - ENVIRONMENT <<<"
echo "Note: This might take some time"
# Note: All Perl libraries might be possible to install to Conda but from compatibility issues we have install them separately. Also, Perl likes to break Conda environments so it's safer to make them separate.

cd $INSTALL/
git clone https://github.com/jizhang/perl-virtualenv.git # commit f931774
cd perl-virtualenv/
git reset f931774 --hard
chmod u+x virtualenv.pl
./virtualenv.pl teraseq
source teraseq/bin/activate
# Make sure you have cpanm working - redownload
curl -L https://cpanmin.us/ -o teraseq/bin/cpanm
chmod +x teraseq/bin/cpanm
which perl && perl -v
which cpanm && cpanm -v

echo ">> INSTALL PERL - MODULES <<"
# Note: If installation of a module fails, try to rerun the installation with `--force`.

cpanm inc::Module::Install
cpanm autodie
cpanm DBI
cpanm Devel::Size
cpanm Getopt::Long::Descriptive
cpanm IO::File
cpanm IO::Interactive
cpanm IO::Uncompress::Gunzip
cpanm Params::Validate
cpanm Params::Util
cpanm Sub::Install
cpanm Modern::Perl
cpanm MooseX::App::Simple
cpanm MooseX::Getopt::Meta::Attribute::Trait::NoGetopt

echo ">> INSTALL PERL - GENOO <<"
# Note: Please use the GitHub version as the CPAN version is not up-to-date:
git clone --recursive https://github.com/genoo/GenOO.git teraseq/lib/perl5/GenOO_git # commit 6527029
cd teraseq/lib/perl5/GenOO_git/
git reset 6527029 --hard
cd ../
mkdir GenOO
cp -r GenOO_git/lib/GenOO/* GenOO/

echo ">> INSTALL PERL - CLIPSeqTools <<"
# Install CLIPSeqTools
cpanm CLIPSeqTools

echo ">> INSTALL PERL - GenOOx minimap2 parser <<"
cp -r $DIR/misc/GenOOx/* $INSTALL/perl-virtualenv/teraseq/lib/perl5/GenOOx/

echo ">>> ALL DONE <<<"

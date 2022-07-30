
Main scripts to setup and run MAJIQ

#-----------
# BACKGROUND
#-----------

# DOCS
# * Tutorial PDF: https://majiq.biociphers.org/docs/TUTORIAL%20-%20V1.5.pdf

# MAJIQ is a resource for quantifying local splice variations (LSVs) from RNA-seq data. MAJIQ has 2 components:
# 
# 1. MAJIQ builder: Use RNA-seq BAMs and GTF/GFF annotations to define splice graphs and known/novel LSVs.
# .. note: majiq doesn't like bash variables; full paths need to be specified
#
# 2. MAJQI quantifier: Use builder outputs to quantify LSV PSI and delta PSI between conditions.
# 
# USING THE BUILDER
# 
# First, prep a config file for your aligned experiments. Provided documentation shows an example config file and builder call.
# 
# For our purposes, we define the following congig file in simple text .INI format.

#-------------------
# MAJIQ SETUP CONDA
#-------------------

# make the new env
conda create -n majiq python=3.8.0; conda activate majiq
# install dependencies with conda
conda config --add channels r
conda config --add channels bioconda
conda install pysam cython h5py

# to download current majiq version
export HTSLIB_LIBRARY_DIR=/path/to/htslib/lib
export HTSLIB_INCLUDE_DIR=/path/to/htslib/include
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/usr/local/lib"
pip install git+https://bitbucket.org/biociphers/majiq_academic.git

#---------------------------
# FILES AND LOCATIONS -- HX1
#---------------------------
# Make the output dir
base_path=/home/metamaden/ri_results/gb_revisions/majiq/SRR2911306/
sudo mkdir $base_path
sudo chmod 777 $base_path

# Make and save the config file
ini_path=$base_path'settings.ini'
sudo touch $ini_path
vi $ini_path

# [info]
# bamdirs=/eternity/data/RI_benchmarking_hx1/
# genome=GRCh38
# strandness=None
# [experiments]
# SRR2911306=SRR2911306.sorted

# Make the console redirect output file
touch $base_path'out'
chmod 777 $base_path'out'

# Run MAJIQ
# Start new screen session
screen -S majiq_run_hx1

# activate the virtual env
conda activate majiq

# run majiq build
majiq build \
	/eternity/data/RI_benchmarking_resources/gencode_v35_annotation_files/gencode.v35.primary_assembly.annotation.gff3 \
	-o /home/metamaden/ri_results/gb_revisions/majiq/SRR2911306/ \
	-c /home/metamaden/ri_results/gb_revisions/majiq/SRR2911306/settings.ini \
	-j 10 2> /home/metamaden/ri_results/gb_revisions/majiq/SRR2911306/out

# Quantify PSI from MAJIQ table
majiq psi /home/metamaden/ri_results/gb_revisions/majiq/SRR2911306/SRR2911306.sorted.majiq \
-j 10 -o /home/metamaden/ri_results/gb_revisions/majiq/SRR2911306/psi -n hx1 \
2> /home/metamaden/ri_results/gb_revisions/majiq/SRR2911306/out_psi
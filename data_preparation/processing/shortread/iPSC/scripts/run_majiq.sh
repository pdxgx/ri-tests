
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

#----------------------------
# FILES AND LOCATIONS -- iPSC
#----------------------------
The BAM file is available in full or separated by chr.

# path to full bam -- iPSC
/eternity/data/RI_benchmarking_BAMs/SRR6026510.sorted.bam 

# ipsc config file
#
# ```
# [info]
# bamdirs=/eternity/data/RI_benchmarking_BAMs/
# genome=GRCh38
# strandness=None
# [experiments]
# SRR6026510=SRR6026510.sorted
# ```

# Make and store the ini file with:
#
# ```
ini_path=/home/metamaden/ri_results/gb_revisions/majiq/SRR6026510/settings.ini
touch $ini_path
vi $ini_path
# ```

# Make the output dir
#
sudo mkdir /home/metamaden/ri_results/gb_revisions/majiq/SRR6026510/
sudo chmod 777 /home/metamaden/ri_results/gb_revisions/majiq/SRR6026510/

# Make the console redirect output file
touch /home/metamaden/ri_results/gb_revisions/majiq/SRR6026510/out
chmod 777 /home/metamaden/ri_results/gb_revisions/majiq/SRR6026510/out

# The command to run the builder is as follows:
# enter new screen session
screen -S majiq_run

# activate the virtual env
conda activate majiq

# define the base path for files
base_path='/home/metamaden/ri_results/gb_revisions/majiq/SRR6026510/'
anno_path='/eternity/data/RI_benchmarking_resources/gencode_v35_annotation_files/gencode.v35.primary_assembly.annotation.gff3'

# do majiq build
# note: majiq doesn't like bash variables; full paths need to be specified
# run build
# majiq build $anno_path -o $base_path -c $ini_path 2> $base_path'/out'
# run build debug
# majiq build $anno_path -o $base_path -c $base_path'/settings.ini' --debug

# run majiq build
majiq build \
	/eternity/data/RI_benchmarking_resources/gencode_v35_annotation_files/gencode.v35.primary_assembly.annotation.gff3 \
	-o /home/metamaden/ri_results/gb_revisions/majiq/SRR6026510/ \
	-c /home/metamaden/ri_results/gb_revisions/majiq/SRR6026510/settings.ini \
	-j 10 2> /home/metamaden/ri_results/gb_revisions/majiq/SRR6026510/out

# Output in `out` file looks like:
# 2022-06-19 15:24:52,948 (PID:18201) - INFO - Majiq Build v2.4.dev3+g85d0781
# 2022-06-19 15:24:52,948 (PID:18201) - INFO - Command: /home/metamaden/anaconda3/envs/majiq/bin/majiq build /eternity/data/RI_benchmarking_resources/gencode_v35_annotation_files/gencode.v35.primary_assembly.annotation.gff3 -o /home/metamaden/ri_results/gb_revisions/majiq/SRR6026510/ -c /home/metamaden/ri_results/gb_revisions/majiq/SRR6026510/settings.ini -j 10
# 2022-06-19 15:24:52,948 (PID:18201) - INFO - Parsing GFF3
# ...
# ...
# 2022-06-19 15:31:35,819 (PID:18201) - INFO - SRR6026510.sorted: 89380 LSVs
# 2022-06-19 15:31:36,469 (PID:18201) - INFO - MAJIQ Builder is ended successfully!

# The resulting files are:
# * .sql -- the splicegraph db that can be used with VIOLA
# * .majiq.log -- logging output from the build
# * . sj -- intermediate files that can be used instead of BAMs to speed up future runs
# * .majiq -- inputs for majiq quantifiers

# Quantify PSI from MAJIQ table
majiq psi /home/metamaden/ri_results/gb_revisions/majiq/SRR6026510/SRR6026510.sorted.majiq \
-j 10 -o /home/metamaden/ri_results/gb_revisions/majiq/SRR6026510/psi -n ipsc \
2> /home/metamaden/ri_results/gb_revisions/majiq/SRR6026510/out_psi

# The `out_psi` output looks like:
# 2022-06-19 16:08:08,316 (PID:25856) - INFO - Majiq psi v2.4.dev3+g85d0781
# 2022-06-19 16:08:08,316 (PID:25856) - INFO - Command: /home/metamaden/anaconda3/envs/majiq/bin/majiq psi /home/metamaden/ri_results/gb_revisions/majiq/SRR6026510/SRR6026510.sorted.majiq -j 10 -o /home/metamaden/ri_results/gb_revisions/majiq/SRR6026510/psi -n ipsc
# 2022-06-19 16:08:08,316 (PID:25856) - INFO - Running Psi ...
# 2022-06-19 16:08:08,316 (PID:25856) - INFO - GROUP: ['/home/metamaden/ri_results/gb_revisions/majiq/SRR6026510/SRR6026510.sorted.majiq']
# 2022-06-19 16:08:08,316 (PID:25856) - INFO - Parsing file: /home/metamaden/ri_results/gb_revisions/majiq/SRR6026510/SRR6026510.sorted.majiq
# 2022-06-19 16:08:11,050 (PID:25856) - INFO - Group ipsc: 78350 LSVs
# 2022-06-19 16:08:32,627 (PID:25856) - INFO - Computation done, saving results....
# 2022-06-19 16:10:31,620 (PID:25856) - INFO - PSI calculation for ipsc ended successfully! Result can be found at /home/metamaden/ri_results/gb_revisions/majiq/SRR6026510/psi

#-----------------------------------
# MAJIQ COMMAND BUILDER TOOL -- IPSC
#-----------------------------------
# This is an example from the command builder tool located at:
# <https://biociphers.bitbucket.io/majiq-docs-academic/commandbuilder.html>

# ini file
# 
# [info]
# bamdirs=/eternity/data/RI_benchmarking_BAMs/
# genome=hg38
# 
# [experiments]
# ipsc=ipsc
# ```

# build command
majiq build /eternity/data/RI_benchmarking_resources/gencode_v35_annotation_files/gencode.v35.primary_assembly.annotation.gff3 -o /home/metamaden/ri_results/gb_revisions/majiq/SRR6026510//build -c /home/metamaden/ri_results/gb_revisions/majiq/SRR6026510/settings.ini --min-experiments 0.00 -j 10 

# quantify command
majiq psi -o /home/metamaden/ri_results/gb_revisions/majiq/SRR6026510/psi -j 10 -n ipsc_ipsc /home/metamaden/ri_results/gb_revisions/majiq/SRR6026510/build/ipsc.majiq

# voila command
voila tsv -j 10 -f /home/metamaden/ri_results/gb_revisions/majiq/SRR6026510//voila_output.tsv /home/metamaden/ri_results/gb_revisions/majiq/SRR6026510/
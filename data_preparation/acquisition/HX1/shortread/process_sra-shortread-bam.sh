#!/usr/bin/env sh

# make the conda env
# conda activate sra_process
# conda install -c bioconda/label/cf201901 sra-tools

# download sratoolkit
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.9.6/sratoolkit.2.9.6-centos_linux64.tar.gz
tar xzvf sratoolkit.2.9.6-centos_linux64.tar.gz
export PATH=$PATH:/home/metamaden/sratoolkit.2.9.6-centos_linux64/bin 

# download the hx1 short read sample
# srrid=SRR2911306
# ./home/metamaden/sratoolkit.2.4.1-ubuntu64/bin/fastq-dump.2.4.1 --help 
# fastq-dump --accession SRR2911306 --outdir ./eternity/data/RI_benchmarking_hx1/ --gzip
sudo ~/sratoolkit.2.9.6-centos_linux64/bin/fastq-dump --accession SRR2911306 \
    --outdir ./eternity/data/RI_benchmarking_hx1/ --gzip
# returns:
# > Read 24463210 spots for SRR2911306
# > Written 24463210 spots for SRR2911306
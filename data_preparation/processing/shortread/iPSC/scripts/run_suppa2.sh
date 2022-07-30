# scripts to install and run SUPPA2

#------
# SETUP
#------
# SUPPA2
# with conda
# conda setup
conda create -n suppa2_env python=3.4
conda activate suppa2_env
conda install -c conda-forge -c bioconda suppa=2.3

# start new env
screen -S suppa2_ipsc
conda activate suppa2_env

# specify paths
suppa_path2=~/anaconda3/envs/suppa2_env/bin/suppa.py
# test suppa
python $suppa2_path --help
# returns docstrings

# Salmon
conda install -c bioconda salmon

# make the salmon index files
# download gencode v35 transcripts
#curl https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_35/gencode.v35.transcripts.fa.gz -o gencode.v35.transcripts.fa.gz
#gzip -d gencode.v35.transcripts.fa.gz
fasta_fpath=eternity/data/RI_benchmarking_resources/gencode_v35_annotation_files/gencode.v35.transcripts.fa
salmonindex_outpath=eternity/data/RI_benchmarking_resources/salmon-index_grch38
salmon index -t $fasta_fpath -i $salmonindex_outpath

# get just retained introns event annotations
suppa2_outdir=/home/metamaden/ri_results/gb_revisions/suppa2/
sudo mkdir $suppa2_outdir
sudo chmod 777 $suppa2_outdir
# get annotation path
gtf_path=/eternity/data/RI_benchmarking_resources/gencode_v35_annotation_files/gencode.v35.annotation.gtf
sudo chmod 777 $gtf_path
# Generate the ioe files: 
rianno_outpath=$suppa2_outdir'gencode-v35-ri.events'
python $suppa_path2 generateEvents -i $gtf_path -o $rianno_outpath -e RI -f ioe
# returns:
#> INFO:eventGenerator:Reading input data.
#> INFO:eventGenerator:Calculating events
#> INFO:eventGenerator:Done

#-----------------------
# get tpm events -- ipsc
#-----------------------
# do tpm quantification
fastq_fpath1=eternity/data/RI_benchmarking_fastqs/SRR6026510_1.fastq
fastq_fpath1=eternity/data/RI_benchmarking_fastqs/SRR6026510_2.fastq
out_dpath=/home/metamaden/ri_results/gb_revisions/suppa2/SRR6026510/
sudo mkdir $out_dpath
sudo chmod 777 $out_dpath
sudo mkdir $out_dpath'logs'
sudo chmod 777 $out_dpath'logs'

# note: salmon needs full path strings, doesn't like bash variable evals
salmon quant -i eternity/data/RI_benchmarking_resources/salmon-index_grch38 \
-l ISF -1 eternity/data/RI_benchmarking_fastqs/SRR6026510_1.fastq \
-2 eternity/data/RI_benchmarking_fastqs/SRR6026510_2.fastq \
 -p 4 -o /home/metamaden/ri_results/gb_revisions/suppa2/SRR6026510/

# extract the tpm counts to a new file
salmonout_dpath=/home/metamaden/ri_results/gb_revisions/suppa2/SRR6026510/
mfspy_fpath=~/anaconda3/envs/suppa2_env/bin/multipleFieldSelection.py
python $mfspy_fpath -i $salmonout_dpath'quant.sf' -k 1 -f 4 -o $salmonout_dpath'iso_tpm.txt'

# format the tpm gene ids
# note: tpm has specific expected format that's expected
# first line is the sample ID, second line starts the data
# data should be formatted: transcript_id value (sep = "\t")
# process the tpm in R to make the new file 'iso_tpm_formatted.txt'
R
srrid <- "SRR6026510"
tpm.fpath <- file.path("home", "metamaden", "ri_results", "gb_revisions", 
                       "suppa2", srrid, "iso_tpm.txt")
tpm <- data.table::fread(tpm.fpath, header = F, sep = "|",
                         data.table = F)
tpm <- tpm[,c(1,9)]
tpm.newfpath <- file.path("home", "metamaden", "ri_results", "gb_revisions", 
                         "suppa2", srrid, "iso_tpm_formatted.txt")
data.table::fwrite(data.table::data.table(srrid), file = tpm.newfpath,
                   row.names = F, col.names = F, sep = "\t", quote = F,
                   append = F)
data.table::fwrite(tpm, file = tpm.newfpath, row.names = F, col.names = F, 
                   sep = "\t", quote = F, append = T)


# run the ri quantification
tpm_fpath=$salmonout_dpath'iso_tpm_formatted.txt'
ioe_fpath=/home/metamaden/ri_results/gb_revisions/suppa2/gencode-v35-ri.events_RI_strict.ioe
ritpm_outfname=ri_tpm_suppa2.txt
ritpm_outpath=/home/metamaden/ri_results/gb_revisions/suppa2/SRR6026510/$ritpm_outfname
python $suppa_path2 psiPerEvent -i $ioe_fpath -e $tpm_fpath -o $ritpm_outpath
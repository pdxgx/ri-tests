# notes_rmats.txt
#
# Documentation
# * GitHub tutorial: https://github.com/Xinglab/rmats-turbo/blob/v4.1.2/README.md

# Conda setup
# ```
# conda setup
conda create -n rmats_manual; conda activate rmats_manual
conda install -c conda-forge -c bioconda python=3.6.12 rmats=4.1.2
# ```

# ```
conda activate rmats
# ```

# Denote paired FASTQ files
# note: sep paired files with ':'
# ```
# /eternity/data/RI_benchmarking_fastqs/SRR6026510_1.fastq:/eternity/data/RI_benchmarking_fastqs/SRR6026510_2.fastq
# ```

# Make new output dirs
# ```
sudo mkdir /home/metamaden/ri_results/gb_revisions/rmats/
sudo chmod 777 /home/metamaden/ri_results/gb_revisions/rmats/

# ipsc
sudo mkdir /home/metamaden/ri_results/gb_revisions/rmats/SRR6026510/
sudo chmod 777 /home/metamaden/ri_results/gb_revisions/rmats/SRR6026510/

# hx1
sudo mkdir /home/metamaden/ri_results/gb_revisions/rmats/SRR2911306/
sudo chmod 777 /home/metamaden/ri_results/gb_revisions/rmats/SRR2911306/

sudo touch /home/metamaden/ri_results/gb_revisions/rmats/SRR6026510/out
sudo chmod 777 /home/metamaden/ri_results/gb_revisions/rmats/SRR6026510/out
# ```

#-------------------
# Run rMATS --- iPSC
#-------------------
# Run rMATs on FASTQs
# ```
vi /home/metamaden/ri_results/gb_revisions/rmats/SRR6026510/s1.txt
# paste:
# /eternity/data/RI_benchmarking_fastqs/SRR6026510_1.fastq:/eternity/data/RI_benchmarking_fastqs/SRR6026510_2.fastq

# specify paths
rmats_path=~/anaconda3/envs/rmats_manual/bin/rmats.py

# run rmats
python $rmats_path --s1 /home/metamaden/ri_results/gb_revisions/rmats/SRR6026510/s1.txt \
	--bi /eternity/data/RI_benchmarking_resources/gencode_v35_annotation_files/STAR_index_gencode_v35 \
	--gtf /eternity/data/RI_benchmarking_resources/gencode_v35_annotation_files/gencode.v35.primary_assembly.annotation.gtf \
	 --readLength 126 --nthread 10 \
	 --od /home/metamaden/ri_results/gb_revisions/rmats/SRR6026510 \
	 --tmp /home/metamaden/ri_results/gb_revisions/rmats/SRR6026510 \
	 2> /home/metamaden/ri_results/gb_revisions/rmats/SRR6026510/out

# ```
# rmats output -- new
# ```
# mapping the first sample
# mapping sample_0, /eternity/data/RI_benchmarking_fastqs/SRR6026510_1.fastq /eternity/data/RI_benchmarking_fastqs/SRR6026510_2.fastq is done with status 0
#         STAR --chimSegmentMin 2 --outFilterMismatchNmax 3 --alignEndsType EndToEnd --runThreadN 10 --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate --alignSJDBoverhangMin 6 --alignIntronMax 299999 --genomeDir /eternity/data/RI_benchmarking_resources/gencode_v35_annotation_files/STAR_index_gencode_v35 --sjdbGTFfile /eternity/data/RI_benchmarking_resources/gencode_v35_annotation_files/gencode.v35.primary_assembly.annotation.gtf --outFileNamePrefix /home/metamaden/ri_results/gb_revisions/rmats/SRR6026510/2022-06-25-17_05_38_380651_bam1_1/ --readFilesIn /eternity/data/RI_benchmarking_fastqs/SRR6026510_1.fastq /eternity/data/RI_benchmarking_fastqs/SRR6026510_2.fastq
#         STAR version: 2.7.10a   compiled: 2022-01-14T18:50:00-05:00 :/home/dobin/data/STAR/STARcode/STAR.master/source
# Jun 25 17:05:38 ..... started STAR run
# ...
# ...
# ==========
# Done processing each gene from dictionary to compile AS events
# Found 67018 exon skipping events
# Found 5703 exon MX events
# Found 18446 alt SS events
# There are 11155 alt 3 SS events and 7291 alt 5 SS events.
# Found 7317 RI events
# ==========
#
# ase: 2.786661148071289
# count: 4.277218341827393
# Processing count files.
# Done processing count files.
# ```
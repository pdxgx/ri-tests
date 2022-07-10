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

#------------------
# Run rMATS --- HX1
#------------------

Run rMATs on FASTQs
# ```
vi /home/metamaden/ri_results/gb_revisions/rmats/SRR2911306/s1.txt
# paste:
# eternity/data/RI_benchmarking_hx1/SRR2911306_1.fastq:eternity/data/RI_benchmarking_hx1/SRR2911306_2.fastq

python $rmats_path --s1 /home/metamaden/ri_results/gb_revisions/rmats/SRR2911306/s1.txt \
	--bi /eternity/data/RI_benchmarking_resources/gencode_v35_annotation_files/STAR_index_gencode_v35 \
	--gtf /eternity/data/RI_benchmarking_resources/gencode_v35_annotation_files/gencode.v35.primary_assembly.annotation.gtf \
	 --readLength 90 --nthread 10 \
	 --od /home/metamaden/ri_results/gb_revisions/rmats/SRR2911306 \
	 --tmp /home/metamaden/ri_results/gb_revisions/rmats/SRR2911306 \
	 2> /home/metamaden/ri_results/gb_revisions/rmats/SRR2911306/out
# ```

# returns
# ```
# Average number of transcripts per gene is 3.782377
# Average number of exons per transcript is 6.090643
# Average number of exons per transcript excluding one-exon tx is 6.718616
# Average number of gene per geneGroup is 8.485059
# statistic: 0.04446983337402344
# 
# read outcome totals across all BAMs
# USED: 38342200
# NOT_PAIRED: 0
# NOT_NH_1: 8361148
# NOT_EXPECTED_CIGAR: 316587
# NOT_EXPECTED_READ_LENGTH: 0
# NOT_EXPECTED_STRAND: 0
# EXON_NOT_MATCHED_TO_ANNOTATION: 681484
# JUNCTION_NOT_MATCHED_TO_ANNOTATION: 89549
# CLIPPED: 0
# total: 47790968
# outcomes by BAM written to: /home/metamaden/ri_results/gb_revisions/rmats/SRR2911306/2022-06-25-18_33_47_128045_read_outcomes_by_bam.txt
#
# novel: 184.88521146774292
# The splicing graph and candidate read have been saved into /home/metamaden/ri_results/gb_revisions/rmats/SRR2911306/2022-06-25-18_33_47_128045_*.rmats
# save: 1.478255033493042
# WARNING: The post step should use the same read length as the prep step.
#          The prep step's read length: 200
#          The post step's read length: 90
#          Please check /home/metamaden/ri_results/gb_revisions/rmats/SRR2911306/2022-06-19-18_52_05_771884_0.rmats
# loadsg: 0.1725466251373291
#
# ==========
# Done processing each gene from dictionary to compile AS events
# Found 58222 exon skipping events
# Found 4361 exon MX events
# Found 17771 alt SS events
# There are 10753 alt 3 SS events and 7018 alt 5 SS events.
# Found 7268 RI events
# ==========
#
# ase: 2.6181411743164062
# count: 3.007810592651367
# Processing count files.
# Done processing count files.
# ```
#
#
# Bottom of console output
# ```
# Jun 19 18:52:50 ..... inserting junctions into the genome indices
# Jun 19 18:54:19 ..... started mapping
# Jun 19 18:58:21 ..... finished mapping
# Jun 19 18:58:26 ..... started sorting BAM
# Jun 19 18:59:45 ..... finished successfully
# gtf: 25.962831735610962
# There are 60715 distinct gene ID in the gtf file
# There are 229647 distinct transcript ID in the gtf file
# There are 36806 one-transcript genes in the gtf file
# There are 1398698 exons in the gtf file
# ...
# ...
# Average number of gene per geneGroup is 8.485059
# statistic: 0.04403233528137207

# read outcome totals across all BAMs
# USED: 0
# NOT_PAIRED: 0
# NOT_NH_1: 8361148
# NOT_EXPECTED_CIGAR: 316587
# NOT_EXPECTED_READ_LENGTH: 39113233
# NOT_EXPECTED_STRAND: 0
# EXON_NOT_MATCHED_TO_ANNOTATION: 0
# JUNCTION_NOT_MATCHED_TO_ANNOTATION: 0
# CLIPPED: 0
# total: 47790968
# outcomes by BAM written to: /home/metamaden/ri_results/gb_revisions/rmats/SRR2911306/2022-06-19-18_52_05_771884_read_outcomes_by_bam.txt
#
# novel: 116.3007583618164
# The splicing graph and candidate read have been saved into /home/metamaden/ri_results/gb_revisions/rmats/SRR2911306/2022-06-19-18_52_05_771884_*.rmats
# save: 0.0001537799835205078
# loadsg: 0.21860122680664062
#
# ==========
# Done processing each gene from dictionary to compile AS events
# Found 50071 exon skipping events
# Found 3514 exon MX events
# Found 17085 alt SS events
# There are 10355 alt 3 SS events and 6730 alt 5 SS events.
# Found 7202 RI events
# ==========
#
# ase: 2.1956496238708496
# count: 0.5246806144714355
# Processing count files.
# Done processing count files.
# ```
# Description

Scripts to download and process short-read (SR) RNA-seq data for the sample HX1 (SRA run id SRR2911306).

# Script descriptions

* `fastqdump.sh` -- Download sample FASTQ files using `fastq-dump` function from the `sra-toolkit` software.
* `kma_annotation.sh` -- Makes the intron annotation files for the `KMA` tool. 
* `process_star_bam.sh` -- Align BAMs with `STAR`.
* `process_bowtie2_bam.sh` -- Align BAMs with `bowtie2`.

# Data download 

Sample FASTQs were downloaded 

## Scripts for data download
* `fastqdump.sh` -- Download sample FASTQ files using `fastq-dump` function from the `sra-toolkit` software.

# BAM preparation

The GENCODE v35 genome annotation files were used for alignments. For 4 of 5 SR RI-detection tools (`IntEREst`, `IRFinder-S`, `superintronic`, and `iREAD`), FASTQs were aligned with `STAR`, while the tool `KMA` uses BAMs aligned with `bowtie2`. In either case, BAMs were sorted using `samtools`. 

## Scripts for BAM preparation
* `kma_annotation.sh` -- Makes the intron annotation files for the `KMA` tool. 
* `process_star_bam.sh` -- Align BAMs with `STAR`.
* `process_bowtie2_bam.sh` -- Align BAMs with `bowtie2`.
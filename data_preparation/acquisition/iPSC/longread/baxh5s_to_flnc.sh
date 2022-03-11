#!/bin/bash

ACCNUM=$1
BARCODES=/eternity/data/RI_benchmarking_pacbio/SRP098984_SAMN07611993_human_iPS/barcode_files/custom_barcode_primers_all.fa

echo 'bax2bam'
time bax2bam -o "${ACCNUM}" *.bax.h5

echo 'starting conda env'
#conda activate isoseq-env
source activate isoseq-env

echo 'run ccs version 3.4'
time ccs --minPasses 1 --minPredictedAccuracy 0.90 "${ACCNUM}".subreads.bam "${ACCNUM}".ccs.bam

echo 'run lima'
time lima "${ACCNUM}".ccs.bam "$BARCODES" "${ACCNUM}"_demux.bam  --isoseq 

echo 'merge isoseq demuxed bams'
ls "${ACCNUM}"_demux.f_*.bam > bamlist.txt
split -l 1000 bamlist.txt splitbam
for bamlist in splitbam*; do
  samtools merge "${bamlist}"_bams.bam -b $bamlist
done
samtools merge "${ACCNUM}".demux.ccs.bam *_bams.bam
rm *_bams.bam
rm splitbam*
rm bamlist.txt
rm "${ACCNUM}"_demux.f_*
pbindex "${ACCNUM}".demux.ccs.bam

##samtools 
echo 'run isoseq refine to remove polyAs'
time isoseq3 refine --require-polya "${ACCNUM}".demux.ccs.bam "${BARCODES}" "${ACCNUM}".flnc.bam

echo 'deactivating venv'
#conda deactivate
source deactivate

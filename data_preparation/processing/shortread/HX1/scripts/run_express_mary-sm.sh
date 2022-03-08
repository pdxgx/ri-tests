#!/bin/bash

cd /eternity/data/RI_benchmarking_BAMs/

samples=`ls -1 *.SJ.out.tab  | sed -E 's/(.+)[.]SJ[.]out[.]tab/\1/'`

outdpath=/eternity/data/RI_benchmarking_BAMs/express_output
outfname=_sm 

fapath=/eternity/data/RI_benchmarking_resources/gencode_v35_annotation_files/trans_and_introns.fa

bamdpath=/eternity/data/RI_benchmarking_BAMs/RI_benchmarking_BAMs
bamfname=.kma.bam

for s in $samples; 
    do echo $s; 
        mkdir /eternity/data/RI_benchmarking_BAMs/express_output/$s; 
        express -o $outdpath/$s$outfname $fapath $bamdpath/$s$bamfname; 
    done
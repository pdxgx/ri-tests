#!/usr/bin/bash

# Author: Sean Maden
#
# Grabs SRR IDs with EDirect utils searches. Reads SRX expreiment IDs from 
# "concat_expts.csv" and iterates on esearch queries to get SRR IDs.
#

annofname=concat_expts.csv
exptidv=$(cat $annofname | cut -d, -f1 | sed 's/Experiment//g')
srrtname=srrtable.txt; rm $srrtname; touch $srrtname
for eid in ${exptidv[@]}
    do echo $eid
        esearch -db sra -query $eid" [Project]" | efetch -format xml | \
            xtract -pattern EXPERIMENT_PACKAGE -element PRIMARY_ID >> $srrtname
    done
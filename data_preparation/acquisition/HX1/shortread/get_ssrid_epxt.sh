#!/usr/bin/env sh

# Author: Sean Maden
#
# Scrape the project ID for an experiment from SRA.

exptid=
esearch -db sra -query $exptid" [Project]" | efetch -format xml | 
xtract -pattern EXPERIMENT_PACKAGE -element PRIMARY_ID
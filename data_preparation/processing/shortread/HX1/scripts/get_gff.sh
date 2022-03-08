#!/usr/bin/env bash

# Author Sean Maden

# Download gnome annotation GTF file

wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_35/gencode.v35.primary_assembly.annotation.gtf.gz
gunzip gencode.v35.primary_assembly.annotation.gtf.gz
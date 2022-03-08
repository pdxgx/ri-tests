#!/usr/bin/env sh

# Author: Sean Maden
#
# Download the intron spliceosomal annotations from IAOD/introndb.
# Target human annotations for genome build GRCh38/hg38.

wget https://introndb.lerner.ccf.org/static/bed/GRCh38_U12.bed
wget https://introndb.lerner.ccf.org/static/bed/GRCh38_U2.bed
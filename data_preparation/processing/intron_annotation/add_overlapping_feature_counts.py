#!/usr/bin/env python3

"""
add_overlapping_feature_counts.py
Python 3 code for adding exons per intron, intron base, and % intron overlap

SAMPLE RUN:
time python
../../intronomer-paper/benchmarking_data/add_overlapping_feature_counts.py
-f called_RIs/called_RI_data_summary_iPSC.tsv
-g ../files/gencode.v35.annotation.gtf

featureannotated files needed for plotting from called and all nonzero
summary files.

"""
import argparse
from datetime import datetime
from intervaltree import IntervalTree
import os
import pandas as pd


_ACCEPTABLE_CHROMS = {
    'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9',
    'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17',
    'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrM', 'chrX', 'chrY'
}


def create_exon_intervaltree(gtf_file):
    exon_intervals = {chrom: IntervalTree() for chrom in _ACCEPTABLE_CHROMS}
    with open(gtf_file) as gtf:
        # Parse valid exon lines from the annotation file into a dict by
        # transcript_id
        for line in gtf:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            if '#' in line:
                line = line.split('#')[0].strip()
            try:
                item = line.split('\t')
                chrom, _, feature, left, right = item[0:5]
            except ValueError:
                continue
            left, right = int(left), int(right)
            if feature == 'exon' and left < right:
                try:
                    exon_intervals[chrom][left:right+1] = ''
                except KeyError:
                    continue
    return exon_intervals


def max_features_per_base(df_row, exon_intervals):
    max_count = 0
    for pos in range(df_row['left'], df_row['right']):
        max_count = max(max_count, len(exon_intervals[df_row['chrom']][pos]))
    return max_count


def percent_of_bases_covered(df_row, exon_intervals):
    total_length = df_row['right'] - df_row['left']
    bases_covered = 0
    for pos in range(df_row['left'], df_row['right']):
        if len(exon_intervals[df_row['chrom']][pos]) > 0:
            bases_covered += 1
    return bases_covered / total_length


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Adds annotated exon overlap information to intron file.'
    )
    parser.add_argument(
        '--file-to-annotate', '-f',
        help='summary or other grange file to annotated with exon overlaps '
    )
    parser.add_argument(
        '--annotation-gtf', '-g',
        help='.gtf file with GENCODE annotation (v35).'
    )

    args = parser.parse_args()
    summary_file = args.file_to_annotate
    gtf_path = args.annotation_gtf

    now = datetime.now().strftime('%m-%d-%Y_%H.%M.%S')

    exon_intervals = create_exon_intervaltree(gtf_path)
    tx_df = pd.read_table(summary_file, sep='\t')

    orig_cols = tx_df.columns.values.tolist()
    tx_df['chrom'] = tx_df['intron'].apply(lambda x: x.split(':')[0])
    tx_df['left'] = tx_df['intron'].apply(
        lambda x: int(x.split(':')[1].split('-')[0])
    )
    tx_df['right'] = tx_df['intron'].apply(
        lambda x: int(x.split(':')[1].split('-')[1]) + 1
    )
    tot_feats = 'total_overlapping_features'
    max_feats = 'max_features_per_base'
    perc_bases = '%_bases_overlapped'
    tx_df['total_overlapping_features'] = tx_df.apply(
        lambda x: len(exon_intervals[x['chrom']][x['left']:x['right']]), axis=1
    )
    tx_df['max_features_per_base'] = tx_df.apply(
        lambda x: max_features_per_base(x, exon_intervals), axis=1
    )
    tx_df['%_bases_overlapped'] = tx_df.apply(
        lambda x: percent_of_bases_covered(x, exon_intervals), axis=1
    )

    new_cols = orig_cols + [tot_feats, max_feats, perc_bases]
    tx_df = tx_df[new_cols]
    out_dir = os.path.dirname(summary_file)
    filename = os.path.basename(summary_file).split('.tsv')[0]
    filename += 'featureannotated'
    outfile = os.path.join(out_dir, '{}.tsv'.format(filename))
    with open(outfile, 'w') as output:
        tx_df.to_csv(output, index=False, sep='\t')

#!/usr/bin/env python3

"""
annotated_shortread_ranges.py
Python 3 code for mapping aligned pacbio reads to annotated transcripts

SAMPLE RUN:
time python ../intronomer-paper/benchmarking_data/annotate_shortread_ranges.py
-l HX1_final/processed_tx_df_HX1_02-21-2022_21.03.36.csv
-s HX1_final/granges-lrmap_sr-5-methods_SRR2911306-hx1.csv
-f ../immunotherapy/files/GRCh38.primary_assembly.genome.fa
-L HX1_final/reads_per_gene_and_transcript_HX1.tsv
-C HX1_final/mpile-sstat_gencode-v35_SRR2911306-hx1.csv


"""
import argparse
from datetime import datetime
import os
import pandas as pd
import subprocess as sp


_ACCEPTABLE_CHROMS = {
    'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9',
    'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17',
    'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrM', 'chrX', 'chrY'
}


def jx_to_motif(jx, reference_fasta, samtools_path, hyphen=False):
    """Given an input junction and reference genome, transcribes RNA sequence.

    Input:
    jx: a junction in 'chr_;left;right;strand' format, with 0-based fully
        closed coordinates, (string)
    reference_fasta: reference genome fasta file, previously sorted and
        indexed by samtools (path to fasta file, string)
    samtools_path: to be called from subprocess to collect the sequence (path
        to samtools executable, string)

    Returns a left-to-right nucleotide sequence on either side of the aberrant
    junction sufficiently long to generate the desired protein sequence.
    """
    chrom = jx.split(':')[0]
    left = jx.split(':')[1].split('-')[0]
    right = jx.split('-')[1]
    if chrom not in _ACCEPTABLE_CHROMS:
        return ''
    left_start = left
    left_stop = str(int(left) + 1)
    left_range = chrom + ':' + left_start + '-' + left_stop
    left_output = sp.check_output(
        ['{}'.format(samtools_path), 'faidx', '{}'.format(reference_fasta),
         '{}'.format(left_range)]
    )
    left_seq = ''.join(left_output.decode("utf-8").splitlines()[1:])
    right_start = str(int(right) - 1)
    right_stop = right
    right_range = chrom + ':' + right_start + '-' + right_stop
    right_output = sp.check_output(
        ['{}'.format(samtools_path), 'faidx', '{}'.format(reference_fasta),
         '{}'.format(right_range)]
    )
    right_seq = ''.join(right_output.decode("utf-8").splitlines()[1:])
    sequence = left_seq + right_seq
    if hyphen:
        sequence = sequence[:2] + '-' + sequence[-2:]
    return sequence


def annotate_shortread_lrmap(tx_df, ranges_file, now, out_dir, ref_fasta,
                             samtools, l_counts, s_counts):
    annotated_outfile = os.path.join(
        output_dir, 'LR_annotated_{}'.format(os.path.basename(all_ranges))
    )
    sr_ct = 'short_read_gene_median_coverage'
    lr_gene_ct = 'long_reads_per_gene'
    lr_tx_ct = 'long_reads_per_transcript'
    count_df = pd.read_table(l_counts, sep='\t')
    count_df.rename(
        {'reads_per_transcript': lr_tx_ct, 'reads_per_gene': lr_gene_ct},
        axis=1, inplace=True
    )
    sr_count_df = pd.read_table(
        s_counts, usecols=['gene.id', 'q50'], sep=',',
    )
    sr_count_df['q50'].fillna(0)
    sr_count_df.rename(
        {'gene.id': 'gene_id', 'q50': sr_ct}, axis=1, inplace=True
    )
    count_df = pd.merge(count_df, sr_count_df, how='outer', on='gene_id')
    target_genes = set(count_df.loc[
        (count_df[sr_ct] >= 2) & (count_df[lr_gene_ct] >= 5)
        & (count_df[lr_tx_ct] >= 5)
    ]['gene_id'].tolist())
    count_df = count_df.loc[(count_df[sr_ct] > 0) | (count_df[lr_gene_ct] > 0)]
    canonical_motifs = {
        'GTAG', 'GCAG', 'ATAC', 'CTAC', 'CTGC', 'GTAT'
    }
    result_df = pd.read_table(ranges_file, sep=',')

    result_df['intron'] = result_df.apply(
        lambda x: '{}:{}-{}'.format(x['seqnames'], x['start'], x['end']),
        axis=1
    )
    result_df = pd.merge(result_df, tx_df, how='outer', on='intron')
    result_df.rename(
        {'position': 'intron_position_in_tx'}, axis=1, inplace=True
    )
    perst = 'persistence'
    result_df['max_intron_persistence'] = result_df['intron'].apply(
        lambda x: result_df.loc[result_df['intron'] == x][perst].max()
    )
    result_df['motif'] = result_df['intron'].apply(
        lambda x: jx_to_motif(x, ref_fasta, samtools)
    )
    result_df['canonical_motif'] = result_df['motif'].apply(
        lambda x: int(x in canonical_motifs)
    )
    result_df.dropna(subset=['transcript'], axis=0, inplace=True)
    result_df = pd.merge(result_df, count_df, how='left', on='transcript')
    with open(annotated_outfile, 'w') as output:
        result_df.to_csv(output, sep=',', index=False)
    result_df = result_df.loc[result_df['gene_id'].isin(target_genes)]
    mini_outfile = os.path.join(
        output_dir, 'target_genes_LR_annotated_{}'.format(
            os.path.basename(all_ranges)
        )
    )
    with open(mini_outfile, 'w') as output:
        result_df.to_csv(output, sep=',', index=False)

    tx_df[perst] = pd.to_numeric(tx_df[perst], errors='coerce').notnull()
    print('reading from LR')

    tx_df = tx_df.loc[tx_df[perst] > 0].copy()
    max_idx = tx_df.groupby(['intron'])[perst].transform(max) == tx_df[perst]
    tx_df = tx_df[max_idx]
    outfile = os.path.join(out_dir, 'unique_max_RIS_{}.csv'.format(now))
    with open(outfile, 'w') as output:
        tx_df.to_csv(output, index=False, sep=',')
    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Extracts retained-intron exons from stringtie-assembled '
                    'transcripts.'
    )
    parser.add_argument(
        '--long-read-isoform-results', '-l',
        help='output file from assign_RI_persistence_metric.py, '
             '"processed_tx_df_[timestamp].csv" '
    )
    parser.add_argument(
        '--range-summarized-shortread-results', '-s',
        help='short-read detection results summarized by intron region'
    )
    parser.add_argument(
        '--samtools-path', '-S', default='samtools',
        help='Give the path to access samtools.'
    )
    parser.add_argument(
        '--reference-genome', '-f',
        help='.fa file containing the reference genome sequence. NOTE: This '
             'fasta file must previously have been indexed by running '
             '"samtools faidx <ref.fasta>".'
    )
    parser.add_argument(
        '--long-read-counts', '-L',
        help='"reads_per_gene_and_transcript_[sample].tsv" file generated by '
             'script pacbio_reads_to_transcripts.py'
    )
    parser.add_argument(
        '--short-read-counts', '-C',
        help='"mpile-sstat_gencode-v35_SRR6026510_[sample].csv" file'
    )

    args = parser.parse_args()
    lr_file = args.long_read_isoform_results
    all_ranges = args.range_summarized_shortread_results
    samtools = args.samtools_path
    ref_fasta = args.reference_genome
    lr_counts = args.long_read_counts
    sr_counts = args.short_read_counts

    now = datetime.now().strftime('%m-%d-%Y_%H.%M.%S')
    output_dir = os.path.dirname(all_ranges)
    tx_df = pd.read_table(lr_file, sep=',')
    annotate_shortread_lrmap(
        tx_df, all_ranges, now, output_dir, ref_fasta, samtools,
        lr_counts, sr_counts
    )

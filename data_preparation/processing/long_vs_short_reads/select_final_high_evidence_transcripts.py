#!/usr/bin/env python3

"""
select_final_high_evidence_transcripts.py
Python 3 code for selecting the stringent length-matched transcripts

SAMPLE RUN:
time python
ri-tests/data_preparation/processing/long_vs_short_reads/
select_high_evidence_txs.py
-t iPSC_8tools/RI_txs_to_read_ids_final_01-06-2022_17.05.22_iPSC.tsv.gz
-m iPSC_8tools/
    target_genes_LR_annotated_granges-lrmap_sr-8-methods_SRR6026510-ipsc.csv
-o iPSC_8tools/subset_genes
"""
import argparse
from datetime import datetime
import os
import matplotlib.pyplot as plt
import pandas as pd
from statistics import median


def histogram_panels(input_data, outdir, log_x=False, bincount='auto'):
    fig, axs = plt.subplots(1, len(input_data))
    plt.subplots_adjust(wspace=.15, hspace=.075)
    fig.set_figheight(4)
    fig.set_figwidth(10)
    ylabel_rot = 90
    ax_label_fontsize = 8
    y_labels = [
        'transcript-read pairs (#)',
        'transcript-read pairs (#)'
    ]
    x_label = 'length difference (# bp)'
    flag = '{}bins'.format(bincount)
    if log_x:
        flag += '_logx'
    legend_title_flags = ['pairs (left ends)', 'pairs (right ends)']
    for i, datalist in enumerate(input_data):
        datalist.sort()
        ax = axs[i]
        plt.sca(ax)
        plt.hist(datalist, bins=bincount, log=False)
        ax.xaxis.grid(False)
        ax.yaxis.grid(True)
        ax.axvline(
            -50, color='red', linewidth=.5,
            label='+/- 50 bases'.format(median(datalist))
        )
        ax.axvline(
            50, color='red', linewidth=.5,
        )
        ax.legend(title='{} {}:'.format(len(datalist), legend_title_flags[i]))
        plt.ylabel(y_labels[i], fontsize=ax_label_fontsize)
        plt.setp(
            ax.yaxis.get_majorticklabels(),
            fontsize=ax_label_fontsize, rotation=ylabel_rot
        )
        plt.xlabel(x_label, fontsize=ax_label_fontsize)
        if log_x:
            plt.xscale('log')
        else:
            ax.set_xlim(-1000, 1000)
    now = datetime.now().strftime('%m-%d-%Y_%H.%M.%S')
    out_file = os.path.join(
        outdir, 'tx_vs_read_lengths_{}_{}.pdf'.format(flag, now)
    )
    fig = plt.gcf()
    fig.savefig(out_file, bbox_inches='tight', pad_inches=0.1)
    return


def scatter_regression_and_size(scatter_df, yval, xval, figfile):
    """

    :param scatter_df:
    :param xval:
    :param yval:
    :param out_path:
    :param now:
    :param y_axis_label:
    :param x_axis_label:
    :param chrom:
    :return:
    """
    x_axis_label = 'left coordinate difference\n(transcript - read)'
    y_axis_label = 'right coordinate difference\n(read - transcript)'

    color_dict = {'HX1': '#d3494e', 'iPSC': '#448ee4'}
    groups = scatter_df.groupby([xval, yval])
    scatter_dict = {'x': [], 'y': [], 'count': [], 'color': []}
    for group_index, group in groups:
        scatter_dict['x'].append(group_index[0])
        scatter_dict['y'].append(group_index[1])
        # scatter_dict['count'].append(6 * (len(group) ** (1.0 / 1.5)))
        scatter_dict['count'].append(1 * (len(group) ** (1.0 / 1.5)))
        # scatter_dict['color'].append(color_dict[group_index[2]])
        scatter_dict['color'].append('#887191')
    scatter_df = pd.DataFrame(scatter_dict)

    plt.subplots()
    plt.rcParams.update({'figure.autolayout': True})
    plt.rcParams['figure.figsize'] = 5.0, 5.0
    plt.scatter(
        x='x', y='y', c='color', data=scatter_df, s='count',
        edgecolors='#d8dcd6', linewidths=0.05, label=None, #alpha=0.9
    )
    ax = plt.gca()
    ax.xaxis.grid(True)
    ax.yaxis.grid(True)
    ax.set_xlim(-250, 250)
    ax.set_ylim(-250, 250)
    # ax.set_xlim(-1000, 1000)
    # ax.set_ylim(-1000, 1000)
    plt.xlabel(x_axis_label)
    # ax.set_yscale('log')
    # ax.set_xscale('log')
    plt.ylabel(y_axis_label)
    fig = plt.gcf()
    fig.savefig(figfile)
    plt.clf()
    return


def select_high_evidence_genes(tx_to_read_id, merged_file, gene_file,
                               output_dir, now, sub_dir):
    gene_df = pd.read_table(gene_file, header=None, names=['gene_id'])
    gene_set = set(gene_df['gene_id'].unique())
    orig_df = pd.read_table(tx_to_read_id)
    print(
        '{} total original genes'.format(orig_df['gene_id'].nunique())
    )
    tx_to_read_id_df = orig_df.loc[
        orig_df['gene_id'].isin(gene_set)
    ]
    print(
        'after original gene set: {} total genes'.format(
            tx_to_read_id_df['gene_id'].nunique()
        )
    )
    full_matches = {
        'skipped_splicing-full_length', 'all_introns-full_length'
    }
    print(len(tx_to_read_id_df))
    tx_to_read_id_df = tx_to_read_id_df.loc[
        (tx_to_read_id_df['match_type'].isin(full_matches))
    ]
    print(
        '{} genes with full length matches genes'.format(
            tx_to_read_id_df['gene_id'].nunique()
        )
    )
    rd_l = 'read_left'
    rd_r = 'read_right'
    tx_l = 'tx_left'
    tx_r = 'tx_right'
    l_diff = 'left_diff'
    r_diff = 'right_diff'
    tx_to_read_id_df[l_diff] = tx_to_read_id_df.apply(
        lambda x: x[tx_l] - x[rd_l], axis=1
    )
    tx_to_read_id_df[r_diff] = tx_to_read_id_df.apply(
        lambda x: x[rd_r] - x[tx_r], axis=1
    )
    histogram_panels(
        [tx_to_read_id_df[l_diff].tolist(), tx_to_read_id_df[r_diff].tolist()],
        output_dir, bincount=10000
    )
    max_diff = 50
    tx_to_read_id_df = tx_to_read_id_df.loc[
        (abs(tx_to_read_id_df[l_diff] <= max_diff))
        & (abs(tx_to_read_id_df[r_diff] <= max_diff))
    ]
    print('{} after filtering on coordinates:'.format(
        tx_to_read_id_df['gene_id'].nunique())
    )
    target_txs = set()
    for tx in tx_to_read_id_df['transcript'].unique():
        mini_df = tx_to_read_id_df.loc[
            tx_to_read_id_df['transcript'] == tx
        ]
        if len(mini_df) > 4:
            target_txs.add(tx)
    print('{} target genes out of {} original genes'.format(
        tx_to_read_id_df['gene_id'].nunique(), orig_df['gene_id'].nunique()
    ))
    final_target_genes = tx_to_read_id_df['gene_id'].unique()
    gene_path = os.path.join(output_dir, 'subset_target_genes.txt')
    with open(gene_path, 'w') as output:
        for gene in final_target_genes:
            output.write('{}\n'.format(gene))
    print(len(target_txs))
    full_merged_df = pd.read_table(merged_file, sep=',')
    target_genes = list(full_merged_df.loc[
        full_merged_df['transcript'].isin(target_txs)
    ]['gene_id'].unique())
    print(full_merged_df['transcript'].nunique())
    print('{} final target genes'.format(len(target_genes)))
    subset_merged_df = full_merged_df.loc[
        full_merged_df['gene_id'].isin(target_genes)
    ]
    out_file = os.path.join(sub_dir, 'subset_' + os.path.basename(merged_file))
    subset_merged_df.to_csv(out_file, sep=',', index=False)
    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Selects final stringent subset of transcripts to study.'
    )
    parser.add_argument(
        '--tx-to-read-id', '-t',
        help='RIs_to_read_ids[timestamp] file generated by script '
             'pacbio_reads_to_transcripts.py.'
    )
    parser.add_argument(
        '--merged-file-to-subset', '-m',
        help='target_genes_LR_annotated_granges-lrmap.... file generated by '
             'merging long and short read data by '
             'annotate_shortread_ranges.py script.'
    )
    parser.add_argument(
        '--orig-gene-file', '-r',
        help='target_genes_[sample].tsv'
    )
    parser.add_argument(
        '--subset-directory', '-o',
        help='directory to write subset dataframe output'
    )

    args = parser.parse_args()
    tx_read_file = args.tx_to_read_id
    merged_file = args.merged_file_to_subset
    sub_dir = args.subset_directory
    gene_file = args.orig_gene_file

    now = datetime.now().strftime('%m-%d-%Y_%H.%M.%S')
    output_dir = os.path.dirname(tx_read_file)

    select_high_evidence_genes(
        tx_read_file, merged_file, gene_file, output_dir, now, sub_dir
    )

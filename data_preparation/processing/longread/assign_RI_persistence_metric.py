#!/usr/bin/env python3

"""
assign_RI_persistence_metric.py
Python 3 code for assigning read-transcript maps with intron persistence values

SAMPLE RUN:
time python intronomer-paper/benchmarking_data/assign_RI_persistence_metric.py
-m HX1_final/RI_txs_to_read_ids_final_01-06-2022_16.26.46.tsv

Add -p for persistence-vs-intron position plotting
"""
import argparse
from datetime import datetime
from math import ceil
import matplotlib.pyplot as plt
from matplotlib import style
import numpy as np
import os
import pandas as pd
from scipy import stats
from scipy.spatial.distance import cdist


_ALL_CHROMS = {
    'chr1': '#3c73a8', 'chr2': '#5c8b15', 'chr3': '#06c2ac',
    'chr4': '#f7879a', 'chr5': '#019529', 'chr6': '#8eab12',
    'chr7': '#c65102', 'chr8': '#ff474c', 'chr9': '#fac205',
    'chr10': '#c0737a', 'chr11': '#5e819d', 'chr12': '#6d5acf',
    'chr13': '#82a67d', 'chr14': '#be0119', 'chr15': '#e17701',
    'chr16': '#028f1e', 'chr17': '#0d75f8', 'chr18': '#fe019a',
    'chr19': '#9e43a2', 'chr20': '#6f828a', 'chr21': '#98568d',
    'chr22': '#01889f', 'chrX': '#ff69af', 'chrY': '#89a0b0',
}


def create_intron_matrix(transcript_df):
    transcript_df.sort_values(
        by="read_introns", key=lambda x: x.str.len(),
        ascending=False, inplace=True
    )
    tx_introns = transcript_df['tx_introns'].unique()[0]
    tx_introns = [intron for intron in tx_introns.split(';') if intron]
    if transcript_df['strand'].unique()[0] == '-':
        tx_introns = tx_introns[::-1]
    plot_list = []
    reads = []
    for index, row in transcript_df.iterrows():
        plot_introns = []
        reads.append(row['read_id'])
        read_ints = row['read_introns'].split(';')
        for intron in tx_introns:
            if intron in read_ints:
                plotval = 1
            else:
                i_left, i_right = map(int, intron.split('-'))
                intron_in_read = (
                        i_left > row['read_left']
                        and i_right < row['read_right']
                )
                if intron_in_read:
                    plotval = 0
                else:
                    plotval = np.NaN
            plot_introns.append(plotval)
        plot_list.append(plot_introns)

    intron_df = pd.DataFrame(
        np.array(plot_list), columns=tx_introns, index=reads
    )
    return intron_df


def twodigit_string(value):
    if value == 0:
        return '0'
    elif value >= 0.1:
        return '{:.2g}'.format(value)
    else:
        return '{:.1g}'.format(value)


def hamming_similarity(df_row, df):
    hs = []
    df_row.rename('orig', inplace=True)
    for index, row in df.iterrows():
        r_df = df_row.to_frame().join(row)
        r_df.dropna(axis=0, inplace=True)
        hs.append(cdist(
            np.array([r_df['orig']]), np.array([r_df[index]]), 'hamming'
        ))
    return 1 - np.mean(hs)


def assign_RI_metrics(read_tx_df, output_dir, now, batch_num, flag=''):
    read_tx_df = read_tx_df.loc[read_tx_df['RI-read_ratio'] > 0].copy()
    all_intron_dict = {
        'intron': [], 'transcript': [], 'persistence': [], 'position': [],
        'chrom': []
    }
    metric_check_dict = {
        'intron': [], 'transcript': [],
        'matched Hamming similarity only': []
    }
    for tx in read_tx_df['transcript'].unique():
        sub_df = read_tx_df.loc[read_tx_df['transcript'] == tx].copy()
        chrom = sub_df['chrom'].unique()[0]
        intron_df = create_intron_matrix(sub_df)
        num_introns = len(intron_df.columns) - 1
        n_x = len(intron_df)
        for int_count, intron in enumerate(intron_df.columns.values):
            int_df = intron_df.dropna(subset=[intron]).copy()
            if len(int_df) == 0:
                persistence = 0
            elif int_df[intron].sum() == int_df[intron].count():
                persistence = 0
            else:
                # Percent of transcript reads with intron information
                info_perc = len(int_df) / n_x
                # Only need to calculate splicing progression and similarity
                # if the intron is retained in the read
                int_df = int_df.loc[int_df[intron] == 0]
                # Splicing progression ignoring the target intron
                int_df['p_xy'] = (
                    int_df.drop([intron], axis=1).sum(axis=1)
                    / int_df.drop([intron], axis=1).count(axis=1)
                )
                # Setup for multiplying in read-specific hamming distance
                int_df['h_sim'] = int_df.apply(
                    lambda x: hamming_similarity(x, int_df), axis=1
                )
                int_df['IR_persistence'] = int_df.apply(
                    lambda x: x['p_xy'] * (1 - x[intron]) * x['h_sim'],
                    axis=1
                )
                persistence = info_perc * int_df['IR_persistence'].sum() / n_x
            metric_check_dict['intron'].append(
                '{}:{}'.format(chrom, intron)
            )
            metric_check_dict['matched Hamming similarity only'].append(
                persistence
            )
            metric_check_dict['transcript'].append(tx)

            all_intron_dict['intron'].append('{}:{}'.format(chrom, intron))
            all_intron_dict['chrom'] = chrom
            all_intron_dict['transcript'].append(tx)
            all_intron_dict['persistence'].append(persistence)
            all_intron_dict['position'].append(int_count / num_introns)

    metric_check_df = pd.DataFrame(metric_check_dict)
    out_name = 'processed_tx_df_{}updatecheck_{}.csv'.format(flag, now)
    print_header=True
    if batch_num is not None:
        out_name = 'batch{}_'.format(batch_num) + out_name
        if batch_num > 0:
            print_header = False
    twometrics_outfile = os.path.join(output_dir, out_name)
    metric_check_df.to_csv(
        twometrics_outfile, index=False, sep=',', header=print_header
    )

    tx_df = pd.DataFrame(all_intron_dict)
    print_df = tx_df[['intron', 'transcript', 'persistence', 'position']]
    out_name = 'processed_tx_df_{}_{}.csv'.format(flag, now)
    if batch_num is not None:
        out_name = 'batch{}_'.format(batch_num) + out_name
    outfile = os.path.join(output_dir, out_name)
    print_df.to_csv(outfile, index=False, sep=',', header=print_header)
    return tx_df


def scatter_with_regression_and_size(scatter_df, xval, yval, figfile, chrom='',
                                     x_axis_label="5' --> 3' intron position",
                                     y_axis_label='IR persistence',
                                     flag=''):
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
    color_dict = {'HX1': '#d3494e', 'iPSC': '#448ee4'}
    x = scatter_df[xval]
    y = scatter_df[yval]
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    line = slope * x + intercept
    print('linear regression stats:')
    print(slope, intercept, r_value, p_value, std_err)
    line_label = (
        'y = {}x + {}\nR = {}'
        ''.format(round(slope, 3), round(intercept, 3), round(r_value, 2))
    )
    groups = scatter_df.groupby([xval, yval])
    scatter_dict = {'x': [], 'y': [], 'count': [], 'color': []}
    for group_index, group in groups:
        scatter_dict['x'].append(group_index[0])
        scatter_dict['y'].append(group_index[1])
        scatter_dict['count'].append(6 * (len(group) ** (1.0 / 1.5)))
        if chrom:
            scatter_dict['color'].append(_ALL_CHROMS[chrom])
        elif flag:
            scatter_dict['color'].append(color_dict[flag])
        else:
            scatter_dict['color'].append('#b9484e')
    scatter_df = pd.DataFrame(scatter_dict)

    ax_max = max(scatter_df.x.max(), scatter_df.y.max()) * 1.1
    point_dict = {'fake_x': [], 'fake_y': [], 'color': []}
    for color in scatter_dict['color']:
        point_dict['fake_x'].append(-ax_max * 10)
        point_dict['fake_y'].append(-ax_max * 10)
        point_dict['color'].append(color)
    fake_df = pd.DataFrame(point_dict)

    plt.rcParams.update({'figure.autolayout': True})
    plt.rcParams['figure.figsize'] = 5.0, 5.0
    style.use('seaborn-whitegrid')
    plt.scatter(
        x='x', y='y', c='color', data=scatter_df, s='count',
        edgecolors='#d8dcd6', linewidths=0.05, label=None,
    )
    if chrom:
        legend_label = chrom
    else:
        legend_label = 'all chromosomes'
    plt.scatter(
        x='fake_x', y='fake_y', c='color', data=fake_df, s=40,
        edgecolors='xkcd:light grey', linewidth=0.01, label=legend_label
    )
    plt.plot(
        x, line, c='xkcd:dark grey', linewidth=0.8, label=line_label
    )
    plt.legend(fontsize='small', frameon=True)

    plt.xlabel(x_axis_label)
    ax = plt.gca()
    y_max = scatter_df.y.max() * 1.1
    x_max = scatter_df.x.max() * 1.1
    ax.set_ylim(ymin=-(0.05 * y_max), ymax=y_max)
    ax.set_xlim(xmin=-(0.05 * x_max), xmax=x_max)

    plt.ylabel(y_axis_label)
    fig = plt.gcf()
    fig.savefig(figfile)
    plt.clf()
    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Calculates intron persistence metric values.'
    )
    parser.add_argument(
        '--read-transcript-map', '-m',
        help='file with name "RI_txs_to_read_ids_final__[timestamp]__.tsv" '
             'generated previously by this script, mapping reads to target '
             'best-match transcripts.'
    )
    parser.add_argument(
        '--plot-persistence-vs-intron-position', '-p', action='store_true',
        help='select this option to create persistence-vs-intron-position '
             'plots'
    )
    parser.add_argument(
        '--batch-number', type=int,
        help='For running in batches on slurm: the slurm array task id to '
             'process correct chunk of the dataframe.'
    )
    parser.add_argument(
        '--total-batches', type=int,
        help='total number of batches for batch job on slurm'
    )

    args = parser.parse_args()
    read_map = args.read_transcript_map
    plot_persistence_vs_pos = args.plot_persistence_vs_intron_position
    batch_num = args.batch_number
    tot_batches = args.total_batches

    now = datetime.now().strftime('%m-%d-%Y_%H.%M.%S')
    output_dir = os.path.dirname(read_map)
    read_tx_df = pd.read_table(read_map)

    if 'iPSC' in output_dir:
        flag = 'iPSC'
    else:
        flag = 'HX1'

    if batch_num is not None:
        all_txs = read_tx_df['transcript'].unique().tolist()
        all_txs.sort()
        chunksize = ceil(len(all_txs) / tot_batches)
        tx_lists = [
            all_txs[i:i + chunksize]
            for i in range(0, len(all_txs), chunksize)
        ]
        target_txs = tx_lists[batch_num]
        print(
            'Batch {}: {} total transcripts, {} in current batch'.format(
                batch_num, len(all_txs), len(target_txs)
            )
        )
        read_tx_df = read_tx_df.loc[read_tx_df['transcript'].isin(target_txs)]
        print(
            'total length of df: {}, {} txs'.format(
                len(read_tx_df), read_tx_df['transcript'].nunique()
            )
        )

    tx_df = assign_RI_metrics(
        read_tx_df, output_dir, now, batch_num, flag=flag
    )

    if plot_persistence_vs_pos:
        if batch_num is not None:
            fig_file = os.path.join(
                output_dir,
                'batch{}_all_chroms_intron_retention_persistence_vs_position'
                '_{}.pdf'.format(batch_num, now)
            )
        else:
            fig_file = os.path.join(
                output_dir,
                'all_chroms_intron_retention_persistence_vs_position_{}.pdf'
                ''.format(now)
            )
        scatter_with_regression_and_size(
            tx_df, 'position', 'persistence', fig_file, flag=flag
        )
        # chromosome-by-chromosome
        for chr in _ALL_CHROMS.keys():
            sub_df = tx_df.loc[tx_df['chrom'] == chr].copy()
            if len(sub_df) > 0:
                # no "chrY" in sample
                if batch_num:
                    fig_file = os.path.join(
                        output_dir,
                        'batch{}_{}_intron_retention_persistence_vs_position'
                        '_{}.pdf'.format(batch_num, chr, now)
                    )
                else:
                    fig_file = os.path.join(
                        output_dir,
                        '{}_intron_retention_persistence_vs_position_{}.pdf'
                        ''.format(chr, now)
                    )
                scatter_with_regression_and_size(
                    tx_df, 'position', 'persistence', fig_file, chrom=chr
                )

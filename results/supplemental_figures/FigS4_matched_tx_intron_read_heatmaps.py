#!/usr/bin/env python3

"""
intron_read_heatmaps_all_matches.py
Python 3 code for plotting persistence heatmaps for iPSC and HX1 matched txs

time python
../intronomer-paper/benchmarking_data/intron_read_heatmaps_all_matches.py
-H HX1_final/RI_txs_to_read_ids_final_01-06-2022_16.26.46.tsv
-i iPSC_final/RI_txs_to_read_ids_final_01-06-2022_17.05.22.tsv
-x HX1_final/reads_per_gene_and_transcript_HX1.tsv
-p iPSC_final/reads_per_gene_and_transcript_iPSC.tsv

"""
import argparse
from datetime import datetime
import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
from scipy.spatial.distance import cdist


def heatmap(data, row_labels, col_labels, ax=None, cbar_kw={},
            cbarlabel="", bot_labels=[], lgnd_info=[], notop=False,
            yaxis_label='reads', xaxis_top_label='introns',
            xaxis_bottom_label='IR persistence', **kwargs):
    """
    Create a heatmap from a numpy array and two lists of labels.

    Parameters
    ----------
    data
        A 2D numpy array of shape (N, M).
    row_labels
        A list or array of length N with the labels for the rows.
    col_labels
        A list or array of length M with the labels for the columns.
    ax
        A `matplotlib.axes.Axes` instance to which the heatmap is plotted.  If
        not provided, use current axes or create a new one.  Optional.
    cbar_kw
        A dictionary with arguments to `matplotlib.Figure.colorbar`.  Optional.
    cbarlabel
        The label for the colorbar.  Optional.
    **kwargs
        All other arguments are forwarded to `imshow`.
    """
    if not ax:
        ax = plt.gca()

    # Plot the heatmap
    im = ax.imshow(data, **kwargs)

    if not lgnd_info:
        # Create colorbar
        cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
        cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")

    ax.set_xticks(np.arange(data.shape[1]))
    ax.set_yticks(np.arange(data.shape[0]))
    ax.set_xticklabels(col_labels, rotation=90, fontsize=8)
    ax.set_yticklabels(row_labels)
    ax.set_ylabel(yaxis_label, fontsize=8)
    ax.set_xlabel(xaxis_bottom_label, fontsize=10)

    if bot_labels:
        if not notop:
            ax_b = ax.secondary_xaxis('bottom')
            ax_b.set_xlabel('retention persistence', fontsize=8)
            ax_b.set_xticks(np.arange(data.shape[1]))
            ax_b.set_xticklabels(bot_labels, rotation=90)
        else:
            ax.tick_params(
                top=False, bottom=True, labeltop=False, labelbottom=True,
                left=False, labelleft=False
            )
            ax_b = ax.secondary_xaxis('top')
            ax_b.set_xlabel(xaxis_top_label, fontsize=8)
            ax_b.tick_params(
                top=False, bottom=False, labeltop=False, labelbottom=False
            )

    # Turn spines off and create white grid.
    try:
        ax.spines[:].set_visible(False)
    except TypeError:
        pass
    if lgnd_info:
        plt.legend(
            handles=lgnd_info, prop={'size': 6},
            ncol=2,
            fancybox=True,
            loc='upper center',
            # for box on top:
            bbox_to_anchor=(0.45, 1.12)
        )
    ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
    ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
    ax.grid(which="minor", color="w", linestyle='-', linewidth=3)
    ax.tick_params(which="minor", bottom=False, left=False)
    if lgnd_info:
        return im
    else:
        return im, cbar


def create_intron_matrix(transcript_df, output_dir, writeplots=False):
    transcript_df.sort_values(
        by="read_introns", key=lambda x: x.str.len(),
        ascending=False, inplace=True
    )
    tx_introns = transcript_df['tx_introns'].unique()[0]
    tx = transcript_df['transcript'].unique()[0]
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

    if writeplots:
        fig, ax = plt.subplots()
        figfile = os.path.join(
            output_dir, 'intronheatmap1_{}_5to3.pdf'.format(tx)
        )
        heatmap(
            np.array(plot_list), reads, tx_introns, figfile, ax=ax,
            cmap="YlGn", cbarlabel="intron spliced/unspliced"
        )
        fig.tight_layout()

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


def calculate_persistences(intron_df):
    persistence_vals = []
    n_x = len(intron_df)
    for int_count, intron in enumerate(intron_df.columns.values):
        int_df = intron_df.dropna(subset=[intron]).copy()
        if len(int_df) == 0:
            persistence_vals.append('0')
        elif int_df[intron].sum() == int_df[intron].count():
            persistence_vals.append('0')
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
            persistence_vals.append(twodigit_string(persistence))
    return persistence_vals


def print_mutual_heatmaps(txdf1, txdf2, output_dir, now):
    # # Set a higher value for min Retained Intron read ratio if desired:
    txdf1 = txdf1.loc[txdf1['RI-read_ratio'] > 0].copy()
    txdf2 = txdf2.loc[txdf2['RI-read_ratio'] > 0].copy()

    new_txs1 = set(txdf1['transcript'].tolist())
    new_txs = set(txdf2['transcript'].tolist()).intersection(new_txs1)
    txdf1 = txdf1.loc[txdf1['transcript'].isin(new_txs)]
    txdf2 = txdf2.loc[txdf2['transcript'].isin(new_txs)]
    final_txs = set()
    for tx in new_txs:
        h_reads = len(txdf1.loc[txdf1['transcript'] == tx])
        i_reads = len(txdf2.loc[txdf2['transcript'] == tx])
        if h_reads > 20 or i_reads > 20:
            continue
        final_txs.add(tx)
    final_txs = list(final_txs)
    final_txs.sort()
    num_txs = 0
    for tx in final_txs:
        sub_dfh = txdf1.loc[txdf1['transcript'] == tx].copy()
        sub_dfi = txdf2.loc[txdf2['transcript'] == tx].copy()
        raw_dfh = create_intron_matrix(sub_dfh, output_dir)
        if len(raw_dfh.columns.values) > 20:
            continue
        else:
            num_txs += 1
        raw_dfi = create_intron_matrix(sub_dfi, output_dir)
        persistence_vals_h = calculate_persistences(raw_dfh)
        persistence_vals_i = calculate_persistences(raw_dfi)

        width = 6
        height = 6
        plt.rcParams['figure.figsize'] = (width, height)
        base_colormap = 'YlGn'
        top_scaled_cmap = cm.get_cmap(base_colormap, 256)
        newcolors = colors.ListedColormap(
            top_scaled_cmap(np.linspace(0.25, .9, 256))
        )
        splice_color = newcolors(255)
        ri_color = newcolors(0)

        legend_list = []
        for p_col, p_lab in zip([splice_color, ri_color],
                                ['spliced', 'unspliced']):
            legend_list.append(mpatches.Patch(color=p_col, label=p_lab))

        fig, (ax_l, ax_r) = plt.subplots(nrows=1, ncols=2)

        figfile = os.path.join(
            output_dir, 'persistence_hm_mean_withint{}_{}.pdf'.format(tx, now)
        )
        heatmap(
            np.array(raw_dfh), [], persistence_vals_h, ax=ax_l,
            bot_labels=persistence_vals_h,
            lgnd_info=legend_list,
            notop=True,
            cmap=newcolors,
            yaxis_label='HX1 reads',
            xaxis_top_label='{} introns'.format(tx),
            xaxis_bottom_label='HX1 IR persistence',
            cbarlabel="intron spliced/unspliced"
        )
        heatmap(
            np.array(raw_dfi), [], persistence_vals_i, ax=ax_r,
            bot_labels=persistence_vals_i,
            lgnd_info=legend_list,
            notop=True,
            cmap=newcolors,
            yaxis_label='iPSC reads',
            xaxis_top_label='{} introns'.format(tx),
            xaxis_bottom_label='iPSC IR persistence',
            cbarlabel="intron spliced/unspliced"
        )
        fig.tight_layout()
        fig = plt.gcf()
        fig.savefig(figfile)
        fig.clf()
        plt.clf()
        plt.close(fig)
    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Prints heatmaps of transcripts with <20 introns and >=5 '
                    'reads in iPSC and HX1.'
    )
    parser.add_argument(
        '--HX1-read-transcript-map', '-H',
        help='file with name "RI_txs_to_read_ids_final__[timestamp]__.tsv" '
             'generated previously by this pacbio_reads_to_transcripts.py for '
             'sample HX1, mapping reads to target best-match transcripts.'
    )
    parser.add_argument(
        '--iPSC-read-transcript-map', '-i',
        help='file with name "RI_txs_to_read_ids_final__[timestamp]__.tsv" '
             'generated previously by this pacbio_reads_to_transcripts.py for '
             'sample iPSC, mapping reads to target best-match transcripts.'
    )
    parser.add_argument(
        '--output-directory', '-o', default='./',
        help='Directory to write output figures and files.'
    )
    parser.add_argument(
        '--readcounts-HX1', '-x',
        help='reads_per_gene_and_transcript_HX1.tsv generated by '
             'pacbio_reads_to_transcripts.py'
    )
    parser.add_argument(
        '--readcounts-iPSC', '-p',
        help='reads_per_gene_and_transcript_iPSC.tsv generated by '
             'pacbio_reads_to_transcripts.py'
    )

    args = parser.parse_args()
    h_read_map = args.HX1_read_transcript_map
    i_read_map = args.iPSC_read_transcript_map
    output_dir = args.output_directory
    h_readcounts = args.readcounts_HX1
    i_readcounts = args.readcounts_iPSC

    now = datetime.now().strftime('%m-%d-%Y_%H.%M.%S')

    h_reads = pd.read_table(h_readcounts, sep='\t')
    i_reads = pd.read_table(i_readcounts, sep='\t')
    rds = 'reads_per_transcript'
    h_txs = set(h_reads.loc[h_reads[rds] > 4]['transcript'].tolist())
    i_txs = set(i_reads.loc[i_reads[rds] > 4]['transcript'].tolist())
    mutual_transcripts = h_txs.intersection(i_txs)

    print(
        '{} HX1 transcripts, {} iPSC, {} mutual'.format(
            len(h_txs), len(i_txs), len(mutual_transcripts)
        )
    )
    h_txdf = pd.read_table(h_read_map)
    i_txdf = pd.read_table(i_read_map)

    h_txdf = h_txdf.loc[h_txdf['transcript'].isin(mutual_transcripts)]
    i_txdf = i_txdf.loc[i_txdf['transcript'].isin(mutual_transcripts)]

    print_mutual_heatmaps(h_txdf, i_txdf, output_dir, now)

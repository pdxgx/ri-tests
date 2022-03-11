#!/usr/bin/env python3

"""
pacbio_reads_to_transcripts.py
Python 3 code for mapping aligned pacbio reads to annotated transcripts

SAMPLE RUN:
time python intronomer-paper/benchmarking_data/pacbio_reads_to_transcripts.py
-g files/gencode.v34.annotation.gtf
-a SRP065930_SAMN04251426.merged.aligned.sorted.bam
-i shortread_results -I shortread_results
-m RI_txs_to_read_ids_final_10-28-2021_19.43.18.tsv

Add -p for confidence-vs-intron position plotting
Add -r for intron-read heatmaps  (currently only working for relatively low
intron and read counts

Only needs:
    -g and -a (original alignments + gtf file)
  OR
    -m (initially processed file)

"""
import argparse
from datetime import datetime
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import os
import pandas as pd
from scipy.stats import mannwhitneyu, ttest_ind


_THRSH = 'thresh'
_LR_COUNT = 'longread_intron_count'
_IR_COUNT = 'all_iREAD_RIs'
_IR_OVERLAP = 'iREAD_RI_overlap'
_IR_TP = 'iREAD_true_positives'
_IR_FP = 'iREAD_false_positives'
_IR_FN = 'iREAD_false_negatives'
_IR_PREC = 'iREAD_precision'
_IR_REC = 'iREAD_recall'
_IR_INTS = 'iREAD_introns'
_INT_COUNT = 'all_IntEREst_RIs'
_INT_OVERLAP = 'IntEREst_RI_overlap'
_INT_TP = 'IntEREst_true_positives'
_INT_FP = 'IntEREst_false_positives'
_INT_FN = 'IntEREst_false_negatives'
_INT_PREC = 'IntEREst_precision'
_INT_REC = 'IntEREst_recall'
_INT_INTS = 'IntEREst_introns'
_SI_COUNT = 'all_superintronic_RIs'
_SI_OVERLAP = 'superintronic_RI_overlap'
_SI_TP = 'superintronic_true_positives'
_SI_FP = 'superintronic_false_positives'
_SI_FN = 'superintronic_false_negatives'
_SI_PREC = 'superintronic_precision'
_SI_REC = 'superintronic_recall'
_SI_INTS = 'superintronic_introns'
_KMA_COUNT = 'all_kma_RIs'
_KMA_OVERLAP = 'kma_RI_overlap'
_KMA_TP = 'kma_true_positives'
_KMA_FP = 'kma_false_positives'
_KMA_FN = 'kma_false_negatives'
_KMA_PREC = 'kma_precision'
_KMA_REC = 'kma_recall'
_KMA_INTS = 'kma_introns'
_IRFS_COUNT = 'all_IRFinder-S_RIs'
_IRFS_OVERLAP = 'IRFinder-S_RI_overlap'
_IRFS_TP = 'IRFinder-S_true_positives'
_IRFS_FP = 'IRFinder-S_false_positives'
_IRFS_FN = 'IRFinder-S_false_negatives'
_IRFS_PREC = 'IRFinder-S_precision'
_IRFS_REC = 'IRFinder-S_recall'
_IRFS_INTS = 'IRFinder-S_introns'

# FOR CONTINUOUS RI METRIC ANALYSIS
_IR_COL_CONT = 'iread_fpkm_allintron_lwm'
_INT_COL_CONT = 'interest_fpkm_allintron_lwm'
_SI_COL_CONT = 'superintronic_score_allintron_lwm'
_KMA_COL_CONT = 'kma_tpm_allintron_lwm'
_IRFS_COL_CONT = 'irfinders_irratio_allintron_lwm'

# FOR CONTINUOUS RI METRIC ANALYSIS
_IR_COL_BIN = 'iread_fpkm_filtintron_lwm'
_INT_COL_BIN = 'interest_fpkm_filtintron_lwm'
_SI_COL_BIN = 'superintronic_score_filtintron_lwm'
_KMA_COL_BIN = 'kma_tpm_filtintron_lwm'
_IRFS_COL_BIN = 'irfinders_irratio_filtintron_lwm'

_COUNT = 'count'
_OVERLAP = 'overlap'
_TP = 'true positive'
_FP = 'false positive'
_FN = 'false negative'
_PREC = 'precision'
_REC = 'recall'
_INTS = 'introns'
_TNAME = 'tool_name'
_BIN_COL = 'binary column'
_CONT_COL = 'continuous column'

_IR = 'iREAD'
_SI = 'superintronic'
_KMA = 'KMA'
_INT = 'IntEREst'
_IRFS = 'IRFinder-S'

_TOOLS = {
    _IR: {
        _COUNT: _IR_COUNT, _TP: _IR_TP, _FP: _IR_FP, _FN: _IR_FN,
        _PREC: _IR_PREC, _REC: _IR_REC, _INTS: _IR_INTS,
        _OVERLAP: _IR_OVERLAP,
        _BIN_COL: _IR_COL_BIN, _CONT_COL: _IR_COL_CONT
    },
    _INT: {
        _COUNT: _INT_COUNT, _TP: _INT_TP, _FP: _INT_FP, _FN: _INT_FN,
        _PREC: _INT_PREC, _REC: _INT_REC, _INTS: _INT_INTS,
        _OVERLAP: _INT_OVERLAP,
        _BIN_COL: _INT_COL_BIN, _CONT_COL: _INT_COL_CONT
    },
    _SI: {
        _COUNT: _SI_COUNT, _TP: _SI_TP, _FP: _SI_FP, _FN: _SI_FN,
        _PREC: _SI_PREC, _REC: _SI_REC, _INTS: _SI_INTS,
        _OVERLAP: _SI_OVERLAP,
        _BIN_COL: _SI_COL_BIN, _CONT_COL: _SI_COL_CONT
    },
    _KMA: {
        _COUNT: _KMA_COUNT, _TP: _KMA_TP, _FP: _KMA_FP, _FN: _KMA_FN,
        _PREC: _KMA_PREC, _REC: _KMA_REC, _INTS: _KMA_INTS,
        _OVERLAP: _KMA_OVERLAP,
        _BIN_COL: _KMA_COL_BIN, _CONT_COL: _KMA_COL_CONT
    },
    _IRFS: {
        _COUNT: _IRFS_COUNT, _TP: _IRFS_TP, _FP: _IRFS_FP, _FN: _IRFS_FN,
        _PREC: _IRFS_PREC, _REC: _IRFS_REC, _INTS: _IRFS_INTS,
        _OVERLAP: _IRFS_OVERLAP,
        _BIN_COL: _IRFS_COL_BIN, _CONT_COL: _IRFS_COL_CONT
    }
}
_TOOL_COLUMNS = [_IRFS, _SI, _IR, _KMA, _INT]

_WIDTH = 'width'
_POS = 'intron_position_in_tx'
_MOTIF = 'motif'
_READS = 'numreads.median'
_TOT_F = 'total_overlapping_features'
_MAX_F = 'max_features_per_base'
_PERBASE_F = '%_bases_overlapped'

_PERS = 'max_intron_persistence'


def grouped_boxplots_axis(ax, data_dict, plot_dict, logscale=False, y_label='',
                          percent=False, right_lim_shift=2, x_label='',
                          legend=False, showxticks=True, yrotation=0,
                          fliersize=4, ytick_pos=[], ytick_label=[],
                          axislabel_fontsize=12, colors='light',
                          ticklabel_fontsize=12, legend_fontsize=8):
    """

    data_dict: a dictionary with one entry per cluster with the keys being
        another dictionary containing the actual boxplot data (“data”) like
        this:
        grouped_data_dict[cancer_abbr] = {}
        grouped_data_dict[cancer_abbr]['data'] = [in_sra, non_sra]
        So you would have
        in_sra -> `precision` and non_sra -> `recall`, and each would be a
        list containing the relevant data point for each sample data.

    plot_dict: still a dictionary containing a “light” and “dark” color for
    each box in the group but nothing else.
        plot_info_dict = {}
        plot_info_dict['light colors'] = ['xkcd:tangerine', 'xkcd:cerulean']
        plot_info_dict['dark colors'] = ['xkcd:pumpkin', 'xkcd:ocean blue']
        (number of light and dark colors should correspond to the number of
        different types of things you’re plotting per cluster - so 2 here,
        since we have precision and recall.


    :param fig_file:
    :param fig_size:
    :param logscale:
    :param y_label:
    :param percent:
    :param right_lim_shift:
    :param x_label:
    :return:
    """
    light_cols = plot_dict['light colors']
    dark_cols = plot_dict['dark colors']
    num_boxes = len(light_cols)
    adjust_val = num_boxes + 1
    if colors == 'light':
        main_color = light_cols
    else:
        main_color=dark_cols

    bp_pos = list(range(1, num_boxes + 1))

    label_locs = []
    curr_label_loc = sum(bp_pos) / len(bp_pos)
    labels = []
    for abbr, values in data_dict.items():
        if showxticks:
            labels.append(abbr)
        data = values['data']

        boxes = plt.boxplot(
            data, positions=bp_pos, widths=0.6, patch_artist=True
        )
        box_elements = zip(boxes['fliers'], boxes['medians'], boxes['boxes'])
        for i, (fly, med, box) in enumerate(box_elements):
            plt.setp(
                fly, color=light_cols[i], markerfacecolor=light_cols[i],
                marker='.', markersize=fliersize, linestyle='none',
                linewidth=0.15, markeredgecolor=dark_cols[i]
            )
            plt.setp(
                med, color='black', linewidth=1.5
            )
            plt.setp(box, facecolor=main_color[i], edgecolor='black')

        bp_pos = [x + adjust_val for x in bp_pos]
        label_locs.append(curr_label_loc)
        curr_label_loc += adjust_val

    ax.set_yticklabels([])
    ax.xaxis.grid(False)
    ax.yaxis.grid(True)
    ax.set_xlim(left=0, right=curr_label_loc - right_lim_shift)
    ax.spines['bottom'].set_color('black')
    ax.spines['top'].set_color('black')
    ax.spines['right'].set_color('black')
    ax.spines['left'].set_color('black')

    if logscale:
        ax.set_yscale('log')

    plt.ylabel(y_label, fontsize=axislabel_fontsize)
    plt.xlabel(x_label, fontsize=axislabel_fontsize)
    ax.set_xticks(label_locs)
    ax.set_xticklabels(labels)
    if showxticks:
        plt.setp(
            ax.xaxis.get_majorticklabels(), rotation=0,
            fontsize=ticklabel_fontsize, color='black'
        )

    if legend:
        legend_list = []
        for p_col, p_lab in zip(main_color, plot_dict['row labels']):
            legend_list.append(mpatches.Patch(color=p_col, label=p_lab))
        plt.legend(handles=legend_list, prop={'size': legend_fontsize})

    if y_label:
        if ytick_label:
            plt.yticks(
                ytick_pos, ytick_label, fontsize=ticklabel_fontsize,
                rotation=yrotation
            )
        else:
            plt.setp(
                ax.yaxis.get_majorticklabels(), fontsize=ticklabel_fontsize,
                color='black', rotation=yrotation
            )
            if percent:
                ax.yaxis.set_major_formatter(
                    ticker.FuncFormatter(
                        lambda y, _: '{}%'.format('{:,g}'.format(100 * y))
                    )
                )
            elif not logscale:
                ax.yaxis.set_major_formatter(
                    ticker.FuncFormatter(lambda y, _: '{:,g}'.format(y))
                )
    if ytick_pos:
        plt.yticks(
            ytick_pos, ytick_label, fontsize=ticklabel_fontsize,
            rotation=yrotation
        )
    return


def calculate_width(intron_coords):
    left, right =intron_coords.split(':')[1].split('-')
    return 1 + int(right) - int(left)


def all_results_features_vs_truth(hx1_bin, hx1_cont, ipsc_bin, ipsc_cont,
                                  output_dir, now):
    cols_to_load = [_POS, 'intron', _PERS, _READS, _PERBASE_F]
    for col in _TOOL_COLUMNS:
        cols_to_load.append(_TOOLS[col][_TP])
        cols_to_load.append(_TOOLS[col][_FP])
        cols_to_load.append(_TOOLS[col][_FN])

    dataframes = {
        'HX1': {
            'bin': pd.read_table(hx1_bin, sep='\t', usecols=cols_to_load),
            'con': pd.read_table(hx1_cont, sep='\t', usecols=cols_to_load)
        },
        'iPSC': {
            'bin': pd.read_table(ipsc_bin, sep='\t', usecols=cols_to_load),
            'con': pd.read_table(ipsc_cont, sep='\t', usecols=cols_to_load)
        }
    }
    for dfs in dataframes.values():
        for df in dfs.values():
            df[_WIDTH] = df['intron'].apply(lambda x: calculate_width(x))

    fig, axs = plt.subplots(3, 5,)
    plt.subplots_adjust(wspace=0.12, hspace=0.1)
    fig.set_figheight(12)
    fig.set_figwidth(21)
    axislabel_size = 15
    ticklabel_size = 12
    legend_fontsize = 8
    title_size = 15
    light_cols = ['#d3494e', '#f29e8e', '#448ee4', '#b1d1fc']
    dark_cols = ['#a90308', '#cf524e', '#0a437a', '#448ee4']
    row_labels = [
        'HX1: all potential RIs', 'HX1: called RIs',
        'iPSC: all potential RIs', 'iPSC: called RIs',
    ]
    plot_info_dict = {
        'light colors': light_cols,
        'dark colors': dark_cols,
        'row labels': row_labels
    }
    data_entries = ['TP', 'FP', 'FN']
    for i, col in enumerate(_TOOL_COLUMNS):
        # top row: intron length
        plot_data_dict = {x: {'data': [[], [], [], []]} for x in data_entries}
        for j, dfs in enumerate(dataframes.values()):
            plot_data_dict['TP']['data'][0+(j*2)] = list(
                dfs['con'].loc[dfs['con'][_TOOLS[col][_TP]] == 1][_WIDTH]
            )
            plot_data_dict['TP']['data'][1+(j*2)] = list(
                dfs['bin'].loc[dfs['bin'][_TOOLS[col][_TP]] == 1][_WIDTH]
            )
            plot_data_dict['FP']['data'][0+(j*2)] = list(
                dfs['con'].loc[dfs['con'][_TOOLS[col][_FP]] == 1][_WIDTH]
            )
            plot_data_dict['FP']['data'][1+(j*2)] = list(
                dfs['bin'].loc[dfs['bin'][_TOOLS[col][_FP]] == 1][_WIDTH]
            )
            plot_data_dict['FN']['data'][0+(j*2)] = list(
                dfs['con'].loc[dfs['con'][_TOOLS[col][_FN]] == 1][_WIDTH]
            )
            plot_data_dict['FN']['data'][1+(j*2)] = list(
                dfs['bin'].loc[dfs['bin'][_TOOLS[col][_FN]] == 1][_WIDTH]
            )
        curr_ax = axs[0, i]
        plt.sca(curr_ax)
        if i == 0:
            ylabel = 'Intron length\n(# bases)'
            ytick_pos = []
            ytick_labels = []
        else:
            ylabel = ''
            ytick_pos = [100, 1000, 10000, 100000]
            ytick_labels = ['', '', '', '']
        if i == 4:
            print_legend = True
        else:
            print_legend = False
        grouped_boxplots_axis(
            curr_ax, plot_data_dict, plot_info_dict, y_label=ylabel,
            legend=print_legend, yrotation=90, logscale=True,
            ytick_pos=ytick_pos, ytick_label=ytick_labels,
            axislabel_fontsize=axislabel_size,
            ticklabel_fontsize=ticklabel_size, legend_fontsize=legend_fontsize
        )
        curr_ax.set_ylim([20, 400000])
        curr_ax.set_title(col, fontsize=title_size)

        # 2nd row: intron position
        plot_data_dict = {x: {'data': [[], [], [], []]} for x in data_entries}
        for j, dfs in enumerate(dataframes.values()):
            plot_data_dict['TP']['data'][0+(j*2)] = list(
                dfs['con'].loc[dfs['con'][_TOOLS[col][_TP]] == 1][_POS]
            )
            plot_data_dict['TP']['data'][1+(j*2)] = list(
                dfs['bin'].loc[dfs['bin'][_TOOLS[col][_TP]] == 1][_POS]
            )
            plot_data_dict['FP']['data'][0+(j*2)] = list(
                dfs['con'].loc[dfs['con'][_TOOLS[col][_FP]] == 1][_POS]
            )
            plot_data_dict['FP']['data'][1+(j*2)] = list(
                dfs['bin'].loc[dfs['bin'][_TOOLS[col][_FP]] == 1][_POS]
            )
            plot_data_dict['FN']['data'][0+(j*2)] = list(
                dfs['con'].loc[dfs['con'][_TOOLS[col][_FN]] == 1][_POS]
            )
            plot_data_dict['FN']['data'][1+(j*2)] = list(
                dfs['bin'].loc[dfs['bin'][_TOOLS[col][_FN]] == 1][_POS]
            )
        curr_ax = axs[1, i]
        plt.sca(curr_ax)
        if i == 0:
            ylabel = "Intron 5'->3'\nposition in transcript"
        else:
            ylabel = ''
        grouped_boxplots_axis(
            curr_ax, plot_data_dict, plot_info_dict, y_label=ylabel,
            yrotation=90, axislabel_fontsize=axislabel_size,
            ticklabel_fontsize=ticklabel_size, legend_fontsize=legend_fontsize
        )

        # 3rd row: % of bases with overlapping features
        plot_data_dict = {x: {'data': [[], [], [], []]} for x in data_entries}
        for j, dfs in enumerate(dataframes.values()):
            plot_data_dict['TP']['data'][0+(j*2)] = list(
                dfs['con'].loc[dfs['con'][_TOOLS[col][_TP]] == 1][_PERBASE_F]
            )
            plot_data_dict['TP']['data'][1+(j*2)] = list(
                dfs['bin'].loc[dfs['bin'][_TOOLS[col][_TP]] == 1][_PERBASE_F]
            )
            plot_data_dict['FP']['data'][0+(j*2)] = list(
                dfs['con'].loc[dfs['con'][_TOOLS[col][_FP]] == 1][_PERBASE_F]
            )
            plot_data_dict['FP']['data'][1+(j*2)] = list(
                dfs['bin'].loc[dfs['bin'][_TOOLS[col][_FP]] == 1][_PERBASE_F]
            )
            plot_data_dict['FN']['data'][0+(j*2)] = list(
                dfs['con'].loc[dfs['con'][_TOOLS[col][_FN]] == 1][_PERBASE_F]
            )
            plot_data_dict['FN']['data'][1+(j*2)] = list(
                dfs['bin'].loc[dfs['bin'][_TOOLS[col][_FN]] == 1][_PERBASE_F]
            )
        curr_ax = axs[2, i]
        plt.sca(curr_ax)
        if i == 0:
            ylabel = 'Bases with\noverlapping exons'
        else:
            ylabel = ''
        grouped_boxplots_axis(
            curr_ax, plot_data_dict, plot_info_dict, y_label=ylabel,
            percent=True, yrotation=90, axislabel_fontsize=axislabel_size,
            ticklabel_fontsize=ticklabel_size, legend_fontsize=legend_fontsize
        )

    # Save figure
    fig = plt.gcf()
    fig_file = os.path.join(
        output_dir, 'features_vs_truth_all_results_{}.pdf'.format(now)
    )
    fig.savefig(fig_file, bbox_inches='tight', pad_inches=.1)
    return


def calculate_stats(hx1_bin, ipsc_bin):
    cols_to_load = [
        _POS, 'intron', _PERS, _READS, _PERBASE_F,
        _IR_COL_BIN, _INT_COL_BIN, _SI_COL_BIN, _KMA_COL_BIN, _IRFS_COL_BIN
    ]
    dataframes = [
        pd.read_table(hx1_bin, sep='\t', usecols=cols_to_load),
        pd.read_table(ipsc_bin, sep='\t', usecols=cols_to_load)
    ]
    for df in dataframes:
        df[_WIDTH] = df['intron'].apply(lambda x: calculate_width(x))

    headers = ['HX1', 'iPSC']
    tools = [_IRFS, _SI, _IR, _KMA, _INT]
    all_data = [
        {'lengths': {}, 'pos': {}, 'overlap': {}},
        {'lengths': {}, 'pos': {}, 'overlap': {}}
    ]
    for j, df in enumerate(dataframes):
        length_data = [list(df.loc[df[_PERS] >= 0.1][_WIDTH])]
        pos_data = [list(df.loc[df[_PERS] >= 0.1][_POS])]
        percent_data = [list(df.loc[df[_PERS] >= 0.1][_PERBASE_F])]
        all_data[j]['pos']['pacbio'] = length_data[0]
        all_data[j]['lengths']['pacbio'] = pos_data[0]
        all_data[j]['overlap']['pacbio'] = percent_data[0]

        for i, tool in enumerate(tools, 1):
            col = _TOOLS[tool][_BIN_COL]
            length_data.append(list(df.loc[df[col] > 0][_WIDTH]))
            pos_data.append(list(df.loc[df[col] > 0][_POS]))
            percent_data.append(list(df.loc[df[col] > 0][_PERBASE_F]))
            all_data[j]['pos'][tool] = length_data[-1]
            all_data[j]['lengths'][tool] = pos_data[-1]
            all_data[j]['overlap'][tool] = percent_data[-1]

    for datadict, sample in zip(all_data, headers):
        print('\n\nfor sample {}'.format(sample))
        for datatype, toolvals in datadict.items():
            print('\nfor intron {} vs. pacbio'.format(datatype))
            for tool in tools:
                U, p = ttest_ind(
                    toolvals[tool], toolvals['pacbio'], equal_var=False
                )
                # U, p = mannwhitneyu(
                #     toolvals[tool], toolvals['pacbio'], alternative='two-sided'
                # )
                print('{}: statistic = {}, pval = {}'.format(tool, U, p))
    return


def features_vs_detection(hx1_bin, ipsc_bin, output_dir, now):
    cols_to_load = [
        _POS, 'intron', _PERS, _READS, _PERBASE_F,
        _IR_COL_BIN, _INT_COL_BIN, _SI_COL_BIN, _KMA_COL_BIN, _IRFS_COL_BIN
    ]
    dataframes = [
        pd.read_table(hx1_bin, sep='\t', usecols=cols_to_load),
        pd.read_table(ipsc_bin, sep='\t', usecols=cols_to_load)
    ]
    for df in dataframes:
        df[_WIDTH] = df['intron'].apply(lambda x: calculate_width(x))

    fig, axs = plt.subplots(3, 2)
    plt.subplots_adjust(wspace=0.1, hspace=0.1)
    fig.set_figheight(10)
    fig.set_figwidth(7)

    light_cols = [
        '#eeeeee', '#fcefee', '#fdf0cf', '#eff4e5', '#e8f4f7', '#efecf5'
    ]
    dark_cols = [
        '#797979', '#e1665d', '#f8b712', '#689404', '#33a2b7', '#745bad'
    ]
    yrotation = 0
    axislabel_size = 16
    ticklabel_size = 14
    title_size = 16
    blank_labels = ['', '', '', '', '', '']
    headers = ['HX1', 'iPSC']
    tool_cols = [
        _IRFS_COL_BIN, _SI_COL_BIN, _IR_COL_BIN, _KMA_COL_BIN, _INT_COL_BIN
    ]
    for j, df in enumerate(dataframes):
        length_data = [list(df.loc[df[_PERS] >= 0.1][_WIDTH])]
        pos_data = [list(df.loc[df[_PERS] >= 0.1][_POS])]
        percent_data = [list(df.loc[df[_PERS] >= 0.1][_PERBASE_F])]
        for i, col in enumerate(tool_cols, 1):
            length_data.append(list(df.loc[df[col] > 0][_WIDTH]))
            pos_data.append(list(df.loc[df[col] > 0][_POS]))
            percent_data.append(list(df.loc[df[col] > 0][_PERBASE_F]))

        # 1st row: intron length
        ax = axs[0, j]
        plt.sca(ax)
        boxes = plt.boxplot(length_data, patch_artist=True)
        plt.xticks([1, 2, 3, 4, 5, 6], blank_labels, fontsize=ticklabel_size)
        box_elements = zip(boxes['fliers'], boxes['medians'], boxes['boxes'])
        for i, (fly, med, box) in enumerate(box_elements):
            plt.setp(
                fly, color=light_cols[i], markerfacecolor=light_cols[i],
                marker='.', markersize=12, linestyle='none', linewidth=0.15,
                markeredgecolor=dark_cols[i]
            )
            plt.setp(med, color='black', linewidth=2)
            plt.setp(box, facecolor=dark_cols[i], edgecolor='black')
        ax.set_ylim([20, 400000])
        ax.set_yscale('log')
        if j == 0:
            ylabel = 'Intron length\n(# bases)'
            plt.yticks(rotation=yrotation, fontsize=ticklabel_size)
        else:
            ylabel = ''
            plt.yticks([100, 1000, 10000, 100000], ['', '', '', ''])

        plt.ylabel(ylabel, fontsize=axislabel_size)
        ax.set_title(headers[j], fontsize=title_size)
        ax.xaxis.grid(False)
        ax.yaxis.grid(True)

        # 2nd row: intron position
        ax = axs[1, j]
        plt.sca(ax)
        boxes = plt.boxplot(pos_data, patch_artist=True)
        box_elements = zip(boxes['fliers'], boxes['medians'], boxes['boxes'])
        for i, (fly, med, box) in enumerate(box_elements):
            plt.setp(
                fly, color=light_cols[i], markerfacecolor=light_cols[i],
                marker='.', markersize=12, linestyle='none', linewidth=0.15,
                markeredgecolor=dark_cols[i]
            )
            plt.setp(med, color='black', linewidth=2)
            plt.setp(box, facecolor=dark_cols[i], edgecolor='black')
        if j == 0:
            ylabel = "Intron transcript position"
            plt.yticks(
                [0, 0.25 ,0.5, 0.75, 1], ["5'", '', 'middle', '', "3'"],
                rotation=yrotation, fontsize=ticklabel_size
            )
        else:
            ylabel = ''
            plt.yticks([0, 0.25 ,0.5, 0.75, 1], ['', '', '', '', ''])
        plt.xticks([1, 2, 3, 4, 5, 6], blank_labels, fontsize=ticklabel_size)
        plt.ylabel(ylabel, fontsize=axislabel_size)
        ax.xaxis.grid(False)
        ax.yaxis.grid(True)

        # 3rd row: % of bases with overlapping features
        ax = axs[2, j]
        plt.sca(ax)
        boxes = plt.boxplot(percent_data, patch_artist=True)
        row_labels = [
            'PacBio', 'IRFinder-S', 'superintronic',
            'iREAD', 'KMA', 'IntEREst',
        ]
        plt.xticks(
            [1, 2, 3, 4, 5, 6], row_labels, fontsize=ticklabel_size,
            rotation=45, ha="right"
        )
        box_elements = zip(boxes['fliers'], boxes['medians'], boxes['boxes'])
        for i, (fly, med, box) in enumerate(box_elements):
            plt.setp(
                fly, color=light_cols[i], markerfacecolor=light_cols[i],
                marker='.', markersize=12, linestyle='none', linewidth=0.15,
                markeredgecolor=dark_cols[i]
            )
            plt.setp(med, color='black', linewidth=2)
            plt.setp(box, facecolor=dark_cols[i], edgecolor='black')
        if j == 0:
            ylabel = 'Bases with\noverlapping exons'
            plt.yticks(rotation=yrotation, fontsize=ticklabel_size)
            ax.yaxis.set_major_formatter(
                ticker.FuncFormatter(
                    lambda y, _: '{}%'.format('{:,g}'.format(100 * y))
                )
            )
        else:
            ylabel = ''
            plt.yticks([0, 0.2, 0.4, 0.6, 0.8, 1], ['', '', '', '', '', ''])
        plt.ylabel(ylabel, fontsize=axislabel_size)
        ax.xaxis.grid(False)
        ax.yaxis.grid(True)

    # Save figure
    fig = plt.gcf()
    fig_file = os.path.join(
        output_dir, 'features_vs_detection_all_{}.pdf'.format(now)
    )
    fig.savefig(fig_file, bbox_inches='tight', pad_inches=.1)
    return


def filtered_groupedbysamp_featvstruth(hx1_bin, ipsc_bin, output_dir, now):
    cols_to_load = [_POS, 'intron', _PERS, _READS, _PERBASE_F]
    for col in _TOOL_COLUMNS:
        cols_to_load.append(_TOOLS[col][_TP])
        cols_to_load.append(_TOOLS[col][_FP])
        cols_to_load.append(_TOOLS[col][_FN])

    dataframes = [
        pd.read_table(hx1_bin, sep='\t', usecols=cols_to_load),
        pd.read_table(ipsc_bin, sep='\t', usecols=cols_to_load)
    ]
    for df in dataframes:
        df[_WIDTH] = df['intron'].apply(lambda x: calculate_width(x))

    fig, axs = plt.subplots(3, 5)
    plt.subplots_adjust(wspace=0.1, hspace=0.1)
    fig.set_figheight(10)
    # fig.set_figwidth(17)
    fig.set_figwidth(15)

    light_cols = [ '#9db92c', '#d8788a', '#6d9cc6']
    dark_cols = ['#4e5c16', '#944a58', '#335c80']
    # row_labels = ['True positives', 'False positives', 'False negatives']
    row_labels = ['TPs', 'FPs', 'FNs']
    plot_info_dict = {
        'light colors': light_cols,
        'dark colors': dark_cols,
        'row labels': row_labels
    }
    data_entries = ['HX1', 'iPSC']
    yrotation = 0
    axislabel_size = 16
    ticklabel_size = 14
    legend_fontsize= 12
    title_size = 16
    for i, col in enumerate(_TOOL_COLUMNS):
        # top row: intron length
        plot_data_dict = {x: {'data': [[], [], []]} for x in data_entries}
        for j, df in enumerate(dataframes):
            if j == 0:
                name = 'HX1'
            else:
                name = 'iPSC'
            plot_data_dict[name]['data'][0] = list(
                df.loc[df[_TOOLS[col][_TP]] == 1][_WIDTH]
            )
            plot_data_dict[name]['data'][1] = list(
                df.loc[df[_TOOLS[col][_FP]] == 1][_WIDTH]
            )
            plot_data_dict[name]['data'][2] = list(
                df.loc[df[_TOOLS[col][_FN]] == 1][_WIDTH]
            )
        curr_ax = axs[0, i]
        plt.sca(curr_ax)
        if i == 0:
            ylabel = 'Intron length\n(# bases)'
            print_legend = True
            ytick_pos = []
            ytick_labels = []
        else:
            print_legend = False
            ylabel = ''
            ytick_pos = [100, 1000, 10000, 100000]
            ytick_labels = ['', '', '', '']

        grouped_boxplots_axis(
            curr_ax, plot_data_dict, plot_info_dict, y_label=ylabel,
            legend=print_legend, yrotation=yrotation, fliersize=9,
            showxticks=False, logscale=True, ytick_label=ytick_labels,
            ytick_pos=ytick_pos, axislabel_fontsize=axislabel_size,
            ticklabel_fontsize=ticklabel_size, legend_fontsize=legend_fontsize
        )
        curr_ax.set_ylim([20, 400000])
        curr_ax.set_title(col, fontsize=title_size)

        # 2nd row: intron position
        plot_data_dict = {x: {'data': [[], [], []]} for x in data_entries}
        for j, df in enumerate(dataframes):
            if j == 0:
                name = 'HX1'
            else:
                name = 'iPSC'
            plot_data_dict[name]['data'][0] = list(
                df.loc[df[_TOOLS[col][_TP]] == 1][_POS]
            )
            plot_data_dict[name]['data'][1] = list(
                df.loc[df[_TOOLS[col][_FP]] == 1][_POS]
            )
            plot_data_dict[name]['data'][2] = list(
                df.loc[df[_TOOLS[col][_FN]] == 1][_POS]
            )
        curr_ax = axs[1, i]
        plt.sca(curr_ax)
        if i == 0:
            ylabel = "Intron transcript position"
            ytick_pos = [0, 0.25, 0.5, 0.75, 1]
            ytick_label = ["5'", '', 'middle', '', "3'"]
        else:
            ylabel = ''
            ytick_pos = []
            ytick_label = [],
        grouped_boxplots_axis(
            curr_ax, plot_data_dict, plot_info_dict, y_label=ylabel,
            yrotation=yrotation, fliersize=9, ytick_pos=ytick_pos,
            ytick_label=ytick_label, showxticks=False,
            axislabel_fontsize=axislabel_size,
            ticklabel_fontsize=ticklabel_size, legend_fontsize=legend_fontsize
        )

        # 3rd row: % of bases with overlapping features
        plot_data_dict = {x: {'data': [[], [], []]} for x in data_entries}
        for j, df in enumerate(dataframes):
            if j == 0:
                name = 'HX1'
            else:
                name = 'iPSC'
            plot_data_dict[name]['data'][0] = list(
                df.loc[df[_TOOLS[col][_TP]] == 1][_PERBASE_F]
            )
            plot_data_dict[name]['data'][1] = list(
                df.loc[df[_TOOLS[col][_FP]] == 1][_PERBASE_F]
            )
            plot_data_dict[name]['data'][2] = list(
                df.loc[df[_TOOLS[col][_FN]] == 1][_PERBASE_F]
            )
        curr_ax = axs[2, i]
        plt.sca(curr_ax)
        if i == 0:
            ylabel = 'Bases with\noverlapping exons'
        else:
            ylabel = ''
        grouped_boxplots_axis(
            curr_ax, plot_data_dict, plot_info_dict, y_label=ylabel,
            percent=True, yrotation=yrotation, fliersize=9,
            axislabel_fontsize=axislabel_size,
            ticklabel_fontsize=ticklabel_size, legend_fontsize=legend_fontsize
        )

    # Save figure
    fig = plt.gcf()
    fig_file = os.path.join(
        output_dir, 'features_vs_truth_samplegrouped_{}.pdf'.format(now)
    )
    fig.savefig(fig_file, bbox_inches='tight', pad_inches=.1)
    return


def shared_FP_featvstruth(hx1_bin, ipsc_bin, output_dir, now):
    cols_to_load = [_POS, 'intron', _PERS, _READS, _PERBASE_F]
    for col in _TOOL_COLUMNS:
        cols_to_load.append(_TOOLS[col][_FP])

    dataframes_init = [
        pd.read_table(hx1_bin, sep='\t', usecols=cols_to_load),
        pd.read_table(ipsc_bin, sep='\t', usecols=cols_to_load)
    ]

    dataframes = []
    dataframes.append(dataframes_init[0].loc[
        (dataframes_init[0][_IRFS_FP] == 1) &
        (dataframes_init[0][_IR_FP] == 1) &
        (dataframes_init[0][_KMA_FP] == 1) &
        (dataframes_init[0][_SI_FP] == 1) &
        (dataframes_init[0][_INT_FP] == 1)
    ].copy())

    dataframes.append(dataframes_init[1].loc[
        (dataframes_init[1][_IRFS_FP] == 1) &
        (dataframes_init[1][_IR_FP] == 1) &
        (dataframes_init[1][_KMA_FP] == 1) &
        (dataframes_init[1][_SI_FP] == 1) &
        (dataframes_init[1][_INT_FP] == 1)
    ].copy())

    for df in dataframes:
        df[_WIDTH] = df['intron'].apply(lambda x: calculate_width(x))

    width_list = []
    pos = []
    overlap = []
    for i, df in enumerate(dataframes):
        width_list.append(df[_WIDTH].tolist())
        pos.append(df[_POS].tolist())
        overlap.append(df[_PERBASE_F].tolist())

    fig, axs = plt.subplots(1, 3)
    fig.set_figheight(3)
    fig.set_figwidth(10)
    plt.subplots_adjust(wspace=0.3, hspace=0.1)
    light_cols = ['#f29e8e', '#b1d1fc']
    dark_cols = ['#cf524e', '#448ee4']
    box_labels = ['HX1 shared FPs', 'iPSC shared FPs']
    yrotation = 90
    axislabel_size = 11
    ticklabel_size = 9

    # Left col: length
    ax = axs[0]
    plt.sca(ax)
    boxes = plt.boxplot(width_list, patch_artist=True)
    plt.xticks([1, 2], box_labels, fontsize=ticklabel_size)
    box_elements = zip(boxes['fliers'], boxes['medians'], boxes['boxes'])
    for i, (fly, med, box) in enumerate(box_elements):
        plt.setp(
            fly, color=light_cols[i], markerfacecolor=light_cols[i],
            marker='.', markersize=12, linestyle='none', linewidth=0.15,
            markeredgecolor=dark_cols[i]
        )
        plt.setp(med, color='black', linewidth=2)
        plt.setp(box, facecolor=light_cols[i], edgecolor='black')

    plt.setp(ax.yaxis.get_majorticklabels(), rotation=90)
    ylabel = 'Intron length (# bases)'
    plt.ylabel(ylabel, fontsize=axislabel_size)
    plt.yticks(rotation=yrotation, fontsize=ticklabel_size)
    # ax.set_ylim([0, 5700])
    ax.set_yscale('log')
    ax.xaxis.grid(False)
    ax.yaxis.grid(True)

    # 2nd col: intron position
    ax = axs[1]
    plt.sca(ax)
    boxes = plt.boxplot(pos, patch_artist=True)
    plt.xticks([1, 2], box_labels, fontsize=ticklabel_size)
    box_elements = zip(boxes['fliers'], boxes['medians'], boxes['boxes'])
    for i, (fly, med, box) in enumerate(box_elements):
        plt.setp(
            fly, color=light_cols[i], markerfacecolor=light_cols[i],
            marker='.', markersize=12, linestyle='none', linewidth=0.15,
            markeredgecolor=dark_cols[i]
        )
        plt.setp(med, color='black', linewidth=2)
        plt.setp(box, facecolor=light_cols[i], edgecolor='black')
    ylabel = "Intron transcript position"
    plt.yticks(
        [0, 0.25, 0.5, 0.75, 1], ["5'", '', 'middle', '', "3'"],
        rotation=yrotation, fontsize=ticklabel_size
    )
    plt.ylabel(ylabel, fontsize=axislabel_size)
    ax.xaxis.grid(False)
    ax.yaxis.grid(True)

    # 3rd col: % of bases with overlapping features
    ax = axs[2]
    plt.sca(ax)
    boxes = plt.boxplot(overlap, patch_artist=True)
    plt.xticks([1, 2], box_labels, fontsize=ticklabel_size)
    box_elements = zip(boxes['fliers'], boxes['medians'], boxes['boxes'])
    for i, (fly, med, box) in enumerate(box_elements):
        plt.setp(
            fly, color=light_cols[i], markerfacecolor=light_cols[i],
            marker='.', markersize=12, linestyle='none', linewidth=0.15,
            markeredgecolor=dark_cols[i]
        )
        plt.setp(med, color='black', linewidth=2)
        plt.setp(box, facecolor=light_cols[i], edgecolor='black')
    plt.setp(
        ax.yaxis.get_majorticklabels(), rotation=90
    )
    plt.yticks(rotation=yrotation, fontsize=ticklabel_size)
    ax.yaxis.set_major_formatter(
        ticker.FuncFormatter(
            lambda y, _: '{}%'.format('{:,g}'.format(100 * y))
        )
    )
    ylabel = 'Bases with\noverlapping features'
    plt.ylabel(ylabel, fontsize=axislabel_size)
    ax.xaxis.grid(False)
    ax.yaxis.grid(True)

    # Save figure
    fig = plt.gcf()
    fig_file = os.path.join(
        output_dir, 'features_vs_truth_sharedFPs_{}.pdf'.format(now)
    )
    fig.savefig(fig_file, bbox_inches='tight', pad_inches=.1)
    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Extracts retained-intron exons from stringtie-assembled '
                    'transcripts.'
    )
    parser.add_argument(
        '--HX1-called-results', '-c',
        help='Truth category assignment & intron features summary file for '
             'HX1, binary RI calls'
    )
    parser.add_argument(
        '--HX1-nonzero-results', '-n',
        help='Truth category assignment & intron features summary file for '
             'HX1, continuous RI calls'
    )
    parser.add_argument(
        '--iPSC-called-results', '-C',
        help='Truth category assignment & intron features summary file for '
             'iPSC, binary RI calls'
    )
    parser.add_argument(
        '--iPSC-nonzero-results', '-N',
        help='Truth category assignment & intron features summary file for '
             'iPSC, continuous RI calls'
    )
    parser.add_argument(
        '--output-directory', '-o',
        help='directory in which to store output files and figures; if no '
             'directory specified, this will be set as the same directory '
             'containing the range summarized results.'
    )

    args = parser.parse_args()
    hx1_call = args.HX1_called_results
    hx1_expr = args.HX1_nonzero_results
    ipsc_call = args.iPSC_called_results
    ipsc_expr = args.iPSC_nonzero_results
    out_dir = args.output_directory

    now = datetime.now().strftime('%m-%d-%Y_%H.%M.%S')

    # calculate_stats(hx1_call, ipsc_call)
    features_vs_detection(hx1_call, ipsc_call, out_dir, now)
    all_results_features_vs_truth(
        hx1_call, hx1_expr, ipsc_call, ipsc_expr, out_dir, now
    )
    filtered_groupedbysamp_featvstruth(hx1_call, ipsc_call, out_dir, now)
    shared_FP_featvstruth(hx1_call, ipsc_call, out_dir, now)

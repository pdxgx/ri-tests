#!/usr/bin/env python3

"""
Figs2B_4A_S12_S13_intron_features.py
Python 3 code for plotting intron feature boxplot sets

SAMPLE RUN:
time python ri-tests/results/figures/features_vs_truth_categories.py
-c final_8tool_HX1/called_RIs/
    called_RI_data_summary_HX1featureannotated_GCcontent.tsv
-n final_8tool_HX1/nonzero_RIs/
    nonzero_RI_data_summary_HX1featureannotated_GCcontent.tsv
-C final_8tool_iPSC/called_RIs/
    called_RI_data_summary_iPSCfeatureannotated_GCcontent.tsv
-N final_8tool_iPSC/nonzero_RIs/
    nonzero_RI_data_summary_iPSCfeatureannotated_GCcontent.tsv
-o paper_results/

"""
import argparse
from datetime import datetime
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import os
import pandas as pd
from scipy.stats import ttest_ind

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
_IR_F = 'iREAD_fscore'

_INT_COUNT = 'all_IntEREst_RIs'
_INT_OVERLAP = 'IntEREst_RI_overlap'
_INT_TP = 'IntEREst_true_positives'
_INT_FP = 'IntEREst_false_positives'
_INT_FN = 'IntEREst_false_negatives'
_INT_PREC = 'IntEREst_precision'
_INT_REC = 'IntEREst_recall'
_INT_INTS = 'IntEREst_introns'
_INT_F = 'IntEREst_fscore'

_SI_COUNT = 'all_superintronic_RIs'
_SI_OVERLAP = 'superintronic_RI_overlap'
_SI_TP = 'superintronic_true_positives'
_SI_FP = 'superintronic_false_positives'
_SI_FN = 'superintronic_false_negatives'
_SI_PREC = 'superintronic_precision'
_SI_REC = 'superintronic_recall'
_SI_INTS = 'superintronic_introns'
_SI_F = 'superintronic_fscore'

_KMA_COUNT = 'all_kma_RIs'
_KMA_OVERLAP = 'kma_RI_overlap'
_KMA_TP = 'kma_true_positives'
_KMA_FP = 'kma_false_positives'
_KMA_FN = 'kma_false_negatives'
_KMA_PREC = 'kma_precision'
_KMA_REC = 'kma_recall'
_KMA_INTS = 'kma_introns'
_KMA_F = 'KMA_fscore'

_IRFS_COUNT = 'all_IRFinder-S_RIs'
_IRFS_OVERLAP = 'IRFinder-S_RI_overlap'
_IRFS_TP = 'IRFinder-S_true_positives'
_IRFS_FP = 'IRFinder-S_false_positives'
_IRFS_FN = 'IRFinder-S_false_negatives'
_IRFS_PREC = 'IRFinder-S_precision'
_IRFS_REC = 'IRFinder-S_recall'
_IRFS_INTS = 'IRFinder-S_introns'
_IRFS_F = 'IRFinder-S_fscore'

_SUP_COUNT = 'all_SUPPA2_RIs'
_SUP_OVERLAP = 'SUPPA2_RI_overlap'
_SUP_TP = 'SUPPA2_true_positives'
_SUP_FP = 'SUPPA2_false_positives'
_SUP_FN = 'SUPPA2_false_negatives'
_SUP_PREC = 'SUPPA2_precision'
_SUP_REC = 'SUPPA2_recall'
_SUP_INTS = 'SUPPA2_introns'
_SUP_F = 'SUPPA2_fscore'

_MAJ_COUNT = 'all_MAJIQ_RIs'
_MAJ_OVERLAP = 'MAJIQ_RI_overlap'
_MAJ_TP = 'MAJIQ_true_positives'
_MAJ_FP = 'MAJIQ_false_positives'
_MAJ_FN = 'MAJIQ_false_negatives'
_MAJ_PREC = 'MAJIQ_precision'
_MAJ_REC = 'MAJIQ_recall'
_MAJ_INTS = 'MAJIQ_introns'
_MAJ_F = 'MAJIQ_fscore'

_RMA_COUNT = 'all_rMATS_RIs'
_RMA_OVERLAP = 'rMATS_RI_overlap'
_RMA_TP = 'rMATS_true_positives'
_RMA_FP = 'rMATS_false_positives'
_RMA_FN = 'rMATS_false_negatives'
_RMA_PREC = 'rMATS_precision'
_RMA_REC = 'rMATS_recall'
_RMA_INTS = 'rMATS_introns'
_RMA_F = 'rMATS_fscore'

# FOR CONTINUOUS RI METRIC ANALYSIS
_IR_COL_CONT = 'iread_fpkm_allintron_lwm'
_INT_COL_CONT = 'interest_fpkm_allintron_lwm'
_SI_COL_CONT = 'superintronic_score_allintron_lwm'
_KMA_COL_CONT = 'kma_tpm_allintron_lwm'
_IRFS_COL_CONT = 'irfinders_irratio_allintron_lwm'
_SUP_COL_CONT = 'suppa2_psi_allintron_lwm'
_MAJ_COL_CONT = 'majiq_meanpsi_allintron_lwm'
_RMA_COL_CONT = 'rmats_inclvl_allintron_lwm'


# FOR CONTINUOUS RI METRIC ANALYSIS
_IR_COL_BIN = 'iread_fpkm_filtintron_lwm'
_INT_COL_BIN = 'interest_fpkm_filtintron_lwm'
_SI_COL_BIN = 'superintronic_score_filtintron_lwm'
_KMA_COL_BIN = 'kma_tpm_filtintron_lwm'
_IRFS_COL_BIN = 'irfinders_irratio_filtintron_lwm'
_SUP_COL_BIN = 'suppa2_psi_filtintron_lwm'
_MAJ_COL_BIN = 'majiq_meanpsi_filtintron_lwm'
_RMA_COL_BIN = 'rmats_inclvl_filtintron_lwm'

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
_FSCORE = 'f-score'

_IR = 'iREAD'
_SI = 'superintronic'
_KMA = 'KMA'
_INT = 'IntEREst'
_IRFS = 'IRFinder-S'
_SUP = 'SUPPA2'
_RMA = 'rMATS'
_MAJ = 'MAJIQ'

_TOOLS = {
    _IR: {
        _COUNT: _IR_COUNT, _TP: _IR_TP, _FP: _IR_FP, _FN: _IR_FN,
        _PREC: _IR_PREC, _REC: _IR_REC, _INTS: _IR_INTS,
        _OVERLAP: _IR_OVERLAP, _FSCORE: _IR_F,
        _BIN_COL: _IR_COL_BIN, _CONT_COL: _IR_COL_CONT
    },
    _INT: {
        _COUNT: _INT_COUNT, _TP: _INT_TP, _FP: _INT_FP, _FN: _INT_FN,
        _PREC: _INT_PREC, _REC: _INT_REC, _INTS: _INT_INTS,
        _OVERLAP: _INT_OVERLAP, _FSCORE: _INT_F,
        _BIN_COL: _INT_COL_BIN, _CONT_COL: _INT_COL_CONT
    },
    _SI: {
        _COUNT: _SI_COUNT, _TP: _SI_TP, _FP: _SI_FP, _FN: _SI_FN,
        _PREC: _SI_PREC, _REC: _SI_REC, _INTS: _SI_INTS,
        _OVERLAP: _SI_OVERLAP, _FSCORE: _SI_F,
        _BIN_COL: _SI_COL_BIN, _CONT_COL: _SI_COL_CONT
    },
    _KMA: {
        _COUNT: _KMA_COUNT, _TP: _KMA_TP, _FP: _KMA_FP, _FN: _KMA_FN,
        _PREC: _KMA_PREC, _REC: _KMA_REC, _INTS: _KMA_INTS,
        _OVERLAP: _KMA_OVERLAP, _FSCORE: _KMA_F,
        _BIN_COL: _KMA_COL_BIN, _CONT_COL: _KMA_COL_CONT
    },
    _IRFS: {
        _COUNT: _IRFS_COUNT, _TP: _IRFS_TP, _FP: _IRFS_FP, _FN: _IRFS_FN,
        _PREC: _IRFS_PREC, _REC: _IRFS_REC, _INTS: _IRFS_INTS,
        _OVERLAP: _IRFS_OVERLAP, _FSCORE: _IRFS_F,
        _BIN_COL: _IRFS_COL_BIN, _CONT_COL: _IRFS_COL_CONT
    },
    _SUP: {
        _COUNT: _SUP_COUNT, _TP: _SUP_TP, _FP: _SUP_FP, _FN: _SUP_FN,
        _PREC: _SUP_PREC, _REC: _SUP_REC, _INTS: _SUP_INTS,
        _OVERLAP: _SUP_OVERLAP, _FSCORE: _SUP_F,
        _BIN_COL: _SUP_COL_BIN, _CONT_COL: _SUP_COL_CONT
    },
    _MAJ: {
        _COUNT: _MAJ_COUNT, _TP: _MAJ_TP, _FP: _MAJ_FP, _FN: _MAJ_FN,
        _PREC: _MAJ_PREC, _REC: _MAJ_REC, _INTS: _MAJ_INTS,
        _OVERLAP: _MAJ_OVERLAP, _FSCORE: _MAJ_F,
        _BIN_COL: _MAJ_COL_BIN, _CONT_COL: _MAJ_COL_CONT
    },
    _RMA: {
        _COUNT: _RMA_COUNT, _TP: _RMA_TP, _FP: _RMA_FP, _FN: _RMA_FN,
        _PREC: _RMA_PREC, _REC: _RMA_REC, _INTS: _RMA_INTS,
        _OVERLAP: _RMA_OVERLAP, _FSCORE: _RMA_F,
        _BIN_COL: _RMA_COL_BIN, _CONT_COL: _RMA_COL_CONT
    }
}
_TOOL_COLUMNS = [_IRFS, _SI, _IR, _KMA, _INT, _RMA, _MAJ, _SUP]

_WIDTH = 'width'
_POS = 'intron_position_in_tx'
_MOTIF = 'motif'
_READS = 'numreads.median'
_TOT_F = 'total_overlapping_features'
_MAX_F = 'max_features_per_base'
_PERBASE_F = '%_bases_overlapped'
_GC_PERC = 'gc_fract'
_PERS = 'max_intron_persistence'

_LIGHT_COLOR_MAP = {
    _MAJ_COL_BIN: '#ffecec', _RMA_COL_BIN: '#f7fded',
    _SUP_COL_BIN: '#eef8ff', _IRFS_COL_BIN: '#fcefee',
    _SI_COL_BIN: '#fdf0cf', _IR_COL_BIN: '#eff4e5',
    _KMA_COL_BIN: '#e8f4f7', _INT_COL_BIN: '#efecf5', _PERS: '#eeeeee'
}
_DARK_COLOR_MAP = {
    _MAJ_COL_BIN: '#FFC0C0', _RMA_COL_BIN: '#DAF7A6',
    _SUP_COL_BIN: '#ABDDFF', _IRFS_COL_BIN: '#e1665d',
    _SI_COL_BIN: '#f8b712', _IR_COL_BIN: '#689404',
    _KMA_COL_BIN: '#33a2b7', _INT_COL_BIN: '#745bad', _PERS: '#797979'
}

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
    cols_to_load = [_POS, 'intron', _PERS, _READS, _PERBASE_F, _GC_PERC]
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

    fig, axs = plt.subplots(4, 8)
    plt.subplots_adjust(wspace=0.12, hspace=0.1)
    fig.set_figheight(18)
    fig.set_figwidth(27)
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
        if i == 7:
            print_legend = True
        else:
            print_legend = False
        grouped_boxplots_axis(
            curr_ax, plot_data_dict, plot_info_dict, y_label=ylabel,
            legend=print_legend, yrotation=90, logscale=True,
            ytick_pos=ytick_pos, ytick_label=ytick_labels, fliersize=8,
            axislabel_fontsize=axislabel_size,
            ticklabel_fontsize=ticklabel_size, legend_fontsize=legend_fontsize
        )
        curr_ax.set_ylim([20, 400000])
        curr_ax.set_title(col, fontsize=title_size)

        # 2nd row: intron GC content
        plot_data_dict = {x: {'data': [[], [], [], []]} for x in data_entries}
        for j, dfs in enumerate(dataframes.values()):
            plot_data_dict['TP']['data'][0 + (j * 2)] = list(
                dfs['con'].loc[dfs['con'][_TOOLS[col][_TP]] == 1][_GC_PERC]
            )
            plot_data_dict['TP']['data'][1 + (j * 2)] = list(
                dfs['bin'].loc[dfs['bin'][_TOOLS[col][_TP]] == 1][_GC_PERC]
            )
            plot_data_dict['FP']['data'][0 + (j * 2)] = list(
                dfs['con'].loc[dfs['con'][_TOOLS[col][_FP]] == 1][_GC_PERC]
            )
            plot_data_dict['FP']['data'][1 + (j * 2)] = list(
                dfs['bin'].loc[dfs['bin'][_TOOLS[col][_FP]] == 1][_GC_PERC]
            )
            plot_data_dict['FN']['data'][0 + (j * 2)] = list(
                dfs['con'].loc[dfs['con'][_TOOLS[col][_FN]] == 1][_GC_PERC]
            )
            plot_data_dict['FN']['data'][1 + (j * 2)] = list(
                dfs['bin'].loc[dfs['bin'][_TOOLS[col][_FN]] == 1][_GC_PERC]
            )
        curr_ax = axs[1, i]
        plt.sca(curr_ax)
        if i == 0:
            ylabel = "Intron GC content\n"
            ytick_pos = [0.1, 0.3, 0.5, 0.7, 0.9]
            ytick_labels = ['10%', '30%', '50%', '70%', '90%']
        else:
            ylabel = ''
            ytick_pos = [0.1, 0.3, 0.5, 0.7, 0.9]
            ytick_labels = ['', '', '', '', '']
        grouped_boxplots_axis(
            curr_ax, plot_data_dict, plot_info_dict, y_label=ylabel,
            percent=True, yrotation=90, axislabel_fontsize=axislabel_size,
            ytick_pos=ytick_pos, ytick_label=ytick_labels, fliersize=8,
            ticklabel_fontsize=ticklabel_size, legend_fontsize=legend_fontsize
        )
        curr_ax.set_ylim([0.075, 0.925])

        # 3rd row: intron position
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
        curr_ax = axs[2, i]
        plt.sca(curr_ax)
        if i == 0:
            ylabel = "Intron 5'->3'\nposition in transcript"
        else:
            ylabel = ''
        grouped_boxplots_axis(
            curr_ax, plot_data_dict, plot_info_dict, y_label=ylabel,
            yrotation=90, axislabel_fontsize=axislabel_size, fliersize=8,
            ticklabel_fontsize=ticklabel_size, legend_fontsize=legend_fontsize
        )

        # 4th row: % of bases with overlapping features
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
        curr_ax = axs[3, i]
        plt.sca(curr_ax)
        if i == 0:
            ylabel = 'Bases with\noverlapping exons'
        else:
            ylabel = ''
        grouped_boxplots_axis(
            curr_ax, plot_data_dict, plot_info_dict, y_label=ylabel,
            percent=True, yrotation=90, axislabel_fontsize=axislabel_size,
            ticklabel_fontsize=ticklabel_size, legend_fontsize=legend_fontsize,
            fliersize=8,
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
        _POS, 'intron', _PERS, _READS, _PERBASE_F, _GC_PERC,
        _IR_COL_BIN, _INT_COL_BIN, _SI_COL_BIN, _KMA_COL_BIN, _IRFS_COL_BIN,
        _RMA_COL_BIN, _MAJ_COL_BIN, _SUP_COL_BIN
    ]
    dataframes = [
        pd.read_table(hx1_bin, sep='\t', usecols=cols_to_load),
        pd.read_table(ipsc_bin, sep='\t', usecols=cols_to_load)
    ]
    for df in dataframes:
        df[_WIDTH] = df['intron'].apply(lambda x: calculate_width(x))

    fig, axs = plt.subplots(4, 2)
    plt.subplots_adjust(wspace=0.1, hspace=0.1)
    # fig.set_figheight(13)
    # fig.set_figwidth(8)
    fig.set_figheight(13)
    fig.set_figwidth(9)

    yrotation = 0
    axislabel_size = 16
    ticklabel_size = 14
    title_size = 16
    blank_labels = ['', '', '', '', '', '', '', '', '']
    headers = ['HX1', 'iPSC']
    tool_cols = [
        _MAJ_COL_BIN, _RMA_COL_BIN, _SUP_COL_BIN, _IRFS_COL_BIN, _SI_COL_BIN,
        _IR_COL_BIN, _KMA_COL_BIN, _INT_COL_BIN,
    ]
    light_cols = [_LIGHT_COLOR_MAP[col] for col in [_PERS] + tool_cols]
    dark_cols = [_DARK_COLOR_MAP[col] for col in [_PERS] + tool_cols]
    for j, df in enumerate(dataframes):
        length_data = [list(df.loc[df[_PERS] >= 0.1][_WIDTH])]
        pos_data = [list(df.loc[df[_PERS] >= 0.1][_POS])]
        percent_data = [list(df.loc[df[_PERS] >= 0.1][_PERBASE_F])]
        gc_data = [list(df.loc[df[_PERS] >= 0.1][_GC_PERC])]
        for i, col in enumerate(tool_cols, 1):
            length_data.append(list(df.loc[df[col] > 0][_WIDTH]))
            pos_data.append(list(df.loc[df[col] > 0][_POS]))
            percent_data.append(list(df.loc[df[col] > 0][_PERBASE_F]))
            gc_data.append(list(df.loc[df[col] > 0][_GC_PERC]))

        # 1st row: intron length
        ax = axs[0, j]
        plt.sca(ax)
        boxes = plt.boxplot(length_data, patch_artist=True)
        plt.xticks(
            [1, 2, 3, 4, 5, 6, 7, 8, 9], blank_labels, fontsize=ticklabel_size
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

        # 2nd row: intron GC content
        ax = axs[1, j]
        plt.sca(ax)
        boxes = plt.boxplot(gc_data, patch_artist=True)
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
            ylabel = "Intron GC content\n"
            plt.yticks(
                [0, 0.2, 0.4, 0.6, 0.8, 1], rotation=yrotation,
                fontsize=ticklabel_size
            )
            ax.yaxis.set_major_formatter(
                ticker.FuncFormatter(
                    lambda y, _: '{}%'.format('{:,g}'.format(100 * y))
                )
            )
        else:
            ylabel = ''
            plt.yticks([0, 0.2, 0.4, 0.6, 0.8, 1],
                       ['', '', '', '', '', ''])
        plt.xticks(
            [1, 2, 3, 4, 5, 6, 7, 8, 9], blank_labels, fontsize=ticklabel_size
        )
        plt.ylabel(ylabel, fontsize=axislabel_size)
        ax.xaxis.grid(False)
        ax.yaxis.grid(True)

        # 3rd row: intron position
        ax = axs[2, j]
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
        plt.xticks(
            [1, 2, 3, 4, 5, 6, 7, 8, 9], blank_labels, fontsize=ticklabel_size
        )
        plt.ylabel(ylabel, fontsize=axislabel_size)
        ax.xaxis.grid(False)
        ax.yaxis.grid(True)

        # 4th row: % of bases with overlapping features
        ax = axs[3, j]
        plt.sca(ax)
        boxes = plt.boxplot(percent_data, patch_artist=True)
        row_labels = [
            'PacBio', _MAJ, _RMA,  _SUP, 'IRFinder-S', 'superintronic',
            'iREAD', 'KMA', 'IntEREst',
        ]
        plt.xticks(
            [1, 2, 3, 4, 5, 6, 7, 8, 9], row_labels, fontsize=ticklabel_size,
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
    cols_to_load = [_POS, 'intron', _PERS, _READS, _PERBASE_F, _GC_PERC]
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

    fig, axs = plt.subplots(4, 8)
    plt.subplots_adjust(wspace=0.1, hspace=0.1)
    fig.set_figheight(13)
    # fig.set_figwidth(17)
    fig.set_figwidth(18)

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

        # 2nd row: intron GC content
        plot_data_dict = {x: {'data': [[], [], []]} for x in data_entries}
        for j, df in enumerate(dataframes):
            if j == 0:
                name = 'HX1'
            else:
                name = 'iPSC'
            plot_data_dict[name]['data'][0] = list(
                df.loc[df[_TOOLS[col][_TP]] == 1][_GC_PERC]
            )
            plot_data_dict[name]['data'][1] = list(
                df.loc[df[_TOOLS[col][_FP]] == 1][_GC_PERC]
            )
            plot_data_dict[name]['data'][2] = list(
                df.loc[df[_TOOLS[col][_FN]] == 1][_GC_PERC]
            )
        curr_ax = axs[1, i]
        plt.sca(curr_ax)
        if i == 0:
            ylabel = "Intron GC content\n"
        else:
            ylabel = ''
        grouped_boxplots_axis(
            curr_ax, plot_data_dict, plot_info_dict, y_label=ylabel,
            percent=True, yrotation=yrotation, fliersize=9,
            axislabel_fontsize=axislabel_size, showxticks=False,
            ticklabel_fontsize=ticklabel_size, legend_fontsize=legend_fontsize
        )
        curr_ax.set_ylim([0.15, 0.85])

        # 3rd row: intron position
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
        curr_ax = axs[2, i]
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

        # 4th row: % of bases with overlapping features
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
        curr_ax = axs[3, i]
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


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Plots intron feature boxplots / tools & truth categories'
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

    features_vs_detection(hx1_call, ipsc_call, out_dir, now)
    all_results_features_vs_truth(
        hx1_call, hx1_expr, ipsc_call, ipsc_expr, out_dir, now
    )
    filtered_groupedbysamp_featvstruth(hx1_call, ipsc_call, out_dir, now)

    print('\nprinting additional info')
    df_dict = {'iPSC': ipsc_call, 'HX1': hx1_call}
    all_spl = 0
    all_unspl = 0
    total_ints = 0
    all_called = 0
    single_tool = 1
    for samp, file in df_dict.items():
        print('\n\n{}'.format(samp))
        full_df = pd.read_table(file, sep='\t')
        all_spl += full_df.loc[full_df[_PERS] == 0]['intron'].nunique()
        all_unspl += full_df.loc[full_df[_PERS] == 1]['intron'].nunique()
        total_ints += full_df['intron'].nunique()
        print('all introns:', full_df['intron'].nunique())
        print(
            'fully spliced:',
            full_df.loc[full_df[_PERS] == 0]['intron'].nunique()
        )
        print(
            full_df.loc[full_df[_PERS] == 0]['intron'].nunique()
            / full_df['intron'].nunique()
        )
        print(
            'fully unspliced:',
            full_df.loc[full_df[_PERS] == 1]['intron'].nunique()
        )
        print(
            full_df.loc[full_df[_PERS] == 1]['intron'].nunique()
            / full_df['intron'].nunique()
        )
        overlap_df = full_df.groupby(['intron'])[[
            _SI_COL_BIN, _INT_COL_BIN, _IR_COL_BIN, _IRFS_COL_BIN,
             _KMA_COL_BIN, _RMA_COL_BIN, _SUP_COL_BIN, _MAJ_COL_BIN, _PERS
        ]].max().reset_index()
        overlap_df['ts'] = overlap_df.apply(
            lambda x:
                int(x[_SI_COL_BIN] > 0)
                + int(x[_INT_COL_BIN] > 0)
                + int(x[_IRFS_COL_BIN] > 0)
                + int(x[_KMA_COL_BIN] > 0)
                + int(x[_IR_COL_BIN] > 0)
                + int(x[_RMA_COL_BIN] > 0)
                + int(x[_SUP_COL_BIN] > 0)
                + int(x[_MAJ_COL_BIN] > 0),
            axis=1
        )
        tool_count_dict = overlap_df.groupby('ts')['intron'].count().to_dict()
        curr_allcalled = 0
        curr_singletool = 0
        for num_tools, num_introns in tool_count_dict.items():
            if num_tools > 0:
                all_called += num_introns
                curr_allcalled += num_introns
            if num_tools == 1:
                curr_singletool = num_introns
                single_tool += num_introns
        print(
            '{} of called introns called by exactly one tool ({}/{})'
            ''.format(
                curr_singletool / curr_allcalled,
                curr_singletool,
                curr_allcalled
            )
        )
        uncalled_persistent_introns = full_df.loc[
            (full_df[_PERS] > 0)
            & (full_df[_SI_COL_BIN] == 0)
            & (full_df[_INT_COL_BIN] == 0)
            & (full_df[_IR_COL_BIN] == 0)
            & (full_df[_IRFS_COL_BIN] == 0)
            & (full_df[_KMA_COL_BIN] == 0)
            & (full_df[_SUP_COL_BIN] == 0)
            & (full_df[_MAJ_COL_BIN] == 0)
            & (full_df[_RMA_COL_BIN] == 0)
        ]['intron'].nunique()
        pers_ints = full_df.loc[full_df[_PERS] > 0]['intron'].nunique()
        print(
            '{} persistent introns not called by any short-read tool ({} / {})'
            ''.format(
                uncalled_persistent_introns / pers_ints,
                uncalled_persistent_introns, pers_ints
            )
        )
        num_genes = full_df['gene_id'].nunique()
        print('\n{} genes with RI'.format(num_genes))

    print('\nall:')
    print('fully spliced:', all_spl)
    print(all_spl / total_ints)
    print('fully unspliced:', all_unspl)
    print(all_unspl / total_ints)
    print(
        '{} of called introns called by exactly one tool ({}/{})'
        ''.format(single_tool / all_called, single_tool, all_called)
    )
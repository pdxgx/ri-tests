#!/usr/bin/env python3

"""
SR_intron_sharedness.py
Python 3 code for plotting truth categories of shared called RIs

SAMPLE RUN:
time python ../../intronomer-paper/benchmarking_data/SR_intron_sharedness.py
-s target_genes_LR_annotated_granges-lrmap_sr-5-methods_SRR6026510-ipsc.csv
-o ../paper_results
"""
import argparse
from datetime import datetime
import matplotlib.pyplot as plt
import os
import pandas as pd
import upsetplot as up


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


def process_full_shortread_results(ranges_file, all_sr_introns):
    cols_to_use = ['intron']
    if all_sr_introns:
        coltype = _CONT_COL
    else:
        coltype = _BIN_COL
    cols_to_use += [tooldict[coltype] for tooldict in _TOOLS.values()]
    result_df = pd.read_table(ranges_file, sep=',', usecols=cols_to_use)
    result_df.fillna(0, inplace=True)
    intron_sets = {}
    for tooldict in _TOOLS.values():
        intron_sets[tooldict[_INTS]] = set(
            result_df.loc[result_df[tooldict[coltype]] > 0]['intron'].tolist()
        )
    return intron_sets


def tp_fp_stacked_barplots(lr_introns, intron_sets, flag, out_dir, now,
                           df_location=None):
    upset_df = up.from_contents({
        'pacbio': lr_introns,
        'IRFinder-S': intron_sets[_TOOLS[_IRFS][_INTS]],
        'IntEREst': intron_sets[_TOOLS[_INT][_INTS]],
        'superintronic': intron_sets[_TOOLS[_SI][_INTS]],
        'kma': intron_sets[_TOOLS[_KMA][_INTS]],
        'iREAD': intron_sets[_TOOLS[_IR][_INTS]],
    })
    if df_location:
        persistence = 'max_intron_persistence'
        IR_df = pd.read_table(
            df_location, usecols=['intron', persistence], sep=','
        )
        intron_ir = IR_df.set_index('intron')[persistence].to_dict()
        upset_df['persistence'] = upset_df['id'].apply(
            lambda x: intron_ir[x]
        )
    if 'HX1' in flag:
        label_color = '#a90308'
        flag = 'HX1'
        label_height = 11000
        label_xpos = 1
    else:
        label_color = '#0a437a'
        flag = 'iPSC'
        label_height = 14000
        label_xpos = 1

    upset_df.reset_index(inplace=True)

    tool_cols = ['IRFinder-S', 'IntEREst', 'superintronic', 'kma', 'iREAD']
    tp_df = upset_df.loc[upset_df['pacbio'] == True][tool_cols].copy()
    fp_df = upset_df.loc[upset_df['pacbio'] == False][tool_cols].copy()
    sr_tps = [len(tp_df[tp_df.sum(axis=1) == 0])]
    sr_fps = [len(fp_df[fp_df.sum(axis=1) == 0])]
    sharedness_levels = [1, 2, 3, 4, 5]
    for val in sharedness_levels:
        sr_tps.append(len(tp_df[tp_df.sum(axis=1) >= val]))
        sr_fps.append(len(fp_df[fp_df.sum(axis=1) >= val]))

    sr_fns = [sr_tps[0], 0, 0, 0, 0, 0]
    sr_tps[0] = 0
    labels = ['0\n(LR only)', '1+', '2+', '3+', '4+', '5']
    shift = 0.2
    label_locs = [1, 2, 3, 4, 5, 6]
    left_locs = [loc - shift for loc in label_locs]
    right_locs = [loc + shift for loc in label_locs]
    tp_col = '#9db92c'
    fp_col = '#d8788a'
    axislabel_size = 18
    ticklabel_size = 16
    legend_fontsize = 16
    fig, ax = plt.subplots()
    ax.yaxis.grid(True)
    ax.set_axisbelow(True)
    width = 0.35
    tp_lab = 'TPs'
    fp_lab = 'FPs'
    fn_lab = 'FNs'
    plt.bar(left_locs, sr_tps, width, color=tp_col, label=tp_lab)
    plt.bar(right_locs, sr_fps, width, color=fp_col, label=fp_lab)
    plt.bar(
        label_locs, sr_fns, width, color='#6d9cc6',
        label=fn_lab, tick_label=labels
    )
    ax.set_ylabel('Intron count', fontsize=axislabel_size)
    ax.set_xlabel('SR tool count', fontsize=axislabel_size)
    plt.yticks(fontsize=ticklabel_size)
    plt.xticks(fontsize=ticklabel_size)

    plt.text(
        label_xpos, label_height, flag,
        # color=label_color,
        color='black',
        fontsize=axislabel_size+2,
        # fontweight='bold',
        ha='center', va='center',
        bbox={'facecolor': 'white', 'alpha': 0, 'edgecolor': 'white', 'pad': 2}
    )
    ax.legend(prop={'size': legend_fontsize})
    ax.set_yscale('log')
    fig = plt.gcf()
    fig_file = os.path.join(
        out_dir, 'barplot_srcount_thresh0.1_{}_{}.pdf'.format(flag, now)
    )
    fig.savefig(fig_file, bbox_inches='tight', pad_inches=.1)
    return


def compare_long_to_short_reads(all_ranges, now, out_dir, all_sr):
    persistence = 'max_intron_persistence'
    tx_df = pd.read_table(
        all_ranges, sep=',',
        usecols=['intron', 'width', persistence, 'motif']
    )
    tx_df = tx_df.loc[tx_df[persistence] > 0]
    if 'hx' in all_ranges:
        flag = 'HX1'
    else:
        flag = 'iPSC'
    if all_sr:
        flag += '_all'
    else:
        flag += '_filtered'
    intron_sets = process_full_shortread_results(all_ranges, all_sr)
    threshold_levels = [0.1]
    for threshold in threshold_levels:
        lr_introns = set(
            tx_df.loc[tx_df[persistence] >= threshold]['intron'].tolist()
        )
        tp_fp_stacked_barplots(
            lr_introns, intron_sets, flag, out_dir, now, df_location=all_ranges
        )
    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Creates truth category barplots of shared called RIs'
    )
    parser.add_argument(
        '--range-summarized-shortread-results', '-s',
        help='short-read detection results summarized by intron region'
    )
    parser.add_argument(
        '--output-directory', '-o',
        help='directory in which to store output files and figures; if no '
             'directory specified, this will be set as the same directory '
             'containing the range summarized results.'
    )
    parser.add_argument(
        '--all-shortread', '-a', action='store_true',
        help='select this option to look at all short read potential RIs, not '
             'filtered RIs only.'
    )

    args = parser.parse_args()
    all_ranges = args.range_summarized_shortread_results
    all_sr = args.all_shortread
    output_dir = args.output_directory

    now = datetime.now().strftime('%m-%d-%Y_%H.%M.%S')
    if not output_dir:
        output_dir = os.path.dirname(all_ranges)
    compare_long_to_short_reads(all_ranges, now, output_dir, all_sr)

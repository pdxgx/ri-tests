#!/usr/bin/env python3

"""
FigS10_called_vs_all_RIs.py
Python 3 code for barplots of potential vs called/filtered RIs

SAMPLE RUN:
time python ri-tests/results/supplemental_figures/plot_filtered_RIs.py
-b HX1_final/target_genes_LR_annotated_granges-lrmap_sr-5-methods
_SRR2911306-hx1.csv
-B iPSC_final/target_genes_LR_annotated_granges-lrmap_sr-5-methods
_SRR6026510-ipsc.csv 
-o .

"""
import argparse
from datetime import datetime
import matplotlib.pyplot as plt
import os
import pandas as pd


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


def plot_filtered_RI_counts(filepath, out_dir):
    if 'HX1' in filepath:
        flag = 'HX1'
        label_height = 11000
    else:
        flag = 'iPSC'
        label_height = 7050
    df = pd.read_table(filepath, sep=',')
    potential_RIs = []
    called_RIs = []
    labels = []
    colors = []
    color_dict = {
        _IRFS: '#e1665d', _SI: '#f8b712', _IR: '#689404', _KMA: '#33a2b7',
        _INT: '#745bad', _SUP: '#ABDDFF', _MAJ: '#FFC0C0', _RMA: '#DAF7A6'
    }
    for tool in _TOOL_COLUMNS:
        potential_RIs.append(
            df.loc[df[_TOOLS[tool][_CONT_COL]] > 0]['intron'].nunique()
        )
        called_RIs.append(
            df.loc[df[_TOOLS[tool][_BIN_COL]] > 0]['intron'].nunique()
        )
        labels.append(tool)
        colors.append(color_dict[tool])

    # hat = '...'
    hat = 'oo'
    shift = 0.2
    label_locs = [1, 2, 3, 4, 5, 6, 7, 8]
    left_locs = [loc - shift for loc in label_locs]
    right_locs = [loc + shift for loc in label_locs]
    axislabel_size = 18
    ticklabel_size = 16
    legend_fontsize = 14
    fig, ax = plt.subplots()
    width = 0.35
    plt.bar(left_locs, potential_RIs, width, color=colors, hatch=hat)
    plt.bar(right_locs, called_RIs, width, color=colors)
    plt.bar(
        [1], [0], [0], color='white', edgecolor='black',
        hatch=hat, label='potential RIs'
        #label='nonzero expression RIs'
    )
    plt.bar(
        [2], [0], [0], color='white', edgecolor='black',
        hatch='', label='called RIs'
    )
    ax.set_ylabel('RI count', fontsize=axislabel_size)
    plt.yticks(fontsize=ticklabel_size)
    plt.xticks(
        label_locs, labels, rotation=45, ha="right" , fontsize=ticklabel_size
    )
    plt.text(
        1, label_height, flag,
        color='black',
        fontsize=axislabel_size + 2,
        ha='center', va='center',
        bbox={'facecolor': 'white', 'alpha': 0, 'edgecolor': 'white', 'pad': 2}
    )
    ax.legend(prop={'size': legend_fontsize})
    # ax.set_yscale('log')
    fig = plt.gcf()
    fig_file = os.path.join(
        out_dir, 'barplot_calledvspotentialRIs_{}_{}.pdf'.format(flag, now)
    )
    fig.savefig(fig_file, bbox_inches='tight', pad_inches=.1)
    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Plots counts of filtered vs all RIs for each SR tool'
    )
    parser.add_argument(
        '--HX1-binary-results', '-b',
        help='range-summarized, long read-annotated short read results for HX1'
    )
    parser.add_argument(
        '--iPSC-binary-results', '-B',
        help='range-summarized, long read-annotated short read results for '
             'iPSC'
    )
    parser.add_argument(
        '--output-directory', '-o',
        help='directory in which to store output files and figures; if no '
             'directory specified, this will be set as the same directory '
             'containing the range summarized results.'
    )

    args = parser.parse_args()
    hx1_bin = args.HX1_binary_results
    ipsc_bin = args.iPSC_binary_results
    out_dir = args.output_directory

    now = datetime.now().strftime('%m-%d-%Y_%H.%M.%S')

    plot_filtered_RI_counts(hx1_bin, out_dir)
    plot_filtered_RI_counts(ipsc_bin, out_dir)

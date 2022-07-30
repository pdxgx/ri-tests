#!/usr/bin/env python3

"""
long_vs_short_read_RIs.py
Python 3 code for comparing LR data w/ SR data, and calculating performance

SAMPLE RUN:
time python 
ri-tests/data_preparation/processing/long_vs_short_reads/
    long_vs_short_read_RIs.py
-s target_genes_LR_annotated_granges-lrmap_sr-5-methods_SRR2911306-hx1.csv
-o HX1_final/called_RIs -a

Add -P for persistence-vs-intron length plotting
Add -a for using all/unfiltered nonzero expression RIs from short read tools
"""
import argparse
from datetime import datetime
import matplotlib.pyplot as plt
from matplotlib import style
import os
import pandas as pd
from scipy import stats
import upsetplot as up


_STRAND_MAP = {True: '-', False: '+'}


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


def process_full_shortread_results(ranges_file, all_ints):
    if all_ints:
        target_col = _CONT_COL
    else:
        target_col = _BIN_COL
    extra_cols = [_TOOLS[col][target_col] for col in _TOOL_COLUMNS]
    cols_to_use = ['intron'] + extra_cols
    result_df = pd.read_table(ranges_file, sep=',', usecols=cols_to_use)
    result_df.fillna(0, inplace=True)
    intron_sets = {}
    for colname in _TOOL_COLUMNS:
        intron_sets[_TOOLS[colname][_INTS]] = set(
            result_df.loc[
                result_df[_TOOLS[colname][target_col]] > 0
            ].copy()['intron'].tolist()
        )
    return intron_sets


def make_upset_plots(lr_introns, intron_sets, flag, thresh, out_dir, now,
                     set_flag, df_location=None):
    upset_df = up.from_contents({
        'pacbio': lr_introns,
        'superintronic': intron_sets[_TOOLS[_SI][_INTS]],
        'IRFinder-S': intron_sets[_TOOLS[_IRFS][_INTS]],
        'IntEREst': intron_sets[_TOOLS[_INT][_INTS]],
        'iREAD': intron_sets[_TOOLS[_IR][_INTS]],
        'kma': intron_sets[_TOOLS[_KMA][_INTS]],
        _RMA: intron_sets[_TOOLS[_RMA][_INTS]],
        _MAJ: intron_sets[_TOOLS[_MAJ][_INTS]],
        _SUP: intron_sets[_TOOLS[_SUP][_INTS]],
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
    # save_data = False
    save_data = True
    if save_data:
        outfile = os.path.join(
            out_dir, 'upset_data_{}_thresh{}_{}_{}.csv'.format(
                flag, thresh, set_flag, now
            )
        )
        with open(outfile, 'w') as output:
            upset_df.to_csv(output, sep=',')
        # return

    if flag == 'HX1':
        edge = '#fcb001'
        fill = '#a90308'
    else:
        edge = '#c1f80a'
        fill = '#0a437a'

    # Full plots with all overlaps
    upset_plot = up.UpSet(
        upset_df, subset_size='count',
        # sort_by=None,
        sort_categories_by=None,
        sort_by='cardinality',
        facecolor=fill, show_counts=True
    )
    upset_plot.style_subsets(
        present="pacbio", edgecolor=edge, linewidth=1.5,
        label="{} true positives".format(flag)
    )
    upset_plot.plot()
    fig = plt.gcf()
    fig_file = os.path.join(
        out_dir, 'upset_overlaps_all_outlines_{}_thresh{}_{}_{}.pdf'.format(
            flag, thresh, set_flag, now
        )
    )
    fig.savefig(fig_file, bbox_inches='tight', pad_inches=.1)

    # Partial plots with pacbio-overlapping data only
    upset_df = upset_df.xs(True, level=0, axis=0, drop_level=False)
    upset_plot = up.UpSet(
        upset_df, subset_size='count',
        # sort_by='cardinality',
        # sort_by=None,
        sort_by='degree',
        sort_categories_by=None,
        facecolor=fill, show_counts=True
    )
    upset_plot.plot()
    fig = plt.gcf()
    fig_file = os.path.join(
        out_dir, 'upset_overlaps_TP_{}_thresh{}_{}_{}.pdf'.format(
            flag, thresh, set_flag, now
        )
    )
    fig.savefig(fig_file, bbox_inches='tight', pad_inches=.1)
    return


def sr_vs_lr_intronsets(tx_df, all_ranges, now, out_dir, all_ints):
    if 'hx' in all_ranges:
        flag = 'HX1'
    else:
        flag = 'iPSC'
    if all_ints:
        set_flag = 'nonzero'
        tool_cols = [_TOOLS[col][_CONT_COL] for col in _TOOL_COLUMNS]
    else:
        set_flag = 'called'
        tool_cols = [_TOOLS[col][_BIN_COL] for col in _TOOL_COLUMNS]

    intron_sets = process_full_shortread_results(all_ranges, all_ints)
    threshold_levels = [0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0]
    threshold_levels = threshold_levels[::-1]
    precision = 'max_intron_persistence'
    tx_df = tx_df.loc[tx_df[precision] > 0]
    pr_dict = {_THRSH: threshold_levels, _LR_COUNT: []}
    for col in _TOOL_COLUMNS:
        pr_dict[_TOOLS[col][_OVERLAP]] = []
        pr_dict[_TOOLS[col][_PREC]] = []
        pr_dict[_TOOLS[col][_REC]] = []

    cols_to_use = tool_cols + [
        'intron', precision, 'width', 'numreads.median',
        'transcript', 'gene_id', 'intron_position_in_tx',
        'motif', 'canonical_motif',
        'long_reads_per_gene', 'long_reads_per_transcript',
        'short_read_gene_median_coverage'
    ]
    summary_df = pd.read_table(all_ranges, sep=',', usecols=cols_to_use)
    for threshold in threshold_levels:
        lr_introns = set(
            tx_df.loc[tx_df[precision] >= threshold]['intron'].tolist()
        )
        relevant_ints = len(lr_introns)
        pr_dict[_LR_COUNT].append(relevant_ints)
        for col in _TOOL_COLUMNS:
            recalled_introns = intron_sets[_TOOLS[col][_INTS]]
            recalled_int_count = len(recalled_introns)
            true_pos = recalled_introns.intersection(lr_introns)
            overlap = len(true_pos)
            pr_dict[_TOOLS[col][_OVERLAP]].append(overlap)
            pr_dict[_TOOLS[col][_PREC]].append(overlap / recalled_int_count)
            pr_dict[_TOOLS[col][_REC]].append(overlap / relevant_ints)
            if threshold == 0.1:
                false_pos = recalled_introns.difference(lr_introns)
                false_neg = lr_introns.difference(recalled_introns)
                summary_df[_TOOLS[col][_TP]] = summary_df['intron'].apply(
                    lambda x: int(x in true_pos)
                )
                summary_df[_TOOLS[col][_FP]] = summary_df['intron'].apply(
                    lambda x: int(x in false_pos)
                )
                summary_df[_TOOLS[col][_FN]] = summary_df['intron'].apply(
                    lambda x: int(x in false_neg)
                )
        print('making upset with {} long read introns'.format(relevant_ints))
        make_upset_plots(
            lr_introns, intron_sets, flag, threshold, out_dir, now, set_flag,
            df_location=all_ranges
        )
    summary_file = os.path.join(
        out_dir, '{}_RI_data_summary_{}.tsv'.format(set_flag, flag)
    )
    summary_df.to_csv(summary_file, index=False, sep='\t')

    threshold_df = pd.DataFrame(pr_dict)
    threshold_df[_IR_COUNT] = len(intron_sets[_TOOLS[_IR][_INTS]])
    threshold_df[_INT_COUNT] = len(intron_sets[_TOOLS[_INT][_INTS]])
    threshold_df[_KMA_COUNT] = len(intron_sets[_TOOLS[_KMA][_INTS]])
    threshold_df[_SI_COUNT] = len(intron_sets[_TOOLS[_SI][_INTS]])
    threshold_df[_IRFS_COUNT] = len(intron_sets[_TOOLS[_IRFS][_INTS]])
    threshold_df[_MAJ_COUNT] = len(intron_sets[_TOOLS[_MAJ][_INTS]])
    threshold_df[_SUP_COUNT] = len(intron_sets[_TOOLS[_SUP][_INTS]])
    threshold_df[_RMA_COUNT] = len(intron_sets[_TOOLS[_RMA][_INTS]])

    threshold_df = threshold_df[[
        _THRSH, _LR_COUNT, _IR_COUNT, _IR_OVERLAP, _IR_PREC, _IR_REC,
        _INT_COUNT, _INT_OVERLAP, _INT_PREC, _INT_REC, _SI_COUNT, _SI_OVERLAP,
        _SI_PREC, _SI_REC, _KMA_COUNT, _KMA_OVERLAP, _KMA_PREC, _KMA_REC,
        _IRFS_COUNT, _IRFS_OVERLAP, _IRFS_PREC, _IRFS_REC,
        _RMA_COUNT, _RMA_OVERLAP, _RMA_PREC, _RMA_REC,
        _MAJ_COUNT, _MAJ_OVERLAP, _MAJ_PREC, _MAJ_REC,
        _SUP_COUNT, _SUP_OVERLAP, _SUP_PREC, _SUP_REC
    ]]
    df_file = os.path.join(
        out_dir, 'precision_recall_df_{}_{}_{}.csv'.format(flag, set_flag, now)
    )
    threshold_df.to_csv(df_file, index=False, sep=',')
    return


def sr_vs_lr_df_filtering(all_ranges, now, out_dir, all_ints):
    if 'hx' in all_ranges:
        flag = 'HX1'
    else:
        flag = 'iPSC'
    persist = 'max_intron_persistence'
    if all_ints:
        target_col = _CONT_COL
        set_flag = 'nonzero'
    else:
        target_col = _BIN_COL
        set_flag = 'called'
    extra_cols = [_TOOLS[col][target_col] for col in _TOOL_COLUMNS]
    cols_to_use = ['intron', persist] + extra_cols
    res_df = pd.read_table(all_ranges, sep=',', usecols=cols_to_use)
    res_df.fillna(0, inplace=True)
    threshold_levels = [0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0]
    threshold_levels = threshold_levels[::-1]
    pr_dict = {_THRSH: threshold_levels, _LR_COUNT: []}
    for col in _TOOL_COLUMNS:
        pr_dict[_TOOLS[col][_TP]] = []
        pr_dict[_TOOLS[col][_FP]] = []
        pr_dict[_TOOLS[col][_FN]] = []
        pr_dict[_TOOLS[col][_PREC]] = []
        pr_dict[_TOOLS[col][_REC]] = []
        pr_dict[_TOOLS[col][_FSCORE]] = []
    for threshold in threshold_levels:
        pr_dict[_LR_COUNT].append(
            res_df.loc[(res_df[persist] >= threshold)
                        & (res_df[persist] > 0)]['intron'].nunique()
        )
        for col in _TOOL_COLUMNS:
            true_pos = res_df.loc[
                (res_df[_TOOLS[col][target_col]] > 0)
                & (res_df[persist] >= threshold)
                & (res_df[persist] > 0)
            ]['intron'].nunique()
            false_pos = res_df.loc[
                (res_df[_TOOLS[col][target_col]] > 0) &
                ((res_df[persist] < threshold) | (res_df[persist] == 0))
            ]['intron'].nunique()
            false_neg = res_df.loc[
                (res_df[_TOOLS[col][target_col]] == 0)
                & (res_df[persist] > 0)
                & (res_df[persist] >= threshold)
            ]['intron'].nunique()
            pr_dict[_TOOLS[col][_TP]].append(true_pos)
            pr_dict[_TOOLS[col][_FP]].append(false_pos)
            pr_dict[_TOOLS[col][_FN]].append(false_neg)
            precision = true_pos / (true_pos + false_pos)
            recall = true_pos / (true_pos + false_neg)
            pr_dict[_TOOLS[col][_PREC]].append(precision)
            pr_dict[_TOOLS[col][_REC]].append(recall)
            try:
                pr_dict[_TOOLS[col][_FSCORE]].append(
                    (2 * precision * recall) / (precision + recall)
                )
            except ZeroDivisionError:
                pr_dict[_TOOLS[col][_FSCORE]].append(0)
    out_df = pd.DataFrame(pr_dict)

    col_order = []
    coltypeorder = [_COUNT, _TP, _FP, _FN, _PREC, _REC, _FSCORE]
    for col in _TOOL_COLUMNS:
        out_df[_TOOLS[col][_COUNT]] = res_df.loc[
            res_df[_TOOLS[col][target_col]] > 0
        ]['intron'].nunique()
        col_order.extend([_TOOLS[col][coltype] for coltype in coltypeorder])

    # out_df[_IR_COUNT] = res_df.loc[
    #     res_df[_TOOLS[_IR][target_col]] > 0
    # ]['intron'].nunique()
    # out_df[_INT_COUNT] = res_df.loc[
    #     res_df[_TOOLS[_INT][target_col]] > 0
    # ]['intron'].nunique()
    # out_df[_KMA_COUNT] = res_df.loc[
    #     res_df[_TOOLS[_KMA][target_col]] > 0
    # ]['intron'].nunique()
    # out_df[_SI_COUNT] = res_df.loc[
    #     res_df[_TOOLS[_SI][target_col]] > 0
    # ]['intron'].nunique()
    # out_df[_IRFS_COUNT] = res_df.loc[
    #     res_df[_TOOLS[_IRFS][target_col]] > 0
    # ]['intron'].nunique()

    # out_df = out_df[[
    #     _THRSH, _LR_COUNT,
    #     _IR_COUNT, _IR_TP, _IR_FP, _IR_FN, _IR_PREC, _IR_REC, _IR_F,
    #     _INT_COUNT, _INT_TP, _INT_FP, _INT_FN, _INT_PREC, _INT_REC, _INT_F,
    #     _SI_COUNT, _SI_TP, _SI_FP, _SI_FN, _SI_PREC, _SI_REC, _SI_F,
    #     _KMA_COUNT, _KMA_TP, _KMA_FP, _KMA_FN, _KMA_PREC, _KMA_REC, _KMA_F,
    #     _IRFS_COUNT, _IRFS_TP, _IRFS_FP, _IRFS_FN, _IRFS_PREC, _IRFS_REC,
    #     _IRFS_F
    # ]]
    out_df = out_df[[_THRSH, _LR_COUNT] + col_order]
    out_df = out_df.transpose()
    df_file = os.path.join(
        out_dir, 'precision_recall_dffilterversion_{}_{}_{}.csv'.format(
            flag, set_flag, now
        )
    )
    with open(df_file, 'w') as output:
        out_df.to_csv(output, header=False, sep=',')
    return


def scatter_with_regression_and_size(scatter_df, xval, yval, figfile, chrom='',
                                     x_axis_label="5' --> 3' intron position",
                                     y_axis_label='retention persistence'):
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
    x = scatter_df[xval]
    y = scatter_df[yval]
    print(min(x), max(x), min(y), max(y))
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
    plt.clf()
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
    print(figfile)
    fig.savefig(figfile)
    plt.clf()
    return


def calculate_width(intron_coords):
    left, right =intron_coords.split(':')[1].split('-')
    return 1 + int(right) - int(left)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='compares long and short read set overlaps.'
    )
    parser.add_argument(
        '--range-summarized-shortread-results', '-s',
        help='short-read detection results summarized by intron region'
    )
    parser.add_argument(
        '--plot-persistence-vs-all', '-P', action='store_true',
        help='create persistence-vs-various scatterplots'
    )
    parser.add_argument(
        '--output-directory', '-o',
        help='directory in which to store output files and figures; if no '
             'directory specified, this will be set as the same directory '
             'containing the range summarized results.'
    )
    parser.add_argument(
        '--all-nonzero-introns', '-a', action='store_true',
        help='Select this option to analyze all short read nonzero-expression '
             'introns vs. (default) filtered, called RIs'
    )

    args = parser.parse_args()
    all_ranges = args.range_summarized_shortread_results
    extra_plots = args.plot_persistence_vs_all
    output_dir = args.output_directory
    all_ints = args.all_nonzero_introns

    now = datetime.now().strftime('%m-%d-%Y_%H.%M.%S')
    if not output_dir:
        output_dir = os.path.dirname(all_ranges)
    persistence = 'max_intron_persistence'
    tx_df = pd.read_table(
        all_ranges, sep=',', usecols=['intron', persistence, 'motif']
    )
    tx_df[_WIDTH] = tx_df['intron'].apply(lambda x: calculate_width(x))
    data_df = pd.read_table(all_ranges, sep=',', usecols=['gene_id'])
    sr_vs_lr_intronsets(tx_df, all_ranges, now, output_dir, all_ints)
    sr_vs_lr_df_filtering(all_ranges, now, output_dir, all_ints)

    if extra_plots:
        fig_file = os.path.join(
            output_dir, 'persistence_vs_intronsize_{}.pdf'.format(now)
        )
        scatter_with_regression_and_size(
            tx_df, 'width', 'max_intron_persistence', fig_file,
            x_axis_label="intron width"
        )

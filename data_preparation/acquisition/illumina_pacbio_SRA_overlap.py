#!/usr/bin/env python3

"""
illumina_pacbio_sra_overlap.py

sample run:
time python intronomer-paper/benchmarking_data/illumina_pacbio_sra_overlap.py
-i intronomer_testfiles/SraRunInfo_illumina_7-20-20.csv
-p intronomer_testfiles/SraRunInfo_pacbio_7-20-20.csv

"""

import argparse
from datetime import datetime
import os
import pandas as pd

_BASE_COLS = ['SRAStudy', 'Sample', 'BioSample', 'SampleName']
_PB_EXPT = 'pacbio_expts'
_IL_EXPT = 'illumina_expts'
_EXPT = 'Experiment'


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Find SRA studies with both illumina and pacbio data.'
    )
    parser.add_argument(
        '--output-path', '-o', default='./',
        help='Give path to store output study file.'
    )
    parser.add_argument(
        '--illumina-studies', '-i', required=True,
        help='"RunInfo" .csv file downloaded from an SRA query for '
             'Organism: human, '
             'Source: transcriptomic, '
             'Strategy: rna seq, '
             'Platform: illumina '
             'Access: public '
             'Selection: polya '
             '(query done 7/13/21).'
    )
    parser.add_argument(
        '--pacbio-studies', '-p', required=True,
        help='"RunInfo" .csv file downloaded from an SRA query for '
             'Organism: human, '
             'Source: transcriptomic, '
             'Strategy: rna seq, '
             'Platform: pacbio smrt '
             'Access: public '
             '(query done 7/13/21).'
    )

    args = parser.parse_args()
    outpath = args.output_path
    illum_file = args.illumina_studies
    pacb_file = args.pacbio_studies

    now = datetime.now().strftime('%m-%d-%Y_%H.%M.%S')

    columns_to_read = _BASE_COLS + [_EXPT]
    pb_df = pd.read_table(pacb_file, sep=',', usecols=columns_to_read)
    pb_df.rename({_EXPT: _PB_EXPT}, axis=1, inplace=True)
    il_df = pd.read_table(illum_file, sep=',', usecols=columns_to_read)
    il_df.rename({_EXPT: _IL_EXPT}, axis=1, inplace=True)
    df = pb_df.merge(il_df, on=_BASE_COLS, how='inner')
    # with open(os.path.join(outpath, 'all_illumina_pacbio.csv'), 'w') as output:
    #     df.to_csv(output, index=False)

    grp_cols = [_IL_EXPT] + _BASE_COLS
    df = df.groupby(grp_cols)[_PB_EXPT].apply(';'.join).reset_index()
    grp_cols = [_PB_EXPT] + _BASE_COLS
    df = df.groupby(grp_cols)[_IL_EXPT].apply(';'.join).reset_index()
    df = df.loc[df['SRAStudy'] != 'SRP012412']
    df['pacbio_ex_count'] = df[_PB_EXPT].apply(
        lambda x: len(x.split(';'))
    )
    df['illumina_ex_count'] = df[_IL_EXPT].apply(
        lambda x: len(x.split(';'))
    )
    out_file = os.path.join(outpath, 'illumina_pacbio_merged.csv')
    with open(out_file, 'w') as output:
        df.to_csv(output, index=False)

    il_expts = df[_IL_EXPT].unique().tolist()
    pb_expts = df[_PB_EXPT].unique().tolist()
    il_df = il_df.loc[il_df[_IL_EXPT].isin(il_expts)]
    il_df.rename({_IL_EXPT: _EXPT}, axis=1, inplace=True)
    il_df['platform'] = 'illumina'
    pb_df = pb_df.loc[pb_df[_PB_EXPT].isin(pb_expts)]
    pb_df.rename({_PB_EXPT: _EXPT}, axis=1, inplace=True)
    pb_df['platform'] = 'pacbio_smrt'
    cc_df = pd.concat([il_df, pb_df])
    with open(os.path.join(outpath, 'concat_expts.csv'), 'w') as output:
        cc_df.to_csv(output, index=False)

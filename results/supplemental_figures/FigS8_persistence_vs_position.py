#!/usr/bin/env python3

"""
persistence_vs_position.py
Python 3 code for plotting calculated intron persistence vs. intron position

SAMPLE RUN:
time python intronomer-paper/benchmarking_data/persistence_vs_position.py
-p iPSC_final/processed_tx_df_iPSC_02-21-2022_21.03.50.csv -o paper_results
"""
import argparse
from datetime import datetime
import matplotlib.pyplot as plt
from matplotlib import style
import os
import pandas as pd
from scipy import stats


def scatter_with_regression_and_size(scatter_df, xval, yval, figfile, flag,
                                     chrom='',
                                     x_axis_label="5' --> 3' intron position",
                                     y_axis_label='intron persistence'):
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
    if flag == 'HX1':
        color = '#d3494e'
    else:
        color = '#448ee4'
    groups = scatter_df.groupby([xval, yval])
    scatter_dict = {'x': [], 'y': [], 'count': [], 'color': []}
    for group_index, group in groups:
        scatter_dict['x'].append(group_index[0])
        scatter_dict['y'].append(group_index[1])
        scatter_dict['count'].append(6 * (len(group) ** (1.0 / 1.5)))
        scatter_dict['color'].append(color)
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
        legend_label = '{}, all chromosomes'.format(flag)
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
    y_max = scatter_df.y.max() * 1.2
    x_max = scatter_df.x.max() * 1.1
    ax.set_ylim(ymin=-(0.05 * y_max), ymax=y_max)
    ax.set_xlim(xmin=-(0.05 * x_max), xmax=x_max)

    plt.ylabel(y_axis_label)
    fig = plt.gcf()
    print(figfile)
    fig.savefig(figfile)
    plt.clf()
    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Plots calculated persistence metric vs intron position.'
    )
    parser.add_argument(
        '--processed-tx-file', '-p',
        help='file with name "processed_tx_df__[timestamp]__.csv" '
             'generated previously by assign_RI_persistence_metric.py.'
    )
    parser.add_argument(
        '--output-directory', '-o',
        help='directory in which to store output files and figures; if no '
             'directory specified, this will be set as the same directory '
             'containing the range summarized results.'
    )

    args = parser.parse_args()
    read_map = args.processed_tx_file
    output_dir = args.output_directory

    now = datetime.now().strftime('%m-%d-%Y_%H.%M.%S')
    if not output_dir:
        output_dir = os.path.dirname(read_map)

    if 'HX1' in os.path.abspath(read_map):
        flag = 'HX1'
    else:
        flag = 'iPSC'
    read_tx_df = pd.read_table(read_map, sep=',')
    fig_name = (
        'intron_persistence_vs_position_allchroms_{}_{}.pdf'.format(flag, now)
    )
    fig_file = os.path.join(output_dir, fig_name)
    scatter_with_regression_and_size(
        read_tx_df, 'position', 'persistence', fig_file, flag
    )

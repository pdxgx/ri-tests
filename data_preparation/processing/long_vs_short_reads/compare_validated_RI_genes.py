#!/usr/bin/env python3

"""
compare_validated_RI_genes.py
Python 3 code for comparing called RIs against experimentally validated RIs

time python intronomer-paper/benchmarking_data/compare_validated_RI_genes.py
-H HX1_final/called_RIs/called_RI_data_summary_HX1featureannotated.tsv
-e HX1_final/target_genes_HX1.txt
-I iPSC_final/called_RIs/called_RI_data_summary_iPSCfeatureannotated.tsv
-l iPSC_final/target_genes_iPSC.txt
-g files/gencode.v35.annotation.gtf

"""
import argparse
from collections import namedtuple
from datetime import datetime
import os
import re

import pandas as pd


_GENE_SOURCE = {
    'AP1G2': 'Jeong2021',
    'CELF1': 'Li2021',
    'CTSD': 'Wong2014',
    'LBR': 'Wong2014',
    'SRSF7': 'Lejeune2001',
    # 2 introns in BAIAP
    'BAIAP2': 'Li2021',
    'CLASRP': 'Li2021',
    'CMTM5': 'Li2021',
    'EXD3': 'Li2021',
    'FAAH': 'Li2021',
    'FAHD2A': 'Li2021',
    'FAHD2B': 'Li2021',
    # GRAMD3 --> GRAMD2B
    'GRAMD2B': 'Li2021',
    'IGSF8': 'Li2021',
    'IL3RA': 'Li2021',
    'MAP1LC3A': 'Li2021',
    'MRAS': 'Li2021',
    # 2 introns in NDUFS7
    'NDUFS7': 'Li2021',
    'NIFK': 'Li2021',
    'POLR2F': 'Li2021',
    'PTGDS': 'Li2021',
    'QTRT1': 'Li2021',
    'RANBP3L': 'Li2021',
    'SEMA3B': 'Li2021',
    'SH3GL3': 'Li2021',
    'SLC6A12': 'Li2021',
    'SNRNP70': 'Li2021',
    # SPG20 --> SPART
    'SPART': 'Li2021',
    'TCFL5': 'Li2021',
    # 2 introns in PDIA2
    'PDIA2': 'Li2021',
    'TUG1': 'Dumbovic2021',
    'TERT': 'Dumbovic2021',
    'ZRSR2': 'Inoue2021',
    'LZTR1': 'Inoue2021',
    'APP': 'Buckley2011',
    'ATF4': 'Wong2014',
    'CACNA1H': 'Zhong2006',
    'CAMK2B': 'Zhong2006',
    'CREB1': 'Buckley2011',
    'EIF1': 'Wong2014',
    'FMR1': 'Buckley2011',
    'FANCA': 'Wong2014',
    'GABRG3': 'Buckley2011',
    'GRIK1': 'Buckley2011',
    'GRIN1': 'Buckley2011',
    'HP': 'Wong2014',
    'ID3': 'Forrest2004',
    'IL1B': 'Denis2005',
    'INS': 'Mansilla2005',
    'KCNMA1': 'Bell2008,Bell2011',
    'LMNB1': 'Wong2014',
    'ROBO3': 'Chen2008(seeBuckley2014ref12)',
}


GTFLine = namedtuple(
    'GTFLine', ['chr', 'left', 'right', 'strand', 'type', 'val_dict']
)


def parse_gtf_line(line):
    """Checks .gtf line validity and returns parsed elements"""
    line = line.strip()
    if not line or line.startswith('#'):
        return None
    if '#' in line:
        line = line.split('#')[0].strip()
    try:
        item = line.split('\t')
        chrom, left, right, strand = item[0], item[3], item[4], item[6]
        feature, values = item[2], item[8]
    except ValueError:
        return None
    left, right = int(left), int(right)
    if feature not in {'transcript', 'exon'} or left >= right:
        return None
    tag_value_pattern = r"(\w+)\s\"(\S+?)\"[;]"
    val_dict = dict(re.findall(tag_value_pattern, values))
    gtf_line = GTFLine(chrom, left, right, strand, feature, val_dict)
    return gtf_line


def gene_name_to_ids(gtf_file):
    """Extracts sorted annotated exon coordinates from .gtf file

    gtf_file (str): path to gencode annotation gtf.

    Creates dictionary mapping transcript ids to their exons & coordinates.
    Sorts and merges exons for each transcript.
    Returns dictionary searchable for exon chromosome, strand, and left and
    right coordinates.
    """
    gene_name_to_id = {}
    with open(gtf_file) as gtf:
        # Parse valid exon lines from the annotation file into a dict by
        # transcript_id
        for gtf_line in gtf:
            line = parse_gtf_line(gtf_line)
            if line is None:
                continue
            gene_id = line.val_dict['gene_id']
            gene_name = line.val_dict['gene_name']
            if gene_name not in gene_name_to_id.keys():
                gene_name_to_id[gene_name] = gene_id
    return gene_name_to_id


def assess_validatedRI_genes(binary_summary, gene_list, validated_genes,
                             gene_name_to_id, now):
    out_dir = os.path.dirname(gene_list)
    if 'HX1' in binary_summary:
        flag = 'HX1'
    else:
        flag = 'iPSC'
    print('\nstarting {}'.format(flag))
    target_genes = set()
    with open(gene_list) as genes:
        for line_string in genes:
            target_genes.add(line_string.strip())

    print('{} target genes for {}'.format(len(target_genes), flag))
    ri_df = pd.read_table(binary_summary, sep='\t')
    print('{} with potential RIs'.format(ri_df['gene_id'].nunique()))
    val_set = set(validated_genes)
    val_in_target = val_set.intersection(target_genes)
    val_notintarget = len(val_set.difference(target_genes))
    print(
        '{} of {} validated genes in target gene list; {} not present'.format(
            len(val_in_target), len(val_set), val_notintarget
        )
    )
    ri_df = ri_df.loc[ri_df['gene_id'].isin(val_set)]
    print(
        '{} validated genes w/ RIs:'.format(ri_df['gene_id'].nunique())
    )
    print(ri_df['gene_id'].unique())
    id_to_name = {val: key for key, val in gene_name_to_id.items()}
    ri_df['gene_name'] = ri_df['gene_id'].apply(lambda x: id_to_name[x])
    ri_df['source'] = ri_df['gene_name'].apply(lambda x: _GENE_SOURCE[x])
    outfile = os.path.join(
        out_dir, '{}_validatedONLY_RI_genes_{}.tsv'.format(flag, now)
    )
    ri_df.to_csv(outfile, index=False, sep='\t')
    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Compares called RIs against experimentally validated RIs.'
    )
    parser.add_argument(
        '--HX1-called-summary', '-H', help='called-RI summmary file for HX1'
    )
    parser.add_argument(
        '--HX1-gene-list', '-e',
        help='target gene list for HX1 generated by '
             'select_and_plot_target_genes.py'
    )
    parser.add_argument(
        '--iPSC-called-summary', '-I', help='called-RI summmary file for iPSC'
    )
    parser.add_argument(
        '--iPSC-gene-list', '-l',
        help='target gene list for iPSC generated by '
             'select_and_plot_target_genes.py'
    )
    parser.add_argument(
        '--gencode-file', '-g',
        help='.gtf file with GENCODE annotation.'
    )

    args = parser.parse_args()
    hx1_callfile = args.HX1_called_summary
    hx1_genelist = args.HX1_gene_list
    ipsc_callfile = args.iPSC_called_summary
    ipsc_genelist = args.iPSC_gene_list
    gtf_path = args.gencode_file

    now = datetime.now().strftime('%m-%d-%Y_%H.%M.%S')
    gene_name_to_id = gene_name_to_ids(gtf_path)
    validated_gene_ids = [
        gene_name_to_id[name] for name in _GENE_SOURCE.keys()
    ]
    print('{} validated genes'.format(len(set(validated_gene_ids))))
    assess_validatedRI_genes(
        hx1_callfile, hx1_genelist, validated_gene_ids, gene_name_to_id, now
    )
    assess_validatedRI_genes(
        ipsc_callfile, ipsc_genelist, validated_gene_ids, gene_name_to_id, now
    )

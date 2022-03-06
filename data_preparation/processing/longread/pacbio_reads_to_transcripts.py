#!/usr/bin/env python3

"""
pacbio_reads_to_transcripts.py
Python 3 code for mapping aligned pacbio reads to annotated transcripts

SAMPLE RUN:
time python intronomer-paper/benchmarking_data/pacbio_reads_to_transcripts.py
-g files/gencode.v34.annotation.gtf
-a SRP065930_SAMN04251426.merged.aligned.sorted.bam
"""
import argparse
from collections import defaultdict
from datetime import datetime
import os
import pandas as pd
import pysam


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


def extract_introns(gtf_file):
    """Extracts splice site annotations from .gtf file

    This function is a modified version of one that is part of HISAT.
    Copyright 2014, Daehwan Kim <infphilo@gmail.com>

    HISAT is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    HISAT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with HISAT.  If not, see <http://www.gnu.org/licenses/>.
    """
    genes = defaultdict(list)
    trans = {}
    intron_to_txs = {}
    tx_ranges = {}
    tx_to_introns = {}
    tx_to_gene = {}
    with open(gtf_file) as gtf:
        # Parse valid exon lines from the annotation file into a dict by
        # transcript_id
        for line in gtf:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            if '#' in line:
                line = line.split('#')[0].strip()

            try:
                items = line.split('\t')
                chrom, _, feature, left, right, _, strand, _, values = items
            except ValueError:
                continue
            left, right = int(left) - 1, int(right) + 1

            if feature != 'exon' or left >= right:
                continue

            vals_dict = {}
            for attr in values.split(';')[:-1]:
                attr, _, val = attr.strip().partition(' ')
                vals_dict[attr] = val.strip('"')

            if 'gene_id' not in vals_dict or 'transcript_id' not in vals_dict:
                continue

            transcript_id = vals_dict['transcript_id']
            if transcript_id not in trans:
                trans[transcript_id] = [chrom, strand, [[left, right]]]
                genes[vals_dict['gene_id']].append(transcript_id)
            else:
                trans[transcript_id][2].append([left, right])
            tx_to_gene[transcript_id] = vals_dict['gene_id']

    # Sort exons and merge where separating introns are <=5 bps
    for tran, [chrom, strand, exons] in trans.items():
        exons.sort()
        tmp_exons = [exons[0]]
        for i in range(1, len(exons)):
            if exons[i][0] - tmp_exons[-1][1] <= 5:
                tmp_exons[-1][1] = exons[i][1]
            else:
                tmp_exons.append(exons[i])
        trans[tran] = [chrom, strand, tmp_exons]

    # Calculate and print the unique junctions
    for tx_name, (chrom, strand, exons) in trans.items():
        tx_to_introns[tx_name] = set()
        if chrom not in intron_to_txs:
            intron_to_txs[chrom] = {}
            intron_to_txs[chrom]['+'] = defaultdict(set)
            intron_to_txs[chrom]['-'] = defaultdict(set)
            intron_to_txs[chrom]['.'] = defaultdict(set)
        for i in range(1, len(exons)):
            tx_to_introns[tx_name].add((exons[i - 1][1], exons[i][0]))
            left = exons[i - 1][1]
            right = exons[i][0]
            intron_to_txs[chrom][strand][(left, right)].add(tx_name)
        tx_ranges[tx_name] = [exons[0][0], exons[-1][1]]
        tx_to_introns[tx_name] = sorted(tx_to_introns[tx_name])

    return tx_ranges, intron_to_txs, tx_to_introns, tx_to_gene


def check_match_type(read_introns, read_left, read_right, tx_introns):
    skipped_splicing = False
    full_length = True
    target_introns = 0
    missing_introns = 0
    for i, intron in enumerate(tx_introns, 1):
        if intron in read_introns:
            target_introns += 1
            continue
        if intron[0] > read_left and intron[1] < read_right:
            skipped_splicing = True
            target_introns += 1
            missing_introns += 1
        else:
            full_length = False
    if skipped_splicing:
        match_type = 'skipped_splicing-'
    else:
        match_type = 'all_introns-'
    if full_length:
        match_type += 'full_length'
    else:
        match_type += 'partial'
    return match_type, target_introns, missing_introns


def longreads_to_isoforms(bam_file, intron_info, output_dir, now):
    now = datetime.now().strftime('%m-%d-%Y_%H.%M.%S')

    tx_ranges, intron_to_txs, tx_to_introns, tx_to_gene = intron_info
    if bam_file.endswith('bam'):
        print('binary')
        samfile = pysam.AlignmentFile(bam_file, 'rb')
    else:
        print('nonbinary')
        samfile = pysam.AlignmentFile(bam_file, 'r')
    nocigartups = 0
    total_nonexact_reads = 0
    matching_intron_set_read_counts = defaultdict(int)
    read_maps_encountered = set()
    all_reads = 0
    read_name_dict = defaultdict(set)
    tx_to_all_reads = defaultdict(set)
    best_read_match = {}
    tx_to_exact_reads = defaultdict(set)
    tx_to_partial_reads = defaultdict(set)
    RI_txs = defaultdict(int)
    RI_tx_set = set()
    read_info_dict = {
        'read_id': [],
        'orig_read_id': [],
        'read_length': [],
        'read_left': [],
        'read_right': [],
        'chrom': [],
        'strand': [],
        'read_introns': [],
    }
    tx_info_dict = {
        'transcript': [],
        'tx_length': [],
        'tx_left': [],
        'tx_right': [],
        'tx_introns': []
    }
    tx_to_read_id = {
        'transcript': [],
        'read_id': [],
        'match_type': [],
        'region_tx_intron_count': [],
        'unspliced_intron_count': [],
    }
    for i, read in enumerate(samfile):
        collect_read = False
        name = read.query_name
        read_str = read.to_string()
        sequence = read.query_alignment_sequence
        read_name_dict[name].add(sequence)
        if read_str in read_maps_encountered:
            continue
        assigned_read_number = i
        read_maps_encountered.add(read_str)
        all_reads += 1
        possible_txs = set()
        cigar_tups = read.cigartuples
        if not cigar_tups:
            nocigartups += 1
            continue
        left_pos = read.reference_start
        right_pos = read.reference_end
        ref_len = read.reference_length
        strand = _STRAND_MAP[read.is_reverse]
        chrom = read.reference_name
        curr_pos = left_pos
        introns = []
        for (c_type, length) in cigar_tups:
            if c_type == 3:
                i_left = curr_pos + 1
                curr_pos += length
                i_right = curr_pos
                introns.append((i_left, i_right))
            elif c_type in (0, 2):
                curr_pos += length

        for intron in introns:
            try:
                possible_txs.update(intron_to_txs[chrom][strand][intron])
            except KeyError:
                continue
        intron_set = set(introns)
        r_intron_str = ''
        for (left, right) in sorted(list(intron_set)):
            r_intron_str += '{}-{};'.format(left, right)
        matching_intron_set_read_counts[r_intron_str] += 1
        shortlist_txs = set()
        exact_txs = set()
        for tx in possible_txs:
            target_intron_set = set(tx_to_introns[tx])
            if intron_set == target_intron_set:
                exact_txs.add(tx)
                collect_read = True

        exact_tx_count = len(exact_txs)
        if exact_txs:
            if exact_tx_count == 1:
                best_tx = list(exact_txs)[0]
            else:
                min_len_diff = 1000000000000
                best_tx = ''
                for tx in exact_txs:
                    tx_left = tx_ranges[tx][0]
                    tx_right = tx_ranges[tx][1]
                    len_diff = abs(
                        (tx_right - tx_left) - (right_pos - left_pos)
                    )
                    if len_diff < min_len_diff:
                        best_tx = tx
                        min_len_diff = len_diff
            tx_to_read_id['transcript'].append(best_tx)
            tx_to_read_id['read_id'].append(assigned_read_number)
            tx_to_read_id['match_type'].append('all_introns-full_length')
            tx_to_read_id['region_tx_intron_count'].append(
                len(set(tx_to_introns[best_tx]))
            )
            tx_to_read_id['unspliced_intron_count'].append(0)
            best_read_match[assigned_read_number] = best_tx
            tx_to_all_reads[best_tx].add(read_str)
            tx_to_exact_reads[best_tx].add(assigned_read_number)
        else:
            total_nonexact_reads += 1
            for tx in possible_txs:
                target_intron_set = set(tx_to_introns[tx])
                if intron_set.issubset(target_intron_set):
                    shortlist_txs.add(tx)
                    collect_read = True

            if shortlist_txs:
                if len(shortlist_txs) == 1:
                    best_tx = list(shortlist_txs)[0]
                else:
                    min_len_diff = 1000000000000
                    best_tx = ''
                    for tx in shortlist_txs:
                        tx_left = tx_ranges[tx][0]
                        tx_right = tx_ranges[tx][1]
                        len_diff = abs(
                            (tx_right - tx_left) - (right_pos - left_pos)
                        )
                        if len_diff < min_len_diff:
                            best_tx = tx
                            min_len_diff = len_diff
                best_read_match[assigned_read_number] = best_tx
                tx_to_all_reads[best_tx].add(assigned_read_number)
                tx_to_partial_reads[best_tx].add(assigned_read_number)
                tx_to_read_id['transcript'].append(best_tx)
                tx_to_read_id['read_id'].append(assigned_read_number)
                match_info = check_match_type(
                    sorted(list(intron_set)), left_pos, right_pos,
                    sorted(list(tx_to_introns[best_tx]))
                )
                match_type, tx_count, missing_introns = match_info
                tx_to_read_id['match_type'].append(match_type)
                tx_to_read_id['region_tx_intron_count'].append(tx_count)
                tx_to_read_id['unspliced_intron_count'].append(missing_introns)
                if match_type.startswith('skipped_splicing'):
                    RI_txs[best_tx] += 1
                    RI_tx_set.add(best_tx)
        if collect_read:
            read_info_dict['read_id'].append(assigned_read_number)
            read_info_dict['orig_read_id'].append(name)
            read_info_dict['read_length'].append(ref_len)
            read_info_dict['read_left'].append(left_pos)
            read_info_dict['read_right'].append(right_pos)
            read_info_dict['chrom'].append(chrom)
            read_info_dict['strand'].append(strand)
            read_info_dict['read_introns'].append(r_intron_str)

    setlen_dict = defaultdict(int)
    over5s = 0
    over5_txs = set()
    total_gene_set = set()
    for tx, read_set in tx_to_all_reads.items():
        total_gene_set.add(tx_to_gene[tx])
        setlen_dict[len(read_set)] += 1
        if len(read_set) >=5:
            over5s += 1
            over5_txs.add(tx)
    print('\ntotal read support per transcript:')
    print('{} transcripts with 5 or more reads'.format(over5s))

    for tx in set(tx_to_read_id['transcript']):
        target_intron_set = set(tx_to_introns[tx])
        t_intron_str = ''
        for (left, right) in sorted(list(target_intron_set)):
            t_intron_str += '{}-{};'.format(left, right)
        tx_info_dict['transcript'].append(tx)
        tx_info_dict['tx_length'].append(tx_ranges[tx][1] - tx_ranges[tx][0])
        tx_info_dict['tx_left'].append(tx_ranges[tx][0])
        tx_info_dict['tx_right'].append(tx_ranges[tx][1])
        tx_info_dict['tx_introns'].append(t_intron_str)

    read_info_df = pd.DataFrame(read_info_dict)
    outfile = os.path.join(
        output_dir, 'read_info_nosequence_{}.tsv'.format(now)
    )
    with open(outfile, 'w') as output:
        read_info_df.to_csv(output, sep='\t', index=False)

    tx_info_df = pd.DataFrame(tx_info_dict)
    outfile = os.path.join(
        output_dir, 'tx_info_{}.tsv'.format(now)
    )
    with open(outfile, 'w') as output:
        tx_info_df.to_csv(output, sep='\t', index=False)

    tx_to_read_id_df = pd.DataFrame(tx_to_read_id)
    outfile = os.path.join(
        output_dir, 'tx_to_read_id_{}.tsv'.format(now)
    )
    tx_to_read_id_df['gene_id'] = tx_to_read_id_df['transcript'].apply(
        lambda x: tx_to_gene.get(x, '')
    )
    with open(outfile, 'w') as output:
        tx_to_read_id_df.to_csv(output, sep='\t', index=False)

    tx_groups = tx_to_read_id_df.groupby('transcript').count()
    tx_groups['reads_per_transcript'] = tx_groups['read_id']
    tx_groups.reset_index(inplace=True)
    tx_to_gene = tx_to_read_id_df.set_index('transcript')['gene_id'].to_dict()
    tx_groups['gene_id'] = tx_groups['transcript'].apply(
        lambda x: tx_to_gene[x]
    )
    tx_groups = tx_groups[['gene_id', 'transcript', 'reads_per_transcript']]
    gene_groups = tx_to_read_id_df.groupby('gene_id').count()
    gene_groups['reads_per_gene'] = gene_groups['read_id']
    gene_groups = gene_groups[['reads_per_gene']]
    gene_groups.reset_index(inplace=True)
    mergedf = pd.merge(gene_groups, tx_groups, how='outer', on='gene_id')
    outfile = os.path.join(output_dir, 'reads_per_gene_and_transcript.tsv')
    with open(outfile, 'w') as output:
        mergedf.to_csv(output, sep='\t', index=False)

    tx_to_read_id_df = pd.merge(
        tx_to_read_id_df, read_info_df, on=['read_id'], how='left'
    )
    tx_to_read_id_df = pd.merge(
        tx_to_read_id_df, tx_info_df, on=['transcript'], how='left'
    )
    ratios_dict = {}
    RI_tx_set = set()
    for tx in tx_to_read_id_df['transcript'].unique():
        mini_df = tx_to_read_id_df.loc[tx_to_read_id_df['transcript'] == tx]
        skips = len(
            mini_df.loc[mini_df['match_type'].isin(
                ['skipped_splicing-full_length', 'skipped_splicing-partial'])]
        )
        if skips > 0 and len(mini_df) > 4:
            RI_tx_set.add(tx)
            ratios_dict[tx] = skips / len(mini_df)
    tx_to_read_id_df = tx_to_read_id_df.loc[
        tx_to_read_id_df['transcript'].isin(RI_tx_set)
    ]
    tx_to_read_id_df['RI-read_ratio'] = tx_to_read_id_df['transcript'].apply(
        lambda x: '{:.2f}'.format(ratios_dict[x])
    )
    tx_to_read_id_df.sort_values(
        by=['RI-read_ratio'], ascending=False, inplace=True
    )
    tx_to_read_id_df.sort_values(
        by=['transcript'], ascending=True, inplace=True
    )
    tx_to_read_id_df.sort_values(
        by=['match_type'], ascending=False, inplace=True
    )
    outfile = os.path.join(
        output_dir, 'RI_txs_to_read_ids_final_{}.tsv'.format(now)
    )
    with open(outfile, 'w') as output:
        tx_to_read_id_df.to_csv(output, sep='\t', index=False)
    print(
        '\n{} transcripts with at least 5 reads and some potential IR'
        ''.format(tx_to_read_id_df['transcript'].nunique())
    )
    setlen_dict = defaultdict(int)
    over5s = 0
    for tx, read_set in tx_to_exact_reads.items():
        setlen_dict[len(read_set)] += 1
        if len(read_set) >=5:
            over5s += 1
    print('\nexact-match read support per transcript:')
    print('{} transcripts with 5 or more reads'.format(over5s))

    setlen_dict = defaultdict(int)
    over5s = 0
    for tx, read_set in tx_to_partial_reads.items():
        setlen_dict[len(read_set)] += 1
        if len(read_set) >=5:
            over5s += 1
    print('\npartial-match read support per transcript:')
    print('{} transcripts with 5 or more reads'.format(over5s))

    print('\ntotal reads collected: {}'.format(len(read_info_dict['read_id'])))
    print(
        'covering {} genes and {} transcripts'
        ''.format(len(total_gene_set), len(tx_to_all_reads))
    )
    return tx_to_read_id_df


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Extracts retained-intron exons from stringtie-assembled '
                    'transcripts.'
    )
    parser.add_argument(
        '--aligned-long-reads', '-a',
        help='Sorted .bam file with full length non-concatamer poly-a '
             'selected long reads.'
    )
    parser.add_argument(
        '--annotation-gtf', '-g',
        help='.gtf file with GENCODE annotation.'
    )

    args = parser.parse_args()
    bam_path = args.aligned_long_reads
    gtf_path = args.annotation_gtf

    now = datetime.now().strftime('%m-%d-%Y_%H.%M.%S')
    output_dir = os.path.dirname(bam_path)

    intron_info = extract_introns(gtf_path)
    read_tx_df = longreads_to_isoforms(
        bam_path, intron_info, output_dir, now
    )

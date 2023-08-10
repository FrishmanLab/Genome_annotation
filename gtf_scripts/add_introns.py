#!/usr/bin/env python

import sys
import pandas as pd
import re

def gtf_to_df_with_genes(gtf_file):
    column_names = ['Seqid', 'Source', 'Type', 'Start', 'End', 'Score', 'Strand', 'Frame', 'Suppl']
    df = pd.read_csv(gtf_file, sep='\t', index_col=False, names=column_names, dtype={'Start': int, 'End': int})
    
    gene_ids = pd.Series([], dtype='object')
    if df['Source'].eq('Helixer').any():
        gene_ids = extract_id(df['Suppl'], r'Parent=([^;]+)')
    else:
        gene_ids = extract_id(df['Suppl'], r'transcript_id "([^"]+)"')
    
    df['Genes'] = gene_ids
    return df

def extract_id(string_series, pattern):
    return string_series.str.extract(pattern, expand=False)

def add_introns(gtf_file):
    df = gtf_to_df_with_genes(gtf_file)
    output_dfs = []
    
    for _, sub_df in df.groupby('Genes'):
        df_exons = sub_df[sub_df['Type'] == 'exon'].sort_values('Start')
        if len(df_exons) > 1:
            end_values = df_exons['End'].values
            start_values = df_exons['Start'].values
            intron_mask = start_values[1:] - end_values[:-1] > 1
            intron_regions = list(zip(end_values[:-1][intron_mask] + 1, start_values[1:][intron_mask] - 1))
        else:
            intron_regions = []
        
        rows = []
        for _, row in sub_df.iterrows():
            if row['Type'] == 'exon':
                current_exon = (row['Start'], row['End'])
                for intron_region in intron_regions:
                    if current_exon[1] + 1 == intron_region[0]:
                        new_row = row.copy()
                        new_row['Type'] = 'intron'
                        new_row['Start'] = intron_region[0]
                        new_row['End'] = intron_region[1]
                        rows.append(row)
                        rows.append(new_row)
            else:
                rows.append(row)
        
        if len(df_exons) > 0:
            last_exon = df_exons.iloc[-1]
            if last_exon['End'] + 1 not in [region[0] for region in intron_regions]:
                rows.append(last_exon)
        
        output_dfs.append(pd.DataFrame(rows))
    
    output_df = pd.concat(output_dfs, ignore_index=True)
    output_df = output_df[['Seqid', 'Source', 'Type', 'Start', 'End', 'Score', 'Strand', 'Frame', 'Suppl']]
    output_df.to_csv(gtf_file + '+introns.gtf', sep='\t', header=False, index=False)
    return output_df

def main(argv):
    input_gtf = argv[0]
    add_introns(input_gtf)

if __name__ == "__main__":
    main(sys.argv[1:])


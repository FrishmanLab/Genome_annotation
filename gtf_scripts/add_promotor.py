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


def add_promotors(gtf_file):
    df = gtf_to_df_with_genes(gtf_file)
    output_df = pd.DataFrame()  # Create an empty DataFrame to store the results

    for _, sub_df in df.groupby('Genes'):
        df_exons = sub_df[sub_df['Type'] == 'CDS'].sort_values('Start')
        if len(df_exons) > 0:
            if df_exons['Strand'].iloc[0] == '+':
                first_exon = df_exons.iloc[0]
                promotor_row = first_exon.copy()
                promotor_row['Type'] = 'promotor'
                if int(promotor_row['Start']) > 500:
                    promotor_row['End'] = int(promotor_row['Start']) + 2
                    promotor_row['Start'] = int(promotor_row['Start']) - 500
                    sub_df = sub_df.append(promotor_row)
                elif int(promotor_row['Start']) < 500 and int(promotor_row['Start']) > 10:
                    promotor_row['End'] = int(promotor_row['Start']) + 2
                    promotor_row['Start'] = 1 
                    sub_df = sub_df.append(promotor_row)
            elif df_exons['Strand'].iloc[0] == '-':
                first_exon = df_exons.iloc[-1]
                promotor_row = first_exon.copy()
                promotor_row['Type'] = 'promotor'
                promotor_row['Start'] = int(promotor_row['End']) - 2
                promotor_row['End'] = int(promotor_row['End']) + 500
                sub_df = sub_df.append(promotor_row)
                
        output_df = pd.concat([output_df, sub_df], ignore_index=True)

    output_df = output_df[['Seqid', 'Source', 'Type', 'Start', 'End', 'Score', 'Strand', 'Frame', 'Suppl']]
    output_df.to_csv(gtf_file + '+promotors.gtf', sep='\t', header=False, index=False)
    return output_df


def main(argv):
    input_gtf = argv[0]
    add_promotors(input_gtf)

if __name__ == "__main__":
    main(sys.argv[1:])


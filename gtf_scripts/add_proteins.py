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


def add_proteins(gtf_file):
    df = gtf_to_df_with_genes(gtf_file)
    output_df = pd.DataFrame()  # Create an empty DataFrame to store the results
    df_groups =  df.groupby("Genes")
    for name, sub_df in df_groups:
            df_cds = sub_df[sub_df['Type'] == 'CDS']
            df_cds = df_cds.sort_values('Start')
            strand = sub_df['Strand'].unique()[0]
            cds_regions = [(int(row['Start']), int(row['End'])) for _, row in df_cds.iterrows()]
            cds_regions.sort(key=lambda x: x[0])

            if cds_regions:
                start = cds_regions[0][0]
            else:
                start = None

            if cds_regions:
                end = cds_regions[-1][1]
            else:
                end = None

            if start is not None and end is not None:
                protein_row = pd.Series([sub_df['Seqid'].unique()[0], sub_df['Source'].unique()[0], 'protein', start, end, 0, strand, 0, sub_df['Suppl'].unique()[0], name], index=df.columns)
                sub_df.loc[len(sub_df)] = protein_row
            output_df = pd.concat([output_df, sub_df], ignore_index=True)

    output_df = output_df[['Seqid', 'Source', 'Type', 'Start', 'End', 'Score', 'Strand', 'Frame', 'Suppl']]
    output_df.to_csv(gtf_file + '+proteins.gtf', sep='\t', header=False, index=False)
    return output_df


def main(argv):
    input_gtf = argv[0]
    add_proteins(input_gtf)

if __name__ == "__main__":
    main(sys.argv[1:])



        

#!/usr/bin/env python

import sys
import pandas as pd


def gtf_to_df(gtf_file):
    column_names = ['Seqid', 'Source', 'Type', 'Start', 'End', 'Score', 'Strand', 'Frame', 'Suppl']
    df = pd.read_csv(gtf_file, sep='\t', index_col=False, names=column_names, dtype={'Start': int, 'End': int}, comment='#')
    return df

def add_non_promotors(gtf_file):
    df = gtf_to_df(gtf_file)
    output_dfs = []
    
    for _, sub_df in df.groupby('Seqid'):
        df_promotors = sub_df[sub_df['Type'] == 'promotor'].sort_values('Start')
        if len(df_promotors) > 1:
            end_values = df_promotors['End'].values
            start_values = df_promotors['Start'].values
            non_promotor_mask = start_values[1:] - end_values[:-1] > 1
            non_promotor_regions = list(zip(end_values[:-1][non_promotor_mask] + 1, start_values[1:][non_promotor_mask] - 1))
        else:
            non_promotor_regions = []
        
        rows = []
        for _, row in sub_df.iterrows():
            if row['Type'] == 'promotor':
                current_promotor = (row['Start'], row['End'])
                for non_promotor_region in non_promotor_regions:
                    if current_promotor[1] + 1 == non_promotor_region[0]:
                        new_row = row.copy()
                        new_row['Type'] = 'non_promotor'
                        new_row['Start'] = non_promotor_region[0]
                        new_row['End'] = non_promotor_region[1]
                        # rows.append(row)
                        rows.append(new_row)
            else:
                rows.append(row)
        
        if len(df_promotors) > 0:
            last_promotor = df_promotors.iloc[-1]
            if last_promotor['End'] + 1 not in [region[0] for region in non_promotor_regions]:
                rows.append(last_promotor)
        
        output_dfs.append(pd.DataFrame(rows))
    
    output_df = pd.concat(output_dfs, ignore_index=True)
    output_df = output_df[['Seqid', 'Source', 'Type', 'Start', 'End', 'Score', 'Strand', 'Frame', 'Suppl']]
    output_df.to_csv(gtf_file + '+non_promotors.gtf', sep='\t', header=False, index=False)
    return output_df

def main(argv):
    input_gtf = argv[0]
    add_non_promotors(input_gtf)

if __name__ == "__main__":
    main(sys.argv[1:])


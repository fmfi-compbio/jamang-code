import argparse
import pandas as pd
import numpy as np
from scipy.stats import fisher_exact

def main(treated_file, control_file):
    names = ['chr', 'pos', 'genome', 'reads']
    treated_df = pd.read_csv(treated_file, sep='\t', names=names, header=0)
    control_df = pd.read_csv(control_file, sep='\t', names=names, header=0)
    treated_sum = (treated_df['reads']/treated_df['genome']).sum()
    control_sum = (control_df['reads']/control_df['genome']).sum()

    joined_df = treated_df.merge(control_df, on=['chr', 'pos'], 
                                 suffixes=('_treated', '_control'))
    assert treated_df.shape[0] == control_df.shape[0] == joined_df.shape[0]    
    # print(joined_df.head())

    for (index, row) in joined_df.iterrows():
        treated_reads = row['reads_treated']
        control_reads = row['reads_control']
        treated_rest = treated_sum - treated_reads
        control_rest = control_sum - control_reads
        table = np.array([[treated_reads, control_reads],
                          [treated_rest, control_rest]])
        oddsratio, pvalue = fisher_exact(table, alternative='less')
        print(f'{row["chr"]}\t{row["pos"]}\t{pvalue}')            


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Fisher Test for k-mers')
    parser.add_argument('file1', type=str, help='Treated input filename')
    parser.add_argument('file2', type=str, help='Control input filename')
    args = parser.parse_args()

    main(args.file1, args.file2)

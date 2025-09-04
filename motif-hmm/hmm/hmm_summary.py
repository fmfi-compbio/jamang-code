#! /usr/bin/env python3

import pandas as pd
#import matplotlib.pyplot as plt
#import seaborn as sns
import sys

tab = pd.read_csv(sys.stdin, sep="\t", names=("chrom","read","group","value"))

tab2 = tab.groupby(['chrom','group'])['value'].describe().reset_index()
tab3 = tab2.loc[:, ['chrom','group','count','25%', '50%', '75%']]
tab3.to_csv(sys.stdout, sep="\t", index=False)

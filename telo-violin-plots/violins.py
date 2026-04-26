"""
This script creates violin plots of telomere lengths
for each species separately and for all four combined.
"""

import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('agg')
matplotlib.rcParams['font.size']=14
import matplotlib.pyplot as plt
import seaborn as sns

def get_data(filename):
    table = pd.read_csv(filename, sep="\t",
                        names=["telo", "read", "type", "len"])
    subst_func = lambda x : re.sub('e','R',re.sub('s','L',x))
    table['telo'] = table['telo'].apply(subst_func)
    table = table[table['type'].str.startswith('sub') | 
                  (table['type'] == 'telo')]
    grouped_table = (table.groupby(['telo', 'read'])['len']
             .sum().reset_index())

    #print(grouped_table.head())
    return grouped_table


def add_counts(table, ax, max_len, rotation=0):
    # add read counts on top of the plot in axes ax
    sizes = table['telo'].value_counts().sort_index()
    for (x_idx, chr) in enumerate(sizes.index):
        ax.text(x_idx, 0.99*max_len, str(sizes[chr]),
                ha="center", va="top", color='gray', rotation=rotation)
    ax.text(-0.55, 0.99*max_len, "n=",
            ha="right", va="top", color="gray")


def draw_data_single(table, species, filename_fig, filename_txt, max_len):
    # count chromosomes
    chrom = table['telo'].nunique() / 2
    
    fig, ax = plt.subplots(figsize=(chrom, 6))
    
    sns.violinplot(x='telo', y='len', data=table, ax=ax,
                   width=0.9, color='lightblue', 
                   bw_adjust=0.9, cut=0, density_norm="count")

    # read counts rotated by 30 degrees
    add_counts(table, ax, max_len, rotation=30)
    
    #ax.set_title(f"{species} telomeric read lengths")
    ax.set_ylabel('Repeat array length')
    ax.set_xlabel(None)
    ax.set_ylim(0, max_len)
    ax.set_xlim(left=-0.5)
    ax.tick_params("x", rotation=90)
    fig.savefig(filename_fig, bbox_inches='tight')

    # 3rd quantile in kb
    estimates = table.groupby('telo')['len'].quantile(0.75) / 1000
    # print with 1 decimal place
    estimates.to_csv(filename_txt,sep="\t", float_format='%.1f')
    
    
def draw_data_all(species, filename, max_len):
    # count chromosomes for each species
    chr_nums = []
    for (long, table) in species:
        chrom = table['telo'].nunique() / 2
        chr_nums.append(chrom)

    mosaic = \
    '''
    AAA
    BCD
    '''
    
    ratios = [chr_nums[1], chr_nums[2], chr_nums[3]]
    # alt mosaic: AAA.,BCDD
    #ratio4 = chr_nums[1]+chr_nums[2]+chr_nums[3]-chr_nums[0]    
    #ratios = [chr_nums[1], chr_nums[2], chr_nums[3]-ratio4, ratio4]

    fig, axes_dict =plt.subplot_mosaic(
        mosaic,
        width_ratios = ratios,
        figsize=(chr_nums[0], 10),
        layout="constrained",
        sharey=True)

    axes = [axes_dict[x] for x in "ABCD"]
    
    for (i, ax) in enumerate(axes):
                           
        sns.violinplot(x='telo', y='len', data=species[i][1], ax=ax,
                       width=0.9, color='lightblue',
                       bw_adjust=0.9, cut=0, density_norm="count")

        # read counts rotated by 30 degrees
        add_counts(species[i][1], ax, max_len, rotation=30)
        
        # other plot settings    
        ax.set_title(f"{species[i][0]}", style='italic')
        ax.set_xlabel(None)
        ax.set_ylabel(None)
        ax.set_xlim(left=-0.5)
        ax.set_ylim(0, max_len)
        ax.tick_params("x", rotation=90)

    for i in [0, 1]:
        axes[i].set_ylabel('Repeat array length')
    fig.savefig(filename, bbox_inches='tight')

    

species = [("jamAng", "J. angkorensis", 35000),
           ("jamPal", "J. pallidilutea", 20000),
           ("jamPhy", "P. phylloscopi", 25000),
           ("symKan", "S. kandeliae", 48000)]

# with just long name and table of data
species2 = []

for (short, long, max_len) in species:
    table = get_data(f"{short}-counts.tsv")
    species2.append((long, table))
    draw_data_single(table, long,
                     f"{short}-counts.pdf", f"{short}-counts.txt", max_len)

draw_data_all(species2, "all-counts.pdf", 48000)

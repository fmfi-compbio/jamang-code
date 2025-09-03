# avoid initialization of graphics winows (before import plt!)
import matplotlib
matplotlib.use('cairo')   # also possibly 'agg' for png or 'pdf'

import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


import colorcet as cc
chrom_pal = sns.color_palette(cc.glasbey, n_colors=16)

def chr2color(chr, all_chr2):
  assert chr in all_chr2
  idx = all_chr2.index(chr)
  return chrom_pal[idx]

def get_num(chr, sizes):
  for (idx, row) in sizes.iterrows():
    if chr == row.loc["chr"]:
      return idx
  return None

def get_y(chr_name, sp_name, sizes, sp_names):
  chr_idx = get_num(chr_name, sizes)
  sp_idx = sp_names.index(sp_name)
  return -chr_idx*(len(sp_names)+0.5)-sp_idx

def draw_alignments(pslview2, sizes, all_chr2, sp_names, axes, bottom, height):
    for (idx, row) in pslview2.iterrows():
      y = get_y(row.loc['chr1'], row.loc['sp2'], sizes, sp_names)
      width = row.loc["end1"]-row.loc["start1"]
      color = chr2color(row.loc["chr2"], all_chr2)
      axes.add_patch(plt.Rectangle((row["start1"], y+bottom),
                                 width, height, fc=color, ec="None"))

def main(sizes_file, pslview2_files, sp_names, output_file):
  columns_pslview2 = ['score', 'strand', 'chr1', 'len1', 'start1', 'end1',
                      'chr2', 'len2', 'start2', 'end2', 'pid']
  columns_sizes = ['chr', 'length']
  sizes = pd.read_csv(sizes_file, sep="\t",  names=columns_sizes)
  sizes = sizes.query("chr != 'mtDNA'")
  max_chr = sizes["length"].max()

  pslview2 = None
  for (sp_idx, sp_name) in enumerate(sp_names):
    sp_pslview2 = pd.read_csv(pslview2_files[sp_idx],
                              sep="\t",  names=columns_pslview2)
    sp_pslview2['sp2'] = sp_name
    #add sp_pslview2 to pslview2 dataFrame
    if pslview2 is None:
      pslview2 = sp_pslview2
    else:
      pslview2 = pd.concat([pslview2, sp_pslview2], ignore_index=True)
  
  # skip mtDNA
  pslview2 = pslview2.query("chr1 != 'mtDNA' and chr2 != 'mtDNA'")
  # rename chr1 to chr01 etc 
  print("Renaming chr1 to chr01 etc")
  new_chr2 = pslview2['chr2'].replace(re.compile(r'^chr(\d)$'), r'chr0\1')
  print("having new_chr2:", new_chr2.unique())
  pslview2['chr2'] = new_chr2
  print("done renaming chr1 to chr01 etc")

  all_chr2 = sorted(pslview2.chr2.unique())
  lowest_y = get_y(sizes.iloc[-1,0], sp_names[-1], sizes, sp_names)
  
  scale = 2
  figure, axes = plt.subplots(1, 1,
                              #layout='constrained',
                              figsize=(scale * max_chr/1e6*2, scale * len(sizes)/5))
  axes.axis('off')
  max_x = max_chr * 1.001
  axes.set_xlim(0, max_x)
  axes.set_ylim(lowest_y-1, 1)

  # chrom
  bottom = 0.2
  height = 0.6

  draw_alignments(pslview2, sizes, all_chr2, sp_names, axes, bottom, height)
  
  text_gap = 1e4
  for (chr_idx, row) in sizes.iterrows():
    chr_name = row.loc['chr']
    size = row.loc['length']
    #text = chr_name + f" {size/1e6:.1f}Mbp"
    y = get_y(chr_name, sp_names[1], sizes, sp_names)
    axes.text(-text_gap, y+0.5, chr_name, 
              ha='right', va='center', fontdict={"size":8})
    for (sp_idx, sp_name) in enumerate(sp_names):
         
      y = get_y(chr_name, sp_name, sizes, sp_names)      
      axes.add_patch(plt.Rectangle((0, y+bottom),
                                   size, height, ec="black", fill=False, 
                                   linewidth=0.3))
      axes.text(size+text_gap, y+0.5, sp_name, 
                ha='left', va='center', fontdict={"size":5})


  start_y = -10
  start_x = 0.8 * max_x
  x_size = 1e4
  for (idx, chr_name) in enumerate(all_chr2):
    y = start_y - idx
    color = chr2color(chr_name, all_chr2)    
    axes.add_patch(plt.Rectangle((start_x, y+bottom),
                                  x_size, height, ec=color, fc=color))    
    axes.text(start_x+x_size+text_gap, y+0.5, chr_name, 
              ha='left', va='center', fontdict={"size":5})

  figure.savefig(output_file, bbox_inches='tight')



assemblies = ["jamPalA1", "jamPhyA2", "symKanA1"]
pslview2_files = [f"{sp}-LASTSPLIT-jamAngA5.psl.view2" for sp in assemblies]
species_names = [sp[3:5] for sp in assemblies]

main("jamAngA5.sizes", pslview2_files, species_names, "jamAng-chrompaint.pdf")


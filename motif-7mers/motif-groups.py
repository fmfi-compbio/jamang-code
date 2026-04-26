import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from collections import defaultdict
import numpy as np

# read matrix with motif position
a = pd.read_csv('groups-manual-matrix9.csv', header=None).reset_index().rename(columns={'index':'row'})
# make into a long table with columns col, row and name
# drop missing values, i.e. empty fields
a = a.melt(id_vars='row', var_name='col', value_name='name').dropna().set_index('name', verify_integrity=True)
print(a)

# read other file with groups
levels = ['0','0.1','0.7','0.8']
df = pd.read_csv('groups.tsv', sep='\t', names=['name']+levels)
df['species'] = df['name'].apply(lambda x: x[0:2])
df.set_index('name', inplace=True)
df['row'] = a['row']
df['col'] = a['col']
# count singletons (all rows in species minus those with row value)
# i.e. included in the figure
missing = df.groupby('species')['row'].size() \
          - df.groupby('species')['row'].count()
print("missing", missing)
# now drop those without row/column values
df.dropna(inplace=True)
df.reset_index(inplace=True)
print(df.head()) 


# for each rectangle the smallest group with that rectangle
rect2group = dict()
# for each group its rectangle, including redundant ones
group2rect = dict()
# list of groups to draw in reversed order
to_draw = []

# for all levels of groups
for level in levels:
    groups = df[level].unique()
    row_min = df.groupby(level)['row'].min()
    row_max = df.groupby(level)['row'].max()
    col_min = df.groupby(level)['col'].min()
    col_max = df.groupby(level)['col'].max()
    for group in groups:
        rect = (row_min[group], row_max[group], col_min[group], col_max[group])
        group_id = (level, group)
        group2rect[group_id] = rect
        # if we found a new rectangle, add it to structures
        if rect not in rect2group:
            rect2group[rect] = group_id
            to_draw.append(group_id)
        # TODO: check that all in rectangle is part of group

# revrse order of the list to start with least significant levels
to_draw.reverse()

def get_species(df, level, group):
    sel = df[df[level]==group]
    assert sel.shape[0]>0
    species = sel['species'].unique()    
    if len(species)==1:
        return species[0]
    else:
        return None


level_width = {'0':0, '0.1':0, '0.7':0, '0.8':1}
#level_style = ['-', '-', '--']
level_alpha = {'0':0.8, '0.1':0.5, '0.7':0.3, '0.8':0}
level_gap = {'0':0, '0.1':1, '0.7':2, '0.8':3}
    
def plot_rectangles(axes, to_draw, group2rect, df, colors, gapx, gapy):
    
    for group_id in to_draw:
        level = group_id[0]
        rect = group2rect[group_id]
        sp = get_species(df, group_id[0], group_id[1])
        if sp is not None:
            color = colors[sp]
        else:
            color = "black"
        gap = level_gap[level]/5
        gapxh = gapx*gap
        gapyh = gapy*gap
        x = rect[2]-gapxh
        y = -rect[1]-gapyh
        sizex = rect[3]-rect[2]+1-gapx+2*gapxh
        sizey = rect[1]-rect[0]+1-gapy+2*gapyh
        linewidth=level_width[level]
        linestyle='-' #level_style[level]
        alpha = level_alpha[level]

        #if sp != 'Ph':
        #    continue

        #print(group_id, rect, x, y, sizex, sizey, 'alpha', alpha, linewidth)
        
        if alpha > 0 :
            axes.add_patch(plt.Rectangle((x+gapx, y+gapy), sizex, sizey,
                                         capstyle='round',
                                         alpha=alpha, facecolor=color,
                                                linestyle=None
                                         ))
        if linewidth > 0:
            axes.add_patch(plt.Rectangle((x+gapx, y+gapy), sizex, sizey,
                                         linewidth=linewidth,
                                         edgecolor=color,
                                         fill=False,
                                         linestyle=linestyle
                                         ))
            

    return
    

# finally the full figure

#colors = {'Pa': 'red', 'Ph': 'blue', 'An': 'green', 'Ka': 'orange'}
colors = {'Pa': 'C0', 'Ph': 'C4', 'An': 'C1', 'Ka': 'C2'}
sp2long = {'Pa': 'J. pallidilutea', 'Ph': 'P. phylloscopi', 'An': 'J. angkorensis', 'Ka': 'S. kandeliae'}
(figure, axes) = plt.subplots(figsize=(11, 8))
max_row = df['row'].max()
max_col = df['col'].max()
gapx = 0.15
gapy = 0.31
sizex = 1-2*gapx
sizey = 1-2*gapy
offsetx = 0.02+0.5*sizex
offsety = 0.015
plot_rectangles(axes, to_draw, group2rect, df, colors, gapx, gapy)
for i in range(len(df)):
    x = df.iloc[i]['col']
    y = -df.iloc[i]['row']
    name = df.iloc[i]['name']
    sp = name[0:2]
    color = colors[sp]
    #axes.add_patch(plt.Rectangle((x+gapx, y+gapy), sizex, sizey, alpha=0.3, facecolor=color))
    axes.text(x+gapx+offsetx, y+gapy+offsety, name[2:],
              fontdict={'fontsize': 14}, va='bottom', ha='center')

    
axes.set_axis_off()
axes.set_xlim(-gapx, max_col+1+gapx)
axes.set_ylim(-max_row-gapy, 1+gapy)
figure.savefig('groups.pdf', bbox_inches='tight')


# plotting legend
rec_size = 0.9
(figure, axes) = plt.subplots(figsize=(4, 1.5))
for (idx1, level) in enumerate(levels):    
        linewidth=level_width[level]
        linestyle='-' #level_style[level]
        alpha = level_alpha[level]

        level_text = level
        if level == '0':
            level_text = 'id.'
        axes.text(idx1+0.5, 4.1, level_text, ha='center',va='bottom')
        

        for (idx2, species) in enumerate(colors.keys()):
            color = colors[species]


            print(idx1, idx2, level, species, color)
        
            if alpha > 0 :
                axes.add_patch(plt.Rectangle((idx1, idx2), rec_size, rec_size,
                                             capstyle='round',
                                             alpha=alpha, facecolor=color,
                                             linestyle=None
                                             ))
            if linewidth > 0:
                axes.add_patch(plt.Rectangle((idx1, idx2), rec_size, rec_size,
                                             linewidth=linewidth,
                                             edgecolor=color,
                                             fill=False,
                                             capstyle='round',
                                             linestyle=linestyle
                                             ))

for (idx2, species) in enumerate(colors.keys()):
    axes.text(4.2, idx2+0.5, sp2long[species],va='center',ha='left',style='italic')
axes.text(2.5, 4.7, 'distance', ha='center',va='bottom')     
            
axes.set_xlim(0,6)
axes.set_ylim(0,5)


    
#    level_value = df.columns[level_index+1]
#    axes.plot([], [], linewidth=level_width[level_index], color='black', 
#              linestyle=level_style[level_index], label=f'{level_value}')
#axes.legend(loc='upper left')
axes.set_axis_off()
figure.savefig('groups-legend.pdf', bbox_inches='tight')


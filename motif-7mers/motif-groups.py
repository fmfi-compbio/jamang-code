import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from collections import defaultdict
import numpy as np

# read matrix with motif position
a = pd.read_csv('groups-manual-matrix12.csv', header=None).reset_index().rename(columns={'index':'row'})
a = a.melt(id_vars='row', var_name='col', value_name='name').set_index('name', verify_integrity=True)
print(a)

# read other file with groups
df = pd.read_csv('groups-manual.tsv', sep='\t', names=['name','0.1','0.7','0.8','row', 'col'])
df['species'] = df['name'].apply(lambda x: x[0:2])
df.set_index('name', inplace=True)
df['row'] = a['row']
df['col'] = a['col']
df.reset_index(inplace=True)
print(df.head()) 

# find groups with only one member
# and also groups with a unique species
singletons = set()
species = dict()
for i in range(1, 4):
    two_cols = df.iloc[:, [i-1, i]].drop_duplicates()
    counts = two_cols.iloc[:,1].value_counts()
    for group in counts[counts == 1].index:
        singletons.add((i-1, group))

    two_cols = df.iloc[:, [i, 6]].drop_duplicates()
    counts = two_cols.iloc[:,0].value_counts()
    good_groups = counts[counts == 1].index
    two_cols.set_index(two_cols.columns[0], inplace=True)
    for group in good_groups:
        rows = two_cols.loc[group]
        species[(i-1, group)] = rows.iloc[0]

# find and simplify boundaries of groups at three levels
def straight(a, b , c):
    if a[0]==b[0] and b[0]==c[0]:
        return True
    if a[1]==b[1] and b[1]==c[1]:
        return True
    return False

def simplify(differences):    
    for chain_id in differences:
        chain = differences[chain_id]
        start = min(chain.keys())
        current = start
        count = 0
        while True:
            next = chain[current]
            next2 = chain[next]
            assert current != next and current != next2
            if straight(current, next, next2):
                chain[current] = next2
                del chain[next]
            else:
                current = next
                count += 1
                if current == start:
                    break
        assert count == len(chain.keys()), f'{count} != {len(chain.keys())} for {chain_id}\n{chain}'
    return differences

def check_neighbors(df, singletons):
    def get_group(row_index, col_index, a, df_indexed, level):
        if row_index < 0 or row_index >= a.shape[0] or col_index < 0 or col_index >= a.shape[1]:
            return -1
        name = a.iloc[row_index, col_index]
        if name not in df_indexed.index:
            return -1
        return df_indexed.loc[name][level]
    
    # a 2d matrix of group names
    a = df.pivot(columns='col', index='row', values='name')
    #a.to_csv('a.csv',index=False, header=False)
    #display(a)
    # df indexed by name
    df_indexed = df.set_index('name')
    #display(df_indexed.head(2))
    vectors = [(-1,0), (1,0), (0,-1), (0,1)]  # row and col offsets
    vector_corners = {(-1, 0): [(-1, 1), (-1, 0)], # up 
                      (1, 0): [(0, 0), (0, 1)], # down
                      (0, -1): [(-1, 0), (0, 0)], # left
                      (0, 1): [(0, 1), (-1, 1)]} # right
    differences = defaultdict(dict)
    for i in range(len(df)):
        row = df.iloc[i]
        name = row['name']
        row_index = row['row']
        col_index = row['col']
        for vector in vectors:
            new_row_index = row_index + vector[0]
            new_col_index = col_index + vector[1]
            for (level_pos, level) in enumerate(['0.1', '0.7', '0.8']):
                group1 = df_indexed.loc[name][level]
                if (level_pos, group1) in singletons:
                    continue
                group2 = get_group(new_row_index, new_col_index, a, df_indexed, level)
                if group1 != group2:
                    corners = vector_corners[vector]
                    start = (row_index+corners[0][0], col_index+corners[0][1])
                    end = (row_index+corners[1][0], col_index+corners[1][1])
                    differences[(level_pos, group1)][start] = end
    return simplify(differences)


differences = check_neighbors(df, singletons)


# functions for plotting lines
def get_coords(start, end, next, level, gapx, gapy):
    def get_sign(a, b):
        return (np.sign(b[0]-a[0]), np.sign(b[1]-a[1]))

    signs = (get_sign(start, end), get_sign(end, next))
    sign_vectors = {
        ((-1,0),(0,-1)): (1,-1), #up left
        ((1,0),(0,1)): (-1,1), #down right
        ((0,-1),(1,0)): (1,1), #left down
        ((0,1),(-1,0)): (-1,-1), #right up
        ((-1,0),(0,1)): (-1,-1), #up right
        ((1,0),(0,-1)): (1,1), #down left
        ((0,-1),(-1,0)): (1,-1), #left up
        ((0,1),(1,0)): (-1,1) #right down
    }


    if signs in sign_vectors:
        (d_row, d_col) = sign_vectors[signs]
    else:
        (d_row, d_col) = (0,0)

    gapx2 = gapx - gapx * (level+1) / 3.8
    gapy2 = gapy - gapy * (level+1) / 3.8
    

    row = end[0] + d_row * gapy2
    col = end[1] + d_col * gapx2
    return (col, -row)
    


level_width = [3, 1, 2]
level_style = ['-', '-', '--']
            
def plot_lines(axes, differences, species, colors, gapx, gapy):
    
    for (level, group) in differences:
        chain = differences[(level, group)]
        for prev in chain:
            start = chain[prev]
            end = chain[start]
            next = chain[end]
            x1, y1 = get_coords(prev, start, end, level, gapx, gapy)
            x2, y2 = get_coords(start, end, next, level, gapx, gapy)
            if (level, group) in species:
                color = colors[species[(level, group)]]
            else:
                color = 'black'
            axes.plot([x1, x2], [y1, y2], linewidth=level_width[level], color=color, alpha=0.8, linestyle=level_style[level])


# finally the full figure

#colors = {'Pa': 'red', 'Ph': 'blue', 'An': 'green', 'Ka': 'orange'}
colors = {'Pa': 'C0', 'Ph': 'C4', 'An': 'C1', 'Ka': 'C2'}
(figure, axes) = plt.subplots(figsize=(18, 8))
max_row = df['row'].max()
max_col = df['col'].max()
gapx = 0.13
gapy = 0.31
sizex = 1-2*gapx
sizey = 1-2*gapy
offsetx = 0.01
offsety = 0.015
plot_lines(axes, differences, species, colors, gapx, gapy)
for i in range(len(df)):
    x = df.iloc[i]['col']
    y = -df.iloc[i]['row']
    name = df.iloc[i]['name']
    sp = name[0:2]
    color = colors[sp]
    axes.add_patch(plt.Rectangle((x+gapx, y+gapy), sizex, sizey, alpha=0.3, facecolor=color))
    axes.text(x+gapx+offsetx, y+gapy+offsety, name, 
              fontdict={'fontsize': 14}, va='bottom', ha='left')

    
axes.set_axis_off()
axes.set_xlim(-gapx, max_col+1+gapx)
axes.set_ylim(-max_row-gapy, 1+gapy)
figure.savefig('groups.pdf', bbox_inches='tight')


# plotting legend
(figure, axes) = plt.subplots(figsize=(1, 1))
for level_index in range(3):
    level_value = df.columns[level_index+1]
    axes.plot([], [], linewidth=level_width[level_index], color='black', 
              linestyle=level_style[level_index], label=f'{level_value}')
axes.legend(loc='upper left')
axes.set_axis_off()
figure.savefig('groups-legend.pdf', bbox_inches='tight')


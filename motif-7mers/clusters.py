import argparse
import math
import re

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import seaborn as sns

import numpy as np
from scipy.cluster.hierarchy import linkage, dendrogram 
from scipy.spatial.distance import squareform
import pandas as pd
import sklearn.manifold as manifold
from sklearn.cluster import AgglomerativeClustering


def main():
    parser = argparse.ArgumentParser(description='Create clusters.')
    parser.add_argument('input_filename', type=str,
                        help='The input fasta file name')
    parser.add_argument('output_prefix', type=str,
                        help='The output file prefix')
    parser.add_argument('--containment', action='store_true',
                        help='Enable containment distance, otherwise Jaccard')
    parser.add_argument('--tsne', action='store_true',
                        help='Enable TSNE embedding, otherwise MDS')
    parser.add_argument('--threshold', type=float, default=0.5,
                        help='Threshold for clustering visualization (multiple of 0.1')
    parser.add_argument('--kmer', type=int,
                        default=7, help='k-mer size')
    parser.add_argument('--dist_list', 
                        default=None,
                        help='distances to use, se[arated by ",".Default all with step 0.1')

    
    args = parser.parse_args()

    if args.dist_list is not None:
        tmp_list = args.dist_list.split(",")
        args.dist_list = [float(x) for x in tmp_list]
    
    motifs = read_fasta(args.input_filename, args.kmer)
    names = list(motifs.keys())

    dist_matrix = get_dist_matrix(motifs, names, args.containment)

    dist_matrix_named = pd.DataFrame(dist_matrix, index=names, columns=names)
    dist_matrix_named.to_csv(f"{args.output_prefix}-dist.tsv", sep='\t')
    
    #np.savetxt(f"{args.output_prefix}-dist-plain.tsv", dist_matrix)
    #np.savetxt(f"{args.output_prefix}-names.tsv",
    #           np.array(names), fmt="%s")
    #dist_matrix = np.loadtxt("dist-plain.tsv")
    #names = np.genfromtxt('names.tsv',dtype='str')

    dendogram(dist_matrix, names, args.threshold,
              f"{args.output_prefix}-dendogram.pdf")

    table = get_embedding(dist_matrix, names, args.tsne)
    group_keys = add_distances(table, dist_matrix, args.dist_list)
    table['species'] = table['name'].apply(lambda x: x.split('-')[0])

    # write just groups
    table2 = table.drop(columns=['x','y','species'])
    table2.to_csv(f"{args.output_prefix}-groups.tsv", sep="\t", index=False)

    # write embedding to file
    which_group = f"group{args.threshold:.1f}"
    assert which_group in group_keys
    table3 = (table.loc[:,['name','species','x','y',which_group]]
              .rename(columns={which_group:'group'}))
    table3.to_csv(f"{args.output_prefix}-embed.tsv", sep="\t", index=False)
    
    write_groupings(table, group_keys,
                    f"{args.output_prefix}-groups.txt",
                    f"{args.output_prefix}-groups-sp.txt")

    draw_embed(table3, f"{args.output_prefix}-embed.pdf")

def draw_embed(table, output_filename):
    group_counts = table['group'].value_counts()
    big_groups = group_counts[group_counts >= 2].index.tolist()

    eps = 0.01
    group_min = table.groupby('group')[['x', 'y']].min().rename(columns={'x': 'x_min', 'y': 'y_min'})
    group_min -= eps
    group_max = table.groupby('group')[['x', 'y']].max().rename(columns={'x': 'x_max', 'y': 'y_max'})
    group_max += eps
    group_coords = pd.concat([group_min, group_max], axis=1)
    group_coords = group_coords.loc[big_groups]
    # Create a scatter plot
    plt.figure(figsize=(10, 6))
    sns.scatterplot(data=table, x='x', y='y', hue='species', s=20)

    for group in big_groups:
        plt.gca().add_patch(plt.Rectangle(
            (group_coords.loc[group, 'x_min'], group_coords.loc[group, 'y_min']),
            group_coords.loc[group, 'x_max'] - group_coords.loc[group, 'x_min'],
            group_coords.loc[group, 'y_max'] - group_coords.loc[group, 'y_min'],
            edgecolor='gray', facecolor='none'
        ))

        # Add labels to each point
        #for i in range(table.shape[0]):
        #    plt.text(table['x'][i], table['y'][i], table['name'][i], fontsize=9)
    plt.savefig(output_filename)

    
def write_groupings(table, group_keys, filename1, filename2):
    # write all non-trivial groupings to files
    # (both species-level and indivodual names)

    file1 = open(filename1, 'w')
    file2 = open(filename2, 'w')

    for group_key in group_keys:
        print(f"\n### {group_key}", file=file1)
        print(f"\n### {group_key}", file=file2)
        group_counts = table[group_key].value_counts()
        big_groups = group_counts[group_counts >= 2].index.tolist()
        
        for group in big_groups:
            group_sp = table[table[group_key] == group]['species'].value_counts()
            print(f"{group}:", end='', file=file2)
            for sp in group_sp.index:
                print('', sp, group_sp[sp], end='', file=file2)
            print(file=file2)

            group_names = table[table[group_key] == group]['name'].sort_values()
            to_print = ' '.join(group_names)
            print(f"{group}: {to_print}", file=file1)    

    file1.close()
    file2.close()
    


def dendogram(dist_matrix, names, threshold, output_filename):
    clustering = linkage(squareform(dist_matrix), method='average')
    plt.figure()
    dendrogram(clustering, labels=names,
               orientation='right',leaf_font_size=2)
    plt.savefig(output_filename)

    
def get_dist_matrix(motifs, names, containment):
    n = len(names)
    dist_matrix = np.zeros((n,n))
    
    for i in range(n):
        for j in range(i+1,n):
            dist_matrix[i,j]=jaccard(motifs[names[i]],
                                     motifs[names[j]],
                                     containment)
            dist_matrix[j,i]=dist_matrix[i,j]
    return dist_matrix

def rewrite_name(name):
    # currently not used
    #jamPal.sub09s1 -> jamPal-chr09L-p1
    
    # sub and telo to d and p
    name = re.sub(r'sub','p',name)
    name = re.sub(r'telo','d',name)
    # s and e to L and R
    name = re.sub(r'([0-9])s',r'\1L',name) 
    name = re.sub(r'([0-9])e',r'\1R',name) 

    # reorder and change to dashes
    name = re.sub(r'^(\w+)\.([a-z])(\d+)([A-Z])(.*)$',
                  r'\1-chr\3\4-\2\5', name)
    
    return name

def read_fasta(file, kmer):
    '''Reads a fasta file and returns a dictionary of sorted k-mer lists'''
    result = dict()
    with open(file) as f:
        for line in f:
            line = line.strip()

            if line[0] == '>':
                name = line[1:]
                #name = rewrite_name(line[1:])
                assert name not in result
                seq = ""
            else:
                seq += line
                telo = seq + seq[0:kmer-1]
                kmers = [telo[i:i+kmer]
                         for i in range(len(telo)-kmer+1)]
                result[name] = sorted(kmers)
    return result

def jaccard(motif1, motif2, containment=False):
    '''Returns the Jaccard similarity between two sorted multisets'''
    m1 = len(motif1)
    m2 = len(motif2)
    i1 = 0
    i2 = 0

    eps = 0.001
    intersection = 0
    union = 0
    while i1 < m1 and i2 < m2:
        if motif1[i1] == motif2[i2]:
            i1 += 1
            i2 += 1
            intersection += 1
            union += 1
        elif motif1[i1] < motif2[i2]:
            i1 += 1
            union += 1
        else:
            i2 += 1
            union += 1
    
    union += (m1 - i1) + (m2 - i2)
    smaller = min(m1, m2)

    # jaccard or symmetric containment
    if containment:
        return 1-intersection/smaller
    else:
        return 1-intersection/union




def get_embedding(dist_matrix, names, tsne=False):

    if tsne:
        embedded = manifold.TSNE(n_components=2, init='random',
                                 metric='precomputed',
                                 perplexity=5).fit_transform(dist_matrix)
    else:        
        embedded = manifold.MDS(
            n_components=2,
            random_state=0,
            dissimilarity='precomputed',
        ).fit_transform(dist_matrix)

    table = pd.DataFrame({'name':names,
                      'x':embedded[:,0],
                      'y':embedded[:,1]})
    return table


def add_distances(table, dist_matrix, dist_list=None):

    eps = 0.0001
    if dist_list is None:
        dist_list = list(np.arange(0,1.05,0.1))

    keys = []
    for d in dist_list:
        clustering = AgglomerativeClustering(
            linkage='single',
            n_clusters=None,
            distance_threshold=d+eps,
            metric='precomputed').fit(dist_matrix)

        # relabel clusters to have bigger groups first
        labels = pd.Series(clustering.labels_)
        counts = labels.value_counts().reset_index()
        replace = {counts.iloc[x]['index']:x for x in counts.index}
        new_labels = labels.replace(replace)
                    
        #np.savetxt("clusters.txt", clustering.labels_, fmt="%d")
        group_name = f"group{d:.1f}"
        table[group_name] = new_labels
        keys.append(group_name)

    
        table.sort_values(list(reversed(keys)), inplace=True)

        def relabel_first(series):
            replace = dict()
            prev = -1
            count = 1
            for val in series:
                if val != prev:
                    assert val not in replace, "val"
                    replace[val] = count
                    count += 1
                    prev = val
            return series.replace(replace)

        def relabel_next(table, col, prev):
            mins = table.groupby(col)[prev].min()
            replace = dict()
            for orig in mins.index:
                replace[orig] = mins[orig]
            return table[col].replace(replace)


        # relabel again group 0 to be 0,1,...
        table[keys[0]] = relabel_first(table[keys[0]])
        # relabel other groups to get minimum of previous group ID's
        # plus check containment
        for i in range(1,len(keys)):
            table[keys[i]] = relabel_next(table, keys[i], keys[i-1])
    return keys
            

if __name__ == '__main__':
    main()

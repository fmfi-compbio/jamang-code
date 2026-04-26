import pandas as pd
import matplotlib
matplotlib.use('agg')
matplotlib.rcParams['font.size']=12
import matplotlib.pyplot as plt
import seaborn as sns


genomes = ["jamAng", "jamPal", "jamPhy", "symKan", "ustMay"]
long = ["J. angkorensis", "J. pallidilutea", "P. phylloscopi", "S. kandeliae", "M. maydis"]

fig, axes = plt.subplots(2, 3, figsize=(15, 6), sharey=True, sharex=True)
axes1D = axes.flatten()
category = pd.CategoricalDtype(categories=["telomeric", "other"], ordered=True)
for i, genome in enumerate(genomes):
    df = pd.read_csv(f"{genome}-tandem-repeats-proc.tsv", sep="\t", 
                     names=["id", "rep", "count", "length", "telo"])
    # in column 'telo' change 0 to "other" and 1 to "telo"
    df["telo"].replace({0: "other", 1: "telomeric"}, inplace=True)
    # change type of "telo" to category with order 
    df["telo"] = df["telo"].astype(category, copy=True)
    # print(df.head(), genome)

    sns.scatterplot(data=df, x="length", y="rep", hue="telo", ax=axes1D[i], s=10, 
                    palette=["red", "gray"], alpha=0.8, legend=(i==4))
    axes1D[i].set_title(long[i], style='italic')
    axes1D[i].set_xscale("log")
    axes1D[i].set_xlabel("total length")
    axes1D[i].set_ylabel("motif length")
    
# change legend title to "repeat type"
axes1D[4].legend(title="repeat type", bbox_to_anchor=(0.735, 0.5, 0.9, 0.5))
axes1D[5].set_axis_off()
#fig.legend()

        
    
fig.savefig("tandem-repeats.pdf", bbox_inches='tight')



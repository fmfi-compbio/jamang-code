import pandas as pd
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import seaborn as sns

genomes = ["jamAng", "jamPal", "jamPhy", "symKan", "ustMay"]
long = ["J. angkorensis", "J. pallidilutea", "J. phylloscopi", "S. kandeliae", "M. maydis"]

fig, axes = plt.subplots(1, 5, figsize=(20, 4), sharey=True, sharex=True)
category = pd.CategoricalDtype(categories=["telomeric", "other"], ordered=True)
for i, genome in enumerate(genomes):
    df = pd.read_csv(f"{genome}-tandem-repeats-proc.tsv", sep="\t", 
                     names=["id", "rep", "count", "length", "telo"])
    # in column 'telo' change 0 to "other" and 1 to "telo"
    df["telo"].replace({0: "other", 1: "telomeric"}, inplace=True)
    # change type of "telo" to category with order 
    df["telo"] = df["telo"].astype(category, copy=True)
    # print(df.head(), genome)

    sns.scatterplot(data=df, x="length", y="rep", hue="telo", ax=axes[i], s=10, 
                    palette=["red", "gray"], alpha=0.8, legend=(i==0))
    axes[i].set_title(long[i])
    axes[i].set_xscale("log")
    axes[i].set_xlabel("total length")
    axes[i].set_ylabel("motif length")
    # change legend title to "tandem repeat type"
    if i==0:
        axes[i].legend(title="repeat type", loc="upper right")
    
fig.savefig("tandem-repeats.pdf", bbox_inches='tight')



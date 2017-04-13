__author__ = 'jgwall'

import argparse
import biom
import math
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re

debug = False


def main():
    args = parse_args()
    print("Plotting heritabilities from", args.infile)

    # Load and sort data
    data = pd.read_csv(args.infile, sep='\t')
    order = np.argsort(data.loc['actual', :])[::-1]
    data = data.iloc[:, order]

    # Subset
    if args.top_n is not None:
        print("\tSubsetting to just the top", args.top_n, "heritable taxa")
        data = data.iloc[:, :args.top_n]

    # Load biom data
    table = biom.load_table(args.biom)
    ids = table.ids(axis='observation')
    metadata = table.metadata(axis='observation')
    taxonomy = [m['taxonomy'] for m in metadata]
    otu_key = make_otu_key(ids, taxonomy)

    # Plot
    fig = plt.figure(figsize=(3 + .08 * len(data.columns), 20))
    ax_top    = fig.add_axes([0.05, 0.7, 0.92, 0.27], ylabel="Heritability (H$^2$)")   #TODO: Change y-values & figure size
    ax_bottom = fig.add_axes([0.05, 0.2, 0.92, 0.27], ylabel="Heritability (H$^2$)")

    #PLot
    split = math.ceil(len(data.columns)/2)
    plot_herits(ax_top, data.iloc[:, :split], args, otu_key)
    plot_herits(ax_bottom, data.iloc[:, split:], args, otu_key)

    # Save
    fig.savefig(args.outprefix + ".png", dpi=100)
    fig.savefig(args.outprefix + ".svg", dpi=600)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile")
    parser.add_argument("-o", "--outprefix")
    parser.add_argument("-b", "--biom", help="Biom file with taxonomic data for each OTU")
    parser.add_argument("-p", "--p-cutoff", type=float, default=0.001)
    parser.add_argument("-n", "--top-n", type=int, help="Only plot the top n heritable OTUs")
    parser.add_argument("--debug", default=False, action="store_true")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug

    return parser.parse_args()


def plot_herits(ax, data, args, otu_key):
    # Violin plots of random permutations
    perms = data.index != "actual"
    xticks, xlabels = list(), list()
    violins = list()
    for x in range(len(data.columns)):
        trait = data.columns[x]
        xticks.append(x)
        tmp_violin = ax.violinplot(data.loc[perms, trait], positions=[x])
        violins.append(tmp_violin)
        new_label = otu_key[trait.replace("trait_", "")]
        xlabels.append(new_label)

    # Add dots for actual heritabilities
    actuals = np.array(data.loc['actual', :])
    perms = np.array(data.iloc[1:, :])
    pvals = [np.nan] * len(actuals)
    for i in range(len(actuals)):
        pvals[i] = np.sum(perms[:, i] >= actuals[i]) / len(perms)
    colors = np.array(['red' if p <= args.p_cutoff else 'darkgray' for p in pvals])
    print("\t", sum(colors == 'red'), "out of", len(colors), "p-values are significant at <=", args.p_cutoff)
    ax.scatter(xticks, data.loc['actual', :], s=80, c=colors, zorder=99)

    # Prettify tick labels and axis labels
    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels, rotation="vertical", weight="bold")
    [t.set_weight('bold') for t in ax.get_yticklabels()]
    ax.tick_params(left='on', top='off', right='off', bottom='on')
    ax.set_ylabel(ax.get_ylabel(), fontsize='x-large', weight='bold')

    # Change violin colors
    violincolor = 'darkblue'
    for v in violins:
        for body in v['bodies']:  # Violin body
            body.set_facecolor(violincolor)
            body.set_linewidth(0)
            body.set_alpha(0.5)
        v['cbars'].set_color(violincolor)
        v['cbars'].set_linewidth(0.25)
        v['cmaxes'].set_color(violincolor)
        v['cmins'].set_visible(False)

    # Axis limts
    pad = 2
    ax.set_xlim(min(xticks) - pad, max(xticks) + pad)
    ylim = ax.get_ylim()
    ylim = [-0.01, ylim[1]]
    ax.set_ylim(ylim)

    # Custom legend
    red_dots = ax.scatter([], [], color='red', label='Significant at p ≤ ' + str(args.p_cutoff))
    gray_dots = ax.scatter([], [], color='gray', label='Not significant')
    ax.legend(handles=[red_dots, gray_dots], scatterpoints=1, markerscale=3, frameon=False,
              prop={'size': 'medium', 'weight': 'bold'})


def make_otu_key(ids, taxonomy):
    key = dict()
    for id, lineage in zip(ids, taxonomy):
        # print(id,"\n\t",lineage)
        myname = "UNKNOWN"
        clades = [re.sub("^.__", string=l, repl="") for l in lineage]
        # If not assigned
        if clades[0] == 'Unassigned':
            key[id] = "Unassigned"
            continue
        # If have genus and species, use
        if (len(clades) > 1) and (clades[-1] != "") and (clades[-2] != ""):
            myname = clades[-2] + " " + clades[-1]
            key[id] = myname
            continue
        # Otherwise, take last level
        for i in range(1, len(clades) + 1):
            if clades[-i] == "": continue
            level = find_level(lineage[-i])
            if level == 'genus':
                myname = clades[-i] + " sp."
            else:
                myname = "Unnamed " + level + " " + clades[-i]
            break  # Break so take the last available one
        key[id] = myname
        # print("\t", myname)
    return key


def find_level(clade):
    if clade.startswith("k__"): return "kingdom"
    if clade.startswith("p__"): return "phylum"
    if clade.startswith("c__"): return "class"
    if clade.startswith("o__"): return "order"
    if clade.startswith("f__"): return "family"
    if clade.startswith("g__"): return "genus"
    if clade.startswith("s__"): return "species"
    return "unkown_level"


if __name__ == '__main__': main()
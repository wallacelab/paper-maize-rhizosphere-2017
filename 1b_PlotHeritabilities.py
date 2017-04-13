__author__ = 'jgwall'

import argparse
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

debug = False


def main():
    args = parse_args()
    print("Plotting heritabilities from",args.infile)

    # Load and sort data
    data=pd.read_csv(args.infile, sep='\t')
    order = np.argsort(data.loc['actual',:])[::-1]
    data = data.iloc[:,order]

    # Plot
    fig = plt.figure(figsize=(5 + .25 * len(data.columns), 5))
    grid = gridspec.GridSpec(nrows=100, ncols=100)
    ax = fig.add_subplot(grid[:80,:], title="Distributions of null heritabilities", xlabel="trait", ylabel="Heritability")

    # Violin plots of random permutations
    perms = data.index != "actual"
    xticks, xlabels = list(), list()
    for x in range(len(data.columns)):
        trait = data.columns[x]
        xticks.append(x)
        ax.violinplot(data.loc[perms, trait], positions=[x])
        xlabels.append(trait.replace("trait_", ""))

    # Add dots for actual heritabilities
    ax.scatter(xticks, data.loc['actual',:])

    # Prettify
    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels, rotation="vertical")

    # Save
    fig.savefig(args.outfile, dpi=100)



def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile")
    parser.add_argument("-o", "--outfile")
    parser.add_argument("--debug", default=False, action="store_true")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug

    return parser.parse_args()


if __name__ == '__main__': main()
__author__ = 'jgwall'

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

debug = False

# Order to plot bars in
plot_order={"location":3, "week":2, "week:location":1, "inbred":7, "location:inbred":6, "week:inbred":5,
            "week:location:inbred":4, "Residuals":8}
colorkey={"location":"lightsalmon",
          "week":"tomato",
          "week:location":"orangered",
          "inbred":"forestgreen",
          "location:inbred":"lightskyblue",
          "week:inbred":"dodgerblue",
          "week:location:inbred":'royalblue',
          "Residuals":'lightgray'}

def main():
    args = parse_args()
    print("Graphing sum of squares divisions in",args.infile)

    # Load data
    data=pd.read_table(args.infile)
    roworder = [plot_order[i] for i in data.index]
    data=data.iloc[np.argsort(roworder),:]

    if args.num_pcs:
        data=data.iloc[:,:args.num_pcs]

    # data.index = [prettify_terms(i) for i in data.index]

    # Make graphic
    fig = plt.figure(figsize=(15, 6))
    grid = gridspec.GridSpec(nrows=1, ncols=2, hspace=0.5, wspace=0.25)
    ax_raw = fig.add_subplot(grid[:,0])
    ax_fract = fig.add_subplot(grid[:,1])

    plot_bars(data, ax_raw, normalize=False)
    plot_bars(data, ax_fract, normalize=True)


    # Save image
    fig.savefig(args.outprefix + ".png", dpi=100)
    fig.savefig(args.outprefix + ".svg", dpi=600)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile")
    parser.add_argument("-o", "--outprefix")
    parser.add_argument("-n", "--num-pcs", type=int, help='Number of PCs to plot')
    parser.add_argument("--debug", default=False, action="store_true")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug

    return parser.parse_args()


def plot_bars(data, ax, normalize=False):
    if normalize:
        for col in data.columns:
            data[col] = data[col] / np.nansum(data[col])
    bottoms = np.cumsum(data, axis=0)

    # Make color key
    colors = [colorkey[i] for i in data.index]
    xvals = np.arange(len(data.columns))


    for r in reversed(range(1, len(data.index))):
        ax.bar(left = xvals, height = data.iloc[r, :], color = colors[r], bottom=bottoms.iloc[r-1,:], label=prettify(data.index[r]))
    ax.bar(left=xvals, height=data.iloc[0, :], color=colors[0], label=prettify(data.index[0]))

    legend = ax.legend(framealpha=0.9, fontsize="x-small")
    [t.set_weight('bold') for t in legend.get_texts()]

    # Axis labels
    ax.set_xlabel("Principal Coordinate #", weight="bold", fontsize='x-large')
    ax.set_ylabel("Proportional Variance Explained" if normalize else "Raw Variance Explained", weight="bold", fontsize='x-large')
    xticks = np.arange(1, len(data.columns)+1) - 0.5
    xlabels = [str(int(x)) for x in np.ceil(xticks)]
    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels, fontsize="small", weight="bold")
    for ytick in ax.get_yticklabels():
        ytick.set_fontsize('small')
        ytick.set_weight('bold')
    ax.tick_params(axis='both', left='on', top='off', right='off', bottom='off')
    ymax=np.max(np.sum(data))
    ax.set_ylim([0, ymax * 1.03])



def prettify(label):
    label = label.replace("week","Week")
    label = label.replace("inbred", "Inbred")
    label = label.replace("location", "Location")
    label = reversed(label.split(":"))
    label = " x ".join(label)
    return(label)

if __name__ == '__main__': main()
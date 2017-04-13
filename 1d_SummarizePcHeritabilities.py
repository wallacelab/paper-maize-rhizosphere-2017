__author__ = 'jgwall'

import argparse
import math
import matplotlib
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re

debug = False


def main():
    args = parse_args()
    if debug: args.infiles=args.infiles[:2]
    print("Summarizing heritabilities from",len(args.infiles),"input files")

    # Load data
    data = [load_herit(i) for i in args.infiles]
    compiled = dict()
    for d in data:
        for trait in d:
            if trait in compiled:
                print("WARNING! Trait",trait,"exists twice in the dataset! Only one will be retained")
            compiled[trait] = d[trait]
    print("\tLoadeed",len(compiled),"traits")
    

    # # Write text output
    # columns=['trait','herit','pval','perm_max','herit_minus_perm_max']
    # print("Writing text summary to",args.outfile,"\n")
    # OUT = open(args.outfile, 'w')
    # OUT.write("\t".join(columns) + "\n")
    # for trait in sorted(compiled.keys()):
    #     values = [str(compiled[trait][c]) for c in columns]
    #     OUT.write("\t".join(values) + "\n")
    # OUT.close()
    
    
    # Make graphical output
    locations = get_uniques(compiled, 'location')
    ages = get_uniques(compiled, 'age')
    pcs = get_uniques(compiled, 'pc')


    ncol = 5
    nrow = int(len(locations) * 2 / ncol)+1
    fig = plt.figure(figsize=(6 * ncol, 0.75 * len(ages) * nrow + 4))
    grid = gridspec.GridSpec(nrows=nrow, ncols=ncol, hspace=0.25, wspace=0.25)

    myrow=0
    mycol=0
    for loc in locations:
        ax_h2 = fig.add_subplot(grid[myrow, mycol], title=loc + " - H2", xlabel="PC", ylabel="Week")
        ax_pval = fig.add_subplot(grid[myrow+1, mycol], title=loc + " - Empirical p-value", xlabel="PC", ylabel="Week")

        matrix_h2 = make_data_matrix(compiled, location=loc, key="herit")
        matrix_pval = make_data_matrix(compiled, location=loc, key="pval")

        make_heatmap(ax_h2, matrix_h2)
        make_heatmap(ax_pval, matrix_pval, reverse=True, log_transform=True, cmap_name="Reds")

        # Prettify
        prettify(ax_h2)
        prettify(ax_pval)

        mycol +=1
        if mycol == ncol:
            mycol=0
            myrow += 2

    ## Save
    fig.savefig(args.outgraphic, dpi=100)



def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infiles", nargs="*")
    parser.add_argument("-o", "--outfile")
    parser.add_argument("-g", "--outgraphic")
    parser.add_argument("--debug", default=False, action="store_true")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug

    return parser.parse_args()

def load_herit(infile):
    # Load data and identify which is the actual H2 value
    data=pd.read_csv(infile, sep='\t')
    perms = np.array(data.index)
    isperm = perms != "actual"

    # Load results into a dictionary for easier handling
    result=dict()
    for trait in data.columns:
        if trait in result:
            print("\tWARNING! Trait",trait,"in",infile,"occurs more than once!")
        values = np.array(data[trait])
        result[trait] = dict()
        result[trait]["trait"] = trait
        result[trait]["herit"] = values[~isperm][0] # Index 0 to make it no longer an array
        result[trait]["perms"] = values[isperm]
        result[trait]["pval"] = np.sum((result[trait]["herit"] < result[trait]["perms"]) / len(result[trait]["perms"])) # Empirical p-value calculation
        result[trait]["perm_max"] = np.max(result[trait]["perms"]) # Empirical p-value calculation
        result[trait]["herit_minus_perm_max"] = result[trait]["herit"] - result[trait]["perm_max"]
        result[trait]['location'], result[trait]['age'], result[trait]['pc'] = parse_trait_name(trait)

    return result

def parse_trait_name(trait):
    #pattern='trait_(.+)_otu_table_(.+)_even[^_]+__AGE_([^_]+)___pc_(PC.+)' # For patterns of "trait_unweighted_unifrac_otu_table_Ithaca_even10000__AGE_4___pc_PC1"
    pattern='trait_(.+)\.week(.+)\.(.+)_(PC.+)'     #for patterns of trait_unweighted.week05.Ithaca_PC2
    matcher = re.search(string=trait, pattern=pattern)
    name =  matcher.group(1) + ":" + matcher.group(3)   # Name = type of distance matrix + location
    return name, int(matcher.group(2)), matcher.group(4)    # Return name, age, and PC number

def get_uniques(data, key):
    myset = {data[trait][key] for trait in data}
    return sorted(myset)

# Make a data matrix to plot, with weeks in rows and PCs in columns
def make_data_matrix(data, location, key):
    # Subset data to just this location
    subdata=dict()
    for trait in data:
        if data[trait]["location"] == location:
            subdata[trait] = data[trait]

    # Gather data for conversion to data frame
    tempdata = dict()
    for trait in subdata:
        mypc = subdata[trait]["pc"]
        myweek = subdata[trait]["age"]
        if mypc not in tempdata:
            tempdata[mypc] = dict()
        if myweek in tempdata[mypc]:    # Check that this is new data; should only be one per dataset
            print("WARNING! Data for location",location,"age",myweek,"and PC",mypc,"has multiple values! Trait:", trait)
            print("\t",trait,":",key,":",subdata[trait][key],"versus",tempdata[mypc][myweek])
        tempdata[mypc][myweek] = subdata[trait][key]

    # Make sure all weeks have data
    weeks = list(range(1,16)) + [20]
    pcs = ["PC1", "PC2", "PC3", "PC4","PC5"]
    for week in weeks:
        for pc in pcs:
            if pc not in tempdata:
                tempdata[pc] = dict()
            if week not in tempdata[pc]:
                tempdata[pc][week] = np.nan
    return pd.DataFrame(tempdata)

def make_heatmap(ax, matrix, reverse=False, min_value=None, log_transform=False, cmap_name = "Blues"):
    nrow, ncol = matrix.shape
    offset=0.5

    # Place text
    for r in range(nrow):
        for c in range(ncol):
            myval = matrix.iloc[r, c]
            myval = round(myval, 4)
            ax.text(x=c + offset, y=r + offset, s=str(myval), fontsize="medium", ha="center", verticalalignment='center',
                    bbox={"color": "white", "alpha": 0.5, "boxstyle": "round"}, zorder=99)
            # Fix nans so doesn't mess up plotting
            if np.isnan(myval):
                matrix.iloc[r, c] = 1 if log_transform else 0


    # Make color plot, including various transformations
    cm = plt.get_cmap(cmap_name)
    if reverse: cm=plt.get_cmap(cmap_name + "_r")
    normalize = None
    if min_value is not None:   # Code to make there be a color cutoff
        max = np.max(np.ravel(matrix))
        if max < min_value: max = min_value + 1
        normalize=matplotlib.colors.Normalize(vmin=min_value, vmax=max, clip=True)
    if log_transform:
        ax.pcolormesh(np.array(np.log(matrix + 0.00001)), cmap=cm, norm = normalize, zorder=20)    # Add a small amount to prevent log(0) errors
    else:
        ax.pcolormesh(np.array(matrix), cmap=cm, norm = normalize, zorder=20)



    #Kludge to fix Ithaca missing week 1
    addone = False if 1 in matrix.index else True

    # Prettify things (axis labels, etc)
    ax.invert_yaxis()  # Flip so goes from top to bottom
    ax.set_xticks(np.arange(ncol) + offset)
    ax.set_xticklabels(matrix.columns, fontsize="x-large")
    ax.set_yticks(np.arange(nrow) + offset)
    ax.set_yticklabels(matrix.index, fontsize="x-large")
    ax.tick_params(bottom='off', top='off', left='off', right='off')
    ax.set_title(label=ax.get_title(), fontsize="x-small")

# Function to make axes pretty
def prettify(ax):
    ax.set_xlabel(ax.get_xlabel(), fontsize="xx-large", weight="bold", zorder=100)
    ax.set_ylabel(ax.get_ylabel(), fontsize="xx-large", weight="bold", zorder=100)

    oldtitle = ax.get_title()
    components = re.search(pattern="(.+):(.+) - (.+)",string=oldtitle)
    type, place, stat = components.group(1), components.group(2), components.group(3)
    if stat == "H2": stat = "H$^2$"
    if stat == "Empirical p-value" : stat = "p-value"
    newtitle=type + " " + stat + '\n' + place
    ax.set_title(newtitle, fontsize="xx-large", weight="bold", zorder=100, verticalalignment='bottom')


if __name__ == '__main__': main()
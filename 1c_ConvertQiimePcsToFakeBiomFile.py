__author__ = 'jgwall'

import argparse
import math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pandas as pd
from os.path import commonprefix
import re

debug = False

matplotlib.rcParams.update({'font.size':10})

def main():
    args = parse_args()
    pcs = load_pcs(args.infiles, args.num_pcs, args.stem)
    output_fake_biom(pcs, args.outfile)
    if args.outgraphic:
        output_pheno_distributions(pcs, args.outgraphic)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infiles", nargs="*", help="Principal components files from QIIME")
    parser.add_argument("-o", "--outfile", help="Output phenotype file in TASSEL format")
    parser.add_argument("-g", "--outgraphic", help="Output graphic of histograms of each phenotype")
    parser.add_argument("-n", "--num-pcs", type=int, default=1, help="Number of principal components to take")
    parser.add_argument("-s", "--stem", default="", help="Stem to append to beginning of PC name")
    parser.add_argument("--debug", default=False, action="store_true")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug

    return parser.parse_args()


def load_pcs(infiles, num_pcs, stem):
    if debug: infiles=infiles[:5]   # Debug mode = take only first 5 input files
    print("Loading principal components from",len(infiles),"input files")
    data = [load_pc_file(f) for f in infiles]
    print("\tTaking the first",num_pcs,"components of each")
    data = [d[:,:num_pcs+1] for d in data]

    # Handle dataset names
    names = [i.replace('/', '_') for i in infiles]
    left_trim = commonprefix(names)
    #right_trim = commonprefix([n[::-1] for n in names])[::-1]    # [::-1] =  Unfortunate quick syntax to reverse a string
    right_trim = "_pc.txt"
    newnames = [re.sub(pattern="^" + left_trim, repl="", string=s) for s in names]
    newnames = [re.sub(pattern=right_trim + "$", repl="", string=s) for s in newnames]
    # for old, new in zip(names, newnames):
    #     print(new,"from",old)

    # Convert to individual DataFrames for easier unification
    colnames = ["PC" + str(p) for p in range(1, num_pcs+1)]
    data = [pd.DataFrame(d[:,1:], index=d[:,0], dtype=float, columns=colnames) for d in data]  # Slice first columns out to be the index
    for mydata in data:
        mydata.columns = [stem + "_" + c for c in mydata.columns] # Add dataset name as part of column name
    # print(data[0].head())

    # Unify all DataFrames into one
    master = pd.concat(data, axis=1)
    # print(master.head())

    # Return a concatenated dataframe
    return master

# Helper function to load a single PC file from QIIME
def load_pc_file(infile):
    IN=open(infile, "r")
    for line in IN:    # Advance through lines until get to right spot
        if line.startswith("Site"): break
    site, nrow, ncol = line.strip().split('\t')
    matrix = list()
    for line in IN:
        if line == "\n": break  # When hit white space, break out
        data = line.strip().split('\t')
        matrix.append(data)
    matrix = np.array(matrix)
    IN.close()
    return(matrix)



def output_fake_biom(pcs, outfile):
    print("Outputting",len(pcs.columns), "traits in a pseudo-BIOM format to",outfile)
    pcs.index.name="Taxon"
    OUT = open(outfile,"w")
    OUT.write("# NOT constructed from a BIOM file\n")
    OUT.write("#OTU ID")
    pcs.transpose().to_csv(OUT, sep='\t', header=True, index=True, na_rep="NA")
    OUT.close()

def output_pheno_distributions(pcs, outgraphic):
    print("Outputting phenotype distributions to",outgraphic)
    # Determine plotting dimensions
    nplots = len(pcs.columns)
    nrow = ncol = math.ceil(math.sqrt(nplots))
    if nrow * (ncol-1) > nplots: ncol-=1    # Remove a row if unneeded

    # Set up figure
    fig = plt.figure(figsize=(ncol*4, nrow*3))
    grid = gridspec.GridSpec(nrows=nrow, ncols=ncol, hspace=0.5, wspace=0.5)

    # Plot each set of PCs as a histogram
    i=0
    for row in range(nrow):
        for col in range(ncol):
            if i >= nplots: break   # Since will have empty canvases at the end, break out once hit the end
            mydata = np.array(pcs.iloc[:,i])
            mydata = mydata[np.isfinite(mydata)]
            mytitle = pcs.columns[i]
            ax = fig.add_subplot(grid[row, col], xlabel="Distribution of values", ylabel="Count")
            ax.hist(mydata, bins=25)
            ax.set_title(mytitle, fontsize="xx-small", weight="bold")
            i+=1


    fig.savefig(outgraphic, dpi=100)


if __name__ == '__main__': main()
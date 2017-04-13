#! /bin/bash

# Make publication-level graphics out of prior data

rawdir=0_Rawdata/
vardir=2_Variances
plotdir=3_PublicationGraphics

if [ ! -e $plotdir ]; then mkdir $plotdir; fi

# # Variance components of PCs
for set in weighted unweighted; do
  python3 3a_PlotSumSquaresPretty.py -i $vardir/2a_${set}_pcs.sum_squares.txt -o $plotdir/3a_sum_squares.$set -n 10
#   break
done

# # PC heritability grids
# srcdir=1_BroadHerit_Update2_pubs/1f_pc_heirtability
# for set in weighted unweighted; do
#   python3 3b_SummarizePcHeritabilities_pretty.py -i $srcdir/1f_${set}*heritabilities.txt -o $plotdir/3b_${set}_pc_heritabilities.summary.txt -g $plotdir/3b_${set}_pc_heritabilities.summary --exclude Columbia Urbana #--debug
# #   break
# done

# Heritability of OTUs
# python3 3c_PlotOtuHeritabilities_two_column.py  -i 1_BroadHerit_Update2_pubs/1a_otu_heritabilities.combined.txt -o $plotdir/3c_broad_heritabilities -p 0.001 \
#   -b 0_RawdataUpdate/80_percShared_table_and_mapping/otu_table_80percShared_relAbundance.biom --top-n 200
#! /bin/bash

rawdir=0_RawData
vardir=2_Variances

if [ ! -e $vardir ] ;then mkdir $vardir; fi

keyfile="0_RawData/mapping_with_abundance.txt"

# Cut out PCs
for set in weighted unweighted; do
  # Slice out first 20 PCs
  tail -n +10 $rawdir/${set}_unifrac_pc.txt | cut -f1-21 > $vardir/2_${set}_pcs.txt    
  
  # Run analysis
  Rscript 2_QuantifyVarianceComponents.r -i $vardir/2_${set}_pcs.txt -o $vardir/2a_${set}_pcs.sum_squares -k $keyfile
done  
 
#! /usr/bin/env bash

# Calculate heritabilities for individual OTUs and also for community metrics (esp. beta diversity)

TASSEL5="perl /home/jason/Software/TASSEL/tassel-5-standalone/run_pipeline.pl -Xms10g -Xmx40g"

rawdir=0_RawData
parsedir=0a_ParsedData
broaddir=1_BroadHerit

if [ ! -e $parsedir ]; then mkdir $parsedir; fi
if [ ! -e $broaddir ]; then mkdir $broaddir; fi

keyfile=$rawdir/mapping_with_abundance.txt
biom=$rawdir/otu_table_80percShared_relAbundance.biom
biom_text=${biom/.biom/.txt}


##################
# Raw data manipulations
##################

# # Get summaries to look at stats
# biom summarize-table -i $biom -o $parsedir/0a_biom_summary.samples.txt
# biom summarize-table -i $biom -o $parsedir/0a_biom_summary.observations.txt --observations

# Convert biom to text
# biom convert -i $biom -o $biom_text --to-tsv


###########
# Broad-sense heritability
###########

# n_runs=10
# perms_per_run=500
# for i in $(seq $n_runs); do
#   Rscript 1a_BroadSenseHeritability.r -i $biom_text -k $keyfile --log-transform --num-cores 47 \
#     --heritfile $broaddir/1a_otu_heritabilities.$i.txt --random-perms $perms_per_run --covariates "Seq_Count_80Perc + (1|AGE) + (1|ENV)" --seed $i --rescale Seq_Count_80Perc
# done
# Rscript 1b_RecombineHeritabilities.r -i $broaddir/1a_otu_heritabilities.*.txt -o $broaddir/1a_otu_heritabilities.combined.txt
# python3 1b_PlotHeritabilities.py -i $broaddir/1a_otu_heritabilities.combined.txt -o $broaddir/1b_otu_heritabilities.combined.png


##################
# Principal Components heritability by week
##################

# splitdir=$broaddir/1c_split_distances
# pcdir=$broaddir/1d_split_pcs
# biomdir=$broaddir/1e_pc_biom_files
# outdir=$broaddir/1f_pc_heirtability
# 
# if [ ! -e $splitdir ]; then mkdir $splitdir; fi
# if [ ! -e $biomdir ]; then mkdir $biomdir; fi
# if [ ! -e $outdir ]; then mkdir $outdir; fi

# # Split distance matrices by weeks
# for set in weighted unweighted; do
#   distances=$rawdir/../dm_pc_files_update2/unweighted_unifrac_dm.txt
#   Rscript 1c_SplitDistanceMatrices.r -i $distances -o $splitdir/1c_distances.$set -k $keyfile --splits AGE_fixed ENV 
# done
# 
# # Get PCs with QIIME
# principal_coordinates.py -i $splitdir -o $pcdir


# # # Principal components by week: 
# for pcs in $pcdir/*.txt; do
#   stem=${pcs/*distances./}  # Remove front part
#   stem=${stem/.txt/}    # Remove back part
#   
#   echo -e "#####\n$stem\n#####"
#   
#   biom=$biomdir/1e_$stem.biom.txt
#   graphic=$biomdir/1e_$stem.pc_dist.png
#   herits=$outdir/1f_${stem}.heritabilities.txt
#   heritgraph=$outdir/1f_${stem}.heritabilities.png
#   python3 1c_ConvertQiimePcsToFakeBiomFile.py -i $pcs -o $biom --num-pcs 5 --outgraphic $graphic --stem $stem
#   Rscript 1a_BroadSenseHeritability.r -i $biom -o $outdir/1e_${stem}.blups.txt -k $keyfile --num-cores 8 \
#     --heritfile $herits --random-perms 1000 --covariates ""  # No covariates because age and environment are fixed for these samples
#   python3 1b_PlotHeritabilities.py -i $herits -o $heritgraph
# #   break
# done

# python3 1d_SummarizePcHeritabilities.py -i $outdir/*heritabilities.txt -o $broaddir/1d_pc_heritabilities.summary.txt -g $broaddir/1d_pc_heritabilities.pretty.png #--debug


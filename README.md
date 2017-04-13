# Wallace scripts for maize rhizosphere analysis
This repository contains all bioinformatic scripts used in analyze the heritability of the maize rhizosphere in Walters et al. 2017 [citation pending].

Script 1_CalculateHeritabilities.sh calculates broad-sense heritabilities of individual OTUs and also principal coordinates (split by week and location).

Script 2_VarianceComponents.sh determines how different potential factors feed into the principal coordinates across the entire dataset.

Script 3_PrettifyGraphics.sh takes the output of the above two scripts and reformats it into publication-ready figures.

Raw data should be released with 0_RawData.zip. The MD5 sums of each file are in rawdata_md5.txt

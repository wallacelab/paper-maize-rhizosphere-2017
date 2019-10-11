# Wallace scripts for maize rhizosphere analysis
This repository contains all bioinformatic scripts used in analyze the heritability of the maize rhizosphere in [Walters et al. 2018, Large-scale replicated field study of maize rhizosphere identifies heritable microbes](https://www.pnas.org/content/115/28/7368).

Script 1_CalculateHeritabilities.sh calculates broad-sense heritabilities of individual OTUs and also principal coordinates (split by week and location).

Script 2_VarianceComponents.sh determines how different potential factors feed into the principal coordinates across the entire dataset.

Script 3_PrettifyGraphics.sh takes the output of the above two scripts and reformats it into publication-ready figures.

Raw data should be released with 0_RawData.zip. The MD5 sums of each file are in rawdata_md5.txt (Unfortunately, this file is ~300 MB, far too large for Github, so it is available [here](https://outlookuga-my.sharepoint.com/:u:/g/personal/jgwall_uga_edu/ETf8IAYpjUFLmJQ1bf_h3pMBseAGvGDq4X3zb4scMSl5Vg?e=MRj6T2). (There should be a Download link near the top to grab the entire archive instead of navigating through the zipped file structure.)

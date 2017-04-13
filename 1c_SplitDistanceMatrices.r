#! /usr/bin/Rscript

# Convert OTU counts (transformed or otherwise) into BLUPs

library(argparse)
parser=ArgumentParser()
options(stringsAsFactors=F)
parser$add_argument("-i", "--infile", help="Input file of traits from a text OTU table")
parser$add_argument("-o", "--outprefix", help="Prefix for output files")
parser$add_argument("-k", "--keyfile", help="QIIME-formatted key file of sample metadata")
parser$add_argument("-s", "--splits", nargs="*", help="Which keyfile columns to split on")
args=parser$parse_args()
#setwd('/home/jgwall/Projects/LeyRhizosphereHeritability/0_RawdataUpdate/dm_pc_files_update2/')
# args=parser$parse_args(c("-i","unweighted_unifrac_dm.txt","-o",'99_tmp','-k','../80_percShared_table_and_mapping/mapping_with_abundance.txt', "-s", "AGE_fixed", "ENV"))

cat("Subsetting distance matrix in",args$infile,"by",args$splits,'\n')
key=read.delim(args$keyfile, check.names=F, row.names=1)
distances=as.matrix(read.delim(args$infile, row.names=1, header=T))

# Make sure row and column names match
if(!identical(rownames(distances), colnames(distances))){
  stop("ERROR! Row and column names do not match!!\n")
}

# Match key
keymatch = match(rownames(distances), rownames(key))
if(sum(is.na(keymatch))>0){
  warn("WARNING! Some samples not found in sample key\n")
}else{
  cat("\tAll samples found in sample key\n")
}
newkey = key[keymatch,]
if(!identical(rownames(newkey), rownames(distances))){
  stop("ERROR! Match between key and distance matrix failed!\n")
}

# Split
splitdata=newkey[args$splits]
splitkey = apply(splitdata, MARGIN=1, FUN=paste, collapse='.')
splits=unique(splitkey)
cat("Writing",length(splits),"subset distance matrices with output prefix",args$outprefix,"\n")
for(s in splits){
  subdist = distances[splitkey==s, splitkey==s]
  outfile=paste(args$outprefix,s,"txt", sep='.')
  write.table(subdist, file=outfile, sep='\t', quote=F, col.names=NA)
}
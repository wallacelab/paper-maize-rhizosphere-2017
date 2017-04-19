#! /usr/bin/Rscript

#Determine which OTUs to keep based on empirical p-values
library(argparse)
parser=ArgumentParser()
parser$add_argument("-i", "--infiles", nargs="*", help="Input files with actual and permuted heritabilities")
parser$add_argument("-o", "--outfile", help="Output file of OTUs to keep")
args=parser$parse_args()
#setwd('/home/jgwall/Projects/LeyRhizosphereHeritability/1_BroadHerit_Update2_pubs/')
# args=parser$parse_args(c('-i', '1a_otu_heritabilities.1.txt', '1a_otu_heritabilities.2.txt', '1a_otu_heritabilities.3.txt', '-o','99_tmp.txt') )


# Load data
cat("Loading data from",length(args$infiles),"input files\n")
results=lapply(args$infiles, read.delim)

#Split out actual values
actuals=lapply(results, function(x){
	return(subset(x, rownames(x)=='actual'))
  })
actuals=do.call(what=rbind, args=actuals)

# Make sure all actual values are the same
for(i in 1:nrow(actuals)){
  mismatch = actuals[i,] != actuals[1,]
  if(sum(mismatch)>0){
	  stop(sum(mismatch)," values of actual heritabilities differ at index ",i," compared to index 1\n")
	}
}
cat("\tAll actual values match\n")

# Combine permutations; also rename row names on my terms (so doesn't do automatically)
perms=lapply(1:length(results), function(i){
	myperms = subset(results[[i]], rownames(results[[i]])!='actual')
	rownames(myperms) = paste(rownames(myperms), i, sep='_')
	return(myperms)
  })
perms=do.call(what=rbind, args=perms)
cat("Combined data has",nrow(perms),"permutations\n")

# Write out
output = rbind(actuals[1,], perms)
cat("Writing combined data to",args$outfile,"\n")
write.table(output, file=args$outfile, sep='\t', quote=F, row.names=T, col.names=T)
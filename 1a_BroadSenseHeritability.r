#! /usr/bin/Rscript

# Convert OTU counts (transformed or otherwise) into BLUPs

library(argparse)
library(lme4)
library(parallel)
options(stringsAsFactors=F)
parser=ArgumentParser()
parser$add_argument("-i", "--infile", help="Input file of traits from a text OTU table")
parser$add_argument("-o", "--outfile", help="Output file of traits as a matrix")
parser$add_argument("-k", "--keyfile", help="QIIME-formatted key file of sample metadata")
parser$add_argument("-l", "--log-transform", default=FALSE, action="store_true", help="Whether to log-transform values")
parser$add_argument("-n", "--num-cores", default=8, type="integer", help="Number of parallel cores to run")
parser$add_argument("-s", "--subset", help="File of otus to keep; others will be ignored")
parser$add_argument("--seed", type="integer", default=1, help="Random seed for permutations")
parser$add_argument("-w", "--write-transformed", help="Write transformed values to this file")
parser$add_argument("-r", "--random-perms", default=0, type="integer", help="Number of randomly scrambled datasets to run")
parser$add_argument("-c", "--covariates", default="", help="Covariates (in model form, eg, 'AGE + (1|ENV)') for the model")
parser$add_argument("--heritfile", help="Output file for heritability (includes random permutation heritabilities if specified")
parser$add_argument("--rescale", help="Rescale this argument to a 0 to 10 range")
parser$add_argument("--debug", default=FALSE, action="store_true", help="Debug mode")
args=parser$parse_args()
#setwd('/home/jgwall/Projects/LeyRhizosphereHeritability/1_BroadHerit_Update2_pubs/')
# args=parser$parse_args(c("-i","../0a_ParsedData_Update2_pubs/0d_biom_top100_otus.txt","-o",'99_tmp.txt','-k','../0_RawdataUpdate/80_percShared_table_and_mapping/mapping_with_abundance.txt',
#   "--heritfile","99_herit.txt","-l","-r","10", "-c", "Seq_Count_80Perc", "--subset", "1b_otus_to_keep.txt"))

cat('Collapsing traits in',args$infile,"to BLUPs and calculating heritability\n")
data=as.matrix(read.delim(args$infile, check.names=F, row.names=1, skip=1))
key=read.delim(args$keyfile, check.names=F, row.names=1)
key$INBRED_nested = droplevels(interaction(key$ENV, key$AGE, key$INBREDS))

if(!is.null(args$rescale)){	# To keep variables in similar ranges
  cat("Rescaling",args$rescale,"to a 0 to 10 range\n")
  library(plotrix)
  key[[args$rescale]] = rescale(key[[args$rescale]], newrange=c(0,10))
}

model = paste(args$covariates, "(1|INBRED_nested)", sep=" + ")
cat("Model is: ",model,"\n")
target_column="INBRED_nested"	# which column making BLUPs on

# Subset out if specified
if(!is.null(args$subset)){
  cat("Subsetting OTUs to those specified in",args$subset,"\n")
  tosub = scan(args$subset, what=character())
  cat("\tLoaded",length(tosub),"OTUs to subset data to; original data has",nrow(data),"OTUs in it\n")
  data=data[rownames(data) %in% tosub,]
  cat("\tResulting data frame has",nrow(data),"OTUs to analyze\n")
}

# Log-transform OTUs if specified
if(args$log_transform){
  cat("Log-transforming OTUs (assuming fractional data, and this empirically makes them more normally distributed)\n")
  data=apply(data, MARGIN=2, FUN=function(x){x[x==0] = min(x[x!=0]/10); return(x)})	# Set anything equal to 0 in a sample to 10x less than the smallest non-zero element for that sample
  data=log(data)
  if(!is.null(args$write_transformed)){
	write("# Transformed from biom file", file=args$write_transformed)	# to maintain format
	write(c("#OTU ID",colnames(data)), file=args$write_transformed,  sep='\t', append=T, ncolumns=ncol(data)+1)	# Again, formatting
	write.table(data, file=args$write_transformed, sep='\t', quote=F, row.names=T, col.names=F, append=T)
  }
}

# Match data to key
mydata=t(data)
key = subset(key, rownames(key) %in% rownames(mydata))
mymatch = match(rownames(key), rownames(mydata))
mydata = mydata[mymatch,, drop=FALSE]

# Check debug mode
if(args$debug){mydata = mydata[,1:10]}

# Build data frame with the appropriate setup
colnames(mydata) = paste("trait",colnames(mydata),sep="_") # So formula won't cause problems with any numerically named traits
modeldata = cbind(key, mydata)


#Helper function to get BLUPs.
run.lmer = function(yname, model, modeldata, target_column){
	myformula=paste(yname, model, sep="~")
	subdata = subset(modeldata, is.finite(modeldata[,yname]))
	mymodel = lmer(myformula, data=subdata)
	
	# Get BLUPs
	values = ranef(mymodel)[[target_column]] # Value for each line
	blups=data.frame(row.names=row.names(values), vals=as.numeric(values[,1]))
	names(blups)[1]=yname
	
	# Get (broad-sense) heritability
	variances = as.data.frame(summary(mymodel)$varcor)
	heritability = variances$sdcor[variances$grp==target_column]/ sum(variances$sdcor)
	
	# Return data
	return(list(blups=blups, H2=heritability))
}

# Helper function to merge blups
merge_all = function(x,y){
	merged = merge(x,y,all=TRUE, by=0)
	rownames(merged) = merged$Row.names
	merged$Row.names=NULL
	return(merged)
}

# Run models
traits=colnames(mydata)
cat("Data has dimensions",dim(modeldata),"\n")
results = mclapply(traits, run.lmer, model=model, modeldata=modeldata, target_column=target_column, mc.cores=args$num_cores)

# Write BLUPs in TASSEL format
if(!is.null(args$outfile)){
  cat("Writing out BLUPs to", args$outfile,"\n")
  blups = Reduce(merge_all, lapply(results, function(x){x$blups}))
  write("<Phenotype>", file=args$outfile)
  write(c("taxa", rep("data", ncol(blups))), file=args$outfile, append=T, ncol=ncol(blups)+1, sep='\t')
  blups = cbind(rownames(blups), blups)
  names(blups)[1]="Taxon"
  write.table(blups, file=args$outfile, sep='\t', quote=F, row.names=F, col.names=T, append=T)
}

# Get heritabilities
herit = sapply(results, function(x){x$H2})
names(herit) = traits
herit=as.data.frame(t(herit))
rownames(herit)[1]="actual"

# Do random permutations if requested
if(args$random_perms >0){
  cat("Performing",args$random_perms,"random permutations for heritability analysis\n")
  set.seed(args$seed)	# For reproducibility
  
  perms = lapply(1:args$random_perms, function(x){
	  scramble=sample(1:nrow(mydata))
	})

  perm_h2 = lapply(perms, function(scramble){
	  permdata = cbind(key, mydata[scramble,], make.row.names=F)
	  myres = mclapply(traits, run.lmer, model=model, modeldata=permdata, target_column=target_column, mc.cores=args$num_cores, mc.preschedule=F )
	  myherit = sapply(myres, function(x){x$H2})
	  names(myherit) = traits
	  return(myherit)
	})
  
  perm_h2 = do.call(rbind, perm_h2)
  final_h2 = rbind(herit, perm_h2)
  rownames(final_h2) = c("actual", paste("perm",1:args$random_perms, sep=''))
}else{
  final_h2 = herit
}

# Write out heritability results
if(!is.null(args$heritfile)){
  cat("Writing heritability to",args$heritfile,"\n")
  write.table(final_h2, file=args$heritfile, sep='\t', quote=F, row.names=T, col.names=T)
}
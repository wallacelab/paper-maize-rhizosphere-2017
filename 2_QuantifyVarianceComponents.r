# Quick script to quantify the variance components aspects of all the data

# Load data 
#Combine OTU data into a single data frame
library(argparse)
library(car)
library(parallel)
options(contrasts = c("contr.sum","contr.poly")) # For type 3 SS
parser=ArgumentParser()
parser$add_argument("-k", "--keyfile", help="QIIME-formatted key file")
parser$add_argument("-i", "--infile", help="Principal component file")
parser$add_argument("-o", "--outprefix", help="Output file prefix")
parser$add_argument("-n", "--numcores", default=7, help="Number of parallel cores to use")
args=parser$parse_args()
#setwd('/home/jgwall/Desktop/jason/Projects/LeyRhizosphereHeritability/2_Variances/')
# args=parser$parse_args(c('-k','../0_RawdataUpdate/maize_mapping_with_run_data.txt','-i','2_weighted_pcs.txt','-o','99_tmp'))

# Load data
pcs = read.delim(args$infile, row.names=1, header=F)
key=read.delim(args$keyfile, stringsAsFactors=F)
names(key)[1]="Sample"
names(pcs)=paste("pc", 1:ncol(pcs),sep='')

# Match up to key
keymatch = match(rownames(pcs), key$Sample)
if(any(rownames(pcs) != key$Sample[keymatch], na.rm=T)){
  stop("Mismatch occured between keys and PCs\n")
}
data = data.frame(sample=key$Sample[keymatch], week=key$AGE_fixed[keymatch], location=key$ENV[keymatch], inbred=key$INBREDS[keymatch], run=key$RUN_PREFIX[keymatch], pcs)
data=subset(data, !is.na(data$week))	# remove ones with NAs

# Make models with various levels of nesting
full_model = mclapply(names(pcs), function(x){
	myformula = paste(x, "~ week*location*inbred", sep="")
# 	myformula = paste(x, "~ week+location+inbred", sep="")
	mymodel = lm(myformula, data=data)
}, mc.cores=args$numcores)

anovas = lapply(full_model, function(x){
  tmp=as.data.frame(anova(x))
  tmp = data.frame(SS=tmp$'Sum Sq', row.names=rownames(tmp))
  return(tmp)
})
ss = do.call(cbind, anovas)
names(ss) = names(pcs)

# Normalize
normalized = ss
for(i in 1:ncol(normalized)){
  normalized[,i] = normalized[,i] / sum(normalized[,i])
}

# Plot
set.seed(1)
colors = rainbow(nrow(ss)-1)
colors = c(sample(colors), "gray")
plotme=function(x, ...){
  barplot(as.matrix(x), legend.text=T, col=colors, las=2, args.legend=list(cex=0.6), ...)
}

png(paste(args$outprefix,"png", sep='.'), width=1000, height=500)
  par(mfrow=c(1,2), cex=1.4)
  plotme(ss, main="Raw SS", ylab = "Sum of Squares (raw)")
  plotme(normalized, main="Proportion of Total SS", ylab = "Fraction of Sum of Squares")
dev.off()

# Write out text
write.table(ss, file=paste(args$outprefix,"txt", sep='.'), sep='\t', quote=F, row.names=T, col.names=T)

#PURPOSE: Run WGCNA on Y introgression expression data

args<-commandArgs(TRUE)
if(length(args) != 1){
        stop("Missing command line arguments.\nArgument 1: dataset (options: all, exp1, exp2, exp1_LZ, exp1_RS, exp2_LZ, exp2_RS)")
}

dataset<-args[1] #all #exp1 #exp2
library(WGCNA)
options(stringsAsFactors = FALSE)
enableWGCNAThreads(nThreads=8)

print("Reading in expression data...")
#expData<-t(read.table("logTransformed_expData.txt", header=TRUE)) #t() so genes are columns and samples are rows
expData<-t(read.table("OLD_logTransformed_expData.noExpLevelFilter.txt", header=TRUE)) #t() so genes are columns and samples are rows
print(dim(expData))
print("Checking to make sure all genes good...")
gsg = goodSamplesGenes(expData, verbose = 3)
if(!(gsg$allOK)){
	# Optionally, print the gene and sample names that were removed:
	if(sum(!gsg$goodGenes)>0){
		printFlush(paste("Removing genes:", paste(names(expData)[!gsg$goodGenes], collapse = ", ")));
	}
	if(sum(!gsg$goodSamples)>0){
		printFlush(paste("Removing samples:", paste(rownames(expData)[!gsg$goodSamples], collapse = ", ")));
	}
	# Remove the offending genes and samples from the data:
	expData = expData[gsg$goodSamples, gsg$goodGenes]
}
print("All samples and genes good? Should be true:")
gsg = goodSamplesGenes(expData, verbose = 3)
print(gsg$allOK)

#NOTE (29 March 2021): Some weird relationships but no "outliers" - going to try running for all samples for now but might try removing weird samples also; should cross-check to see if they are the ones with likely contamination from other cell types
#Would remove samples here if I was going to
if(dataset=="all"){
	final_expData<-expData
} else if(dataset=="exp1"){
	print("Looking at experiment 1 data only...")
	final_expData<-expData[grepl("CP|WL", rownames(expData)),]
} else if(dataset=="exp2"){
	print("Looking at experiment 2 data only...")
	final_expData<-expData[grepl("LLPP|LPLY|PPLL|PLPY", rownames(expData)),]
} else if(dataset=="exp1_LZ"){
	print("Looking at experiment 1 LZ data only...")
        final_expData<-expData[grepl("CP.*LZ|WL.*LZ", rownames(expData)),]
} else if(dataset=="exp1_RS"){
        print("Looking at experiment 1 RS data only...")
        final_expData<-expData[grepl("CP.*RS|WL.*RS", rownames(expData)),]
} else if(dataset=="exp2_LZ"){
        print("Looking at experiment 2 LZ data only...")
        final_expData<-expData[grepl("LLPP.*LZ|LPLY.*LZ|PPLL.*LZ|PLPY.*LZ", rownames(expData)),]
} else if(dataset=="exp2_RS"){
        print("Looking at experiment 2 RS data only...")
        final_expData<-expData[grepl("LLPP.*RS|LPLY.*RS|PPLL.*RS|PLPY.*RS", rownames(expData)),]
} else{
	stop("Invalid dataset: must be 'all', 'exp1', or 'exp2'")
}

print("Dimensions of dataset:")
print(dim(final_expData))

print("Generating sample cluster to check for outliers:")
sampleTree = hclust(dist(final_expData), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
pdf(file = paste("sampleClustering",dataset,"pdf", sep="."), width = 12, height = 9)
#sizeGrWindow(12,9)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
dev.off()

print("Choosing soft-thresholding power...")
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=30, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(final_expData, powerVector = powers, verbose = 5)
print("Print sft power estimate:")
print(sft)
# Plot the results:
pdf(paste("softThresholding_power_curves",dataset,"pdf", sep="."))
#sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n", ylim=c(0,1),
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

#NOTE (29 March 2021): None of them get to 90%, seems to level off a little over 80% at power=12

save(final_expData, file=paste("03_WGCNA_dataChecks",dataset,"RData",sep="."))

print("Done with 03_WGCNA_dataChecks.r") 

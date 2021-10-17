#PURPOSE: Run automatic clustering and network analysis with WGCNA

args<-commandArgs(TRUE)
if(length(args) != 2){
        stop("Missing command line arguments.\nArgument 1: dataset (options: all, exp1, exp2, exp1_LZ, exp1_RS, exp2_LZ, exp1_RS)\nArgument 2: signed or unsigned? (options: TRUE, FALSE)")
}

dataset<-args[1] #all #exp1 #exp2 #exp1_LZ #exp1_RS #exp2_LZ #exp2_RS
signed<-args[2] #FALSE #Perform signed network analysis? SIGN RECOMMENDED

print(paste("Performing WGCNA analysis for",dataset,"; signed =",signed, sep="."))

library(WGCNA)

options(stringsAsFactors = FALSE)
enableWGCNAThreads(nThreads=8)

lnames<-load(file=paste("03_WGCNA_dataChecks",dataset,"RData",sep=".")) #Load data from previous script
#allTraits<-read.table("traits.txt",header=TRUE)
allTraits<-read.table("traits_noCross.txt",header=TRUE)
if(dataset=="exp1"){
	myTraits<-allTraits[grepl("CP|WL", rownames(allTraits)),]
	stPower<-12
} else if(dataset=="exp2"){
	myTraits<-allTraits[grepl("LLPP|LPLY|PPLL|PLPY", rownames(allTraits)),]
	stPower<-22
} else if(dataset=="all"){
	myTraits<-allTraits
} else if(dataset=="exp1_LZ"){
	myTraits<-allTraits[grepl("CP.*LZ|WL.*LZ", rownames(allTraits)),]
	stPower<-24
} else if(dataset=="exp1_RS"){
        myTraits<-allTraits[grepl("CP.*RS|WL.*RS", rownames(allTraits)),]
	stPower<-6
} else if(dataset=="exp2_LZ"){
        myTraits<-allTraits[grepl("LLPP.*LZ|LPLY.*LZ|PPLL.*LZ|PLPY.*LZ", rownames(allTraits)),]
	stPower<-10
} else if(dataset=="exp2_RS"){
        myTraits<-allTraits[grepl("LLPP.*RS|LPLY.*RS|PPLL.*RS|PLPY.*RS", rownames(allTraits)),]
	stPower<-3
} else{
	stop("Invalid dataset: must be 'all', 'exp1', 'exp1_LZ', 'exp1_RS', 'exp2', 'exp2_LZ', or 'exp2_RS'")
}

print("Performing network analysis...")
#power=stPower (where R^2 levels off, different for each dataset; see softThresholding_power_curves.<dataset>.pdf from 03_WGCNA_datachecks)
#maxBlockSize = number of genes; can make smaller if taking too much memory
if(signed){
	myNetwork<-blockwiseModules(final_expData, corType="pearson", power = stPower, networkType = "signed", minModuleSize = 30, reassignThreshold = 0, mergeCutHeight = 0.25, numericLabels = TRUE, pamRespectsDendro = FALSE, saveTOMs = TRUE, saveTOMFileBase = paste0(dataset,"Data"), verbose = 3, maxBlockSize = 29629)
} else{
	myNetwork<-blockwiseModules(final_expData, corType="pearson", power = stPower, networkType = "unsigned", minModuleSize = 30, reassignThreshold = 0, mergeCutHeight = 0.25, numericLabels = TRUE, pamRespectsDendro = FALSE, saveTOMs = TRUE, saveTOMFileBase = paste0(dataset,"Data"), verbose = 3, maxBlockSize = 29629)
}

if(signed){
	pdf(paste0(dataset,"Data.network.signed.pdf"))
} else{
	pdf(paste0(dataset,"Data.network.pdf"))
}
mergedColors<-labels2colors(myNetwork$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(myNetwork$dendrograms[[1]], mergedColors[myNetwork$blockGenes[[1]]], "Module colors", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()

#Get some info from network analysis
moduleLabels<-myNetwork$colors
moduleColors<-labels2colors(myNetwork$colors)
MEs<-myNetwork$MEs
geneTree<-myNetwork$dendrograms[[1]]
#save(MEs, moduleLabels, moduleColors, geneTree, file=paste("05-WGCNA",dataset,"RData",sep="."))
#Get relationship between label (number) and color
n_mods<-length(unique(moduleLabels))
labelOrder<-c(0:(n_mods-1))
colorOrder<-c()
for(i in labelOrder){
	colorOrder<-c(colorOrder, unique(moduleColors[which(moduleLabels==i)]))
}
label_to_color<-as.data.frame(cbind(lab=labelOrder, col=colorOrder))

#Get per-gene module membership (also known as kME or signed eigengene-based connectivity)
myKME<-signedKME(final_expData, MEs, outputColumnName = "kME", corFnc = "cor", corOptions = "use = 'p'")
#Get column names as module colors for kME
colname_nums<-sapply(colnames(myKME), function(x) gsub("kME","",x))
colname_cols<-sapply(colname_nums, function(x) label_to_color$col[which(label_to_color$lab==x)])
colnames(myKME)<-colname_cols

#Save for future scripts/analyses
save(MEs, moduleLabels, moduleColors, geneTree, label_to_color, myKME, file=paste("05-WGCNA",dataset,"RData",sep="."))

print("Testing for associations with trait data...")
# Define numbers of genes and samples
nGenes<-ncol(final_expData)
nSamples<-nrow(final_expData)
# Recalculate MEs with color labels
MEs0<-moduleEigengenes(final_expData, moduleColors)$eigengenes
MEs<-orderMEs(MEs0)
#Convert categorical traits to binary
myTraitsBin<-binarizeCategoricalColumns.forRegression(myTraits)
print(head(myTraitsBin))
#Calculate significant associations between trait and expression module
moduleTraitCor<-cor(MEs, myTraitsBin, use = "p")
moduleTraitPvalue_raw<-corPvalueStudent(moduleTraitCor, nSamples)
moduleTraitPvalue<-p.adjust(moduleTraitPvalue_raw, method="fdr")

if(signed){
	pdf(paste0(dataset,"Data.traitHeatmap.signed.pdf"), height=11, width=8.5)
} else{
	pdf(paste0(dataset,"Data.traitHeatmap.pdf"), height=11, width=8.5)
}
if( (dataset=="exp1_LZ") || (dataset=="exp1_RS") || (dataset=="exp2_LZ") || (dataset=="exp2_RS") ){
        myTraitsBin<-myTraitsBin[,which(colnames(myTraitsBin) != "cell_type.RS.vs.all")]
}
# Will display correlations and their p-values
textMatrix<-paste(signif(moduleTraitCor, 2), " (", signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix)<-dim(moduleTraitCor)
par(mar = c(9, 8.5, 3, 3))
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor, xLabels = names(myTraitsBin), yLabels = names(MEs), ySymbols = names(MEs), colorLabels = FALSE, colors = blueWhiteRed(50), textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1,1), main = paste("Module-trait relationships"))
dev.off()

#Look at cell type and XY mismatch only (not cross type)
if( (dataset=="exp1_LZ") || (dataset=="exp1_RS") || (dataset=="exp2_LZ") || (dataset=="exp2_RS") ){
	myTraitsBin_noCross<-myTraitsBin[,c("XYmismatch.yes.vs.all","XYdirection.musXdomY.vs.all","XYdirection.XYmatch.vs.all")]
} else{
	myTraitsBin_noCross<-myTraitsBin[,c("cell_type.RS.vs.all","XYmismatch.yes.vs.all","XYdirection.musXdomY.vs.all","XYdirection.XYmatch.vs.all")]
}
moduleTraitCor<-cor(MEs, myTraitsBin_noCross, use = "p")
moduleTraitPvalue_raw<-corPvalueStudent(moduleTraitCor, nSamples)
moduleTraitPvalue<-p.adjust(moduleTraitPvalue_raw, method="fdr")
if(signed){
	pdf(paste0(dataset,"Data.traitHeatmap.noCross.signed.pdf"), height=11, width=8.5)
} else{
	pdf(paste0(dataset,"Data.traitHeatmap.noCross.pdf"), height=11, width=8.5)
}
# Will display correlations and their p-values
textMatrix<-paste(signif(moduleTraitCor, 2), " (", signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix)<-dim(moduleTraitCor)
par(mar = c(9, 8.5, 3, 3))
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor, xLabels = names(myTraitsBin_noCross), yLabels = names(MEs), ySymbols = names(MEs), colorLabels = FALSE, colors = blueWhiteRed(50), textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1,1), main = paste("Module-trait relationships"))
dev.off()

#Look at cell type and geno (not cross type or XY mismatch as its own category)
if( (dataset=="exp1_LZ") || (dataset=="exp1_RS") ){
	myTraitsBin_ctAndGeno<-myTraitsBin[,c("autos.mus.vs.all","Xchr.mus.vs.all","Ychr.mus.vs.all")]
} else if( (dataset=="exp2_LZ") || (dataset=="exp2_RS") ){
	myTraitsBin_ctAndGeno<-myTraitsBin[,c("autos.PL.vs.all","Xchr.mus.vs.all","Ychr.mus.vs.all")]
}else if(dataset=="exp1"){
	myTraitsBin_ctAndGeno<-myTraitsBin[,c("cell_type.RS.vs.all","autos.mus.vs.all","Xchr.mus.vs.all","Ychr.mus.vs.all")]        
} else{
	myTraitsBin_ctAndGeno<-myTraitsBin[,c("cell_type.RS.vs.all","autos.PL.vs.all","Xchr.mus.vs.all","Ychr.mus.vs.all")]
}
moduleTraitCor<-cor(MEs, myTraitsBin_ctAndGeno, use = "p")
moduleTraitPvalue_raw<-corPvalueStudent(moduleTraitCor, nSamples)
moduleTraitPvalue<-p.adjust(moduleTraitPvalue_raw, method="fdr")
if(signed){
        pdf(paste0(dataset,"Data.traitHeatmap.cellTypeAndGeno.signed.pdf"), height=11, width=8.5)
} else{
        pdf(paste0(dataset,"Data.traitHeatmap.cellTypeAndGeno.pdf"), height=11, width=8.5)
}
# Will display correlations and their p-values
textMatrix<-paste(signif(moduleTraitCor, 2), " (", signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix)<-dim(moduleTraitCor)
par(mar = c(9, 8.5, 3, 3))
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor, xLabels = names(myTraitsBin_ctAndGeno), yLabels = names(MEs), ySymbols = names(MEs), colorLabels = FALSE, colors = blueWhiteRed(50), textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1,1), main = paste("Module-trait relationships"))
dev.off()

#Correlation between eigengene connectivity and phenotype (ex: RTM, sperm count, sperm morpho, sex ratio)
#SHOULD PROBABLY DO THIS FOR CELL TYPES SEPARATELY, OTHERWISE ALL PHENO DATA IS DUPLICATED
if( (dataset=="exp1_LZ") || (dataset=="exp1_RS") || (dataset=="exp2_LZ") || (dataset=="exp2_RS") ){
	phenoTraits<-read.table("traits_pheno.txt",header=TRUE)
	if(dataset=="exp1_LZ"){
	        myPhenoTraits<-phenoTraits[grepl("CP.*LZ|WL.*LZ", rownames(phenoTraits)),]
	} else if(dataset=="exp1_RS"){
		myPhenoTraits<-phenoTraits[grepl("CP.*RS|WL.*RS", rownames(phenoTraits)),]
	} else if(dataset=="exp2_LZ"){
	        myPhenoTraits<-phenoTraits[grepl("LLPP.*LZ|LPLY.*LZ|PPLL.*LZ|PLPY.*LZ", rownames(phenoTraits)),]
	} else if(dataset=="exp2_RS"){
	        myPhenoTraits<-phenoTraits[grepl("LLPP.*RS|LPLY.*RS|PPLL.*RS|PLPY.*RS", rownames(phenoTraits)),]
	} else{
	        stop("Invalid dataset: must be 'exp1_LZ', 'exp1_RS', 'exp2_LZ', or 'exp2_RS'")
	}
	#print(dim(MEs))
	#print(dim(myPhenoTraits))
#	nSamples<-nrow(myPhenoTraits)
#	colname_nums<-sapply(colnames(MEs), function(x) gsub("ME","",x))
#	colname_cols<-sapply(colname_nums, function(x) label_to_color$col[which(label_to_color$lab==x)])
#	colnames(MEs)<-colname_cols
	moduleTraitCor<-cor(MEs, myPhenoTraits, use = "p")
	moduleTraitPvalue_raw<-corPvalueStudent(moduleTraitCor, nSamples)
	moduleTraitPvalue<-p.adjust(moduleTraitPvalue_raw, method="fdr")
	if(signed){
	        pdf(paste0(dataset,"Data.traitHeatmap.reproPheno.signed.pdf"), height=11, width=8.5)
	} else{
	        pdf(paste0(dataset,"Data.traitHeatmap.reproPheno.pdf"), height=11, width=8.5)
	}
	# Will display correlations and their p-values
	textMatrix<-paste(signif(moduleTraitCor, 2), " (", signif(moduleTraitPvalue, 1), ")", sep = "")
	dim(textMatrix)<-dim(moduleTraitCor)
	par(mar = c(9, 8.5, 3, 3))
	# Display the correlation values within a heatmap plot
	labeledHeatmap(Matrix = moduleTraitCor, xLabels = names(myPhenoTraits), yLabels = names(MEs), ySymbols = names(MEs), colorLabels = FALSE, colors = blueWhiteRed(50), textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1,1), main = paste("Module-trait relationships"))
	dev.off()
}

print("Done with 05_WGCNA.r")

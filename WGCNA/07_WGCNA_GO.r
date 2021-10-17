#PURPOSE: Perform GO enrichment analysis for interesting modules

args<-commandArgs(TRUE)
if(length(args) != 1){
        stop("Missing command line arguments.\nArgument 1: dataset")
}

dataset<-args[1] #all #exp1 #exp2

print(paste0("Running GO enrichment analysis for ", dataset))

library(WGCNA)
library(anRichment)
library(biomaRt)

options(stringsAsFactors = FALSE)
enableWGCNAThreads(nThreads=8)

#Load data from previous script
lnames<-load(file=paste("03_WGCNA_dataChecks",dataset,"RData",sep="."))
lnames<-c(lnames, load(file=paste0(dataset,"Data-block.1.RData")))
lnames<-c(lnames, load(file=paste("05-WGCNA",dataset,"RData",sep=".")))
allTraits<-read.table("traits.txt",header=TRUE)
if(dataset=="exp1"){
        interesting_mods<-c("blue","turquoise","black","green","yellow","yellowgreen")
} else if(dataset=="exp2"){
        interesting_mods<-c("blue","turquoise","orangered4","yellow","plum1") #grey
} else if(dataset=="exp1_LZ"){
        interesting_mods<-c("green","black","brown")
} else if(dataset=="exp1_RS"){
        interesting_mods<-c("blue","magenta","darkorange","brown") #grey
} else if(dataset=="exp2_LZ"){
        interesting_mods<-c("turquoise","brown")
} else if(dataset=="exp2_RS"){
        interesting_mods<-c("greenyellow","pink","skyblue","white","orange")
} else{
        stop("Invalid dataset: must be 'exp1', 'exp1_LZ', 'exp1_RS', 'exp2', 'exp2_LZ', or 'exp2_RS'")
}

#Get Entrez IDs
mouse_mart<-useMart("ENSEMBL_MART_ENSEMBL", "mmusculus_gene_ensembl")
myIDs<-getBM(attributes=c("ensembl_gene_id","entrezgene_id","chromosome_name"), mart=mouse_mart)
entrez<-c()
for(gene in names(moduleLabels)){
	e<-myIDs$entrezgene_id[which(myIDs$ensembl_gene_id==gene)]
	entrez<-c(entrez, e)
}
table(is.finite(entrez))

#GO analysis
GOcollection<-buildGOcollection(organism = "mouse")
GOenrichment<-enrichmentAnalysis(classLabels = moduleColors, identifiers = entrez, refCollection = GOcollection, useBackground = "given", threshold = 1e-4, thresholdType = "Bonferroni", getOverlapEntrez = TRUE, getOverlapSymbols = TRUE, ignoreLabels = "grey")
collectGarbage()
write.csv(GOenrichment$enrichmentTable, file=paste("GOenrichment-enrichmentTable",dataset,"csv", sep="."), row.names=FALSE)

print("Done with 07_WGCNA_GO.r") 

#PURPOSE: Identify Overlap between modules; early and late as a sanity check but also interesting to see what is associated w/ background vs sex chr (ex: mus X and dom Y regardless of background vs mus X and dom Y but different if background is hybrid vs non-hybrid)
#TRY CONSENSUS ANALYSIS: See this tutorial https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/Consensus-RelateToFemMods.pdf

args<-commandArgs(TRUE)
if(length(args) != 1){
        stop("Missing command line arguments.\nArgument 1: Cell type (options: both, LZ, RS)")
}
ct<-args[1]

print(paste0("Comparing experiments for cell type: ", ct))

library(VennDiagram)
library(gridExtra)
library(ggplot2)
library(biomaRt)
library(WGCNA)

#load(paste("05-WGCNA",dataset,"RData", sep="."))
#load(paste0(dataset,"Data-block.1.RData"))

if(ct=="both"){
	load("05-WGCNA.exp1.RData")
	exp1_moduleLabels<-moduleLabels
	exp1_moduleColors<-moduleColors
	load("05-WGCNA.exp2.RData")
	exp2_moduleLabels<-moduleLabels
        exp2_moduleColors<-moduleColors
        exp1_interesting_mods<-c("blue","turquoise","black","green","yellow","yellowgreen")
	exp2_interesting_mods<-c("blue","turquoise","orangered4","yellow","plum1") #grey
	#Table where each column is an experiment and each row is an interesting comparison (i.e. modules with similar trait associations between the two experiments)
	myContrasts<-as.data.frame(cbind(exp1_mods=c("blue","turquoise"), exp2_mods=c("blue","turquoise")))
} else if(ct=="LZ"){
	load("05-WGCNA.exp1_LZ.RData")
        exp1_moduleLabels<-moduleLabels
        exp1_moduleColors<-moduleColors
        load("05-WGCNA.exp2_LZ.RData")
        exp2_moduleLabels<-moduleLabels
        exp2_moduleColors<-moduleColors
	exp1_interesting_mods<-c("green","black","brown")
	exp2_interesting_mods<-c("turquoise","brown")
	myContrasts<-as.data.frame(cbind(exp1_mods=c("green","black","brown"), exp2_mods=c("brown","turquoise","turquoise")))
} else if(ct=="RS"){
	load("05-WGCNA.exp1_RS.RData")
        exp1_moduleLabels<-moduleLabels
        exp1_moduleColors<-moduleColors
        load("05-WGCNA.exp2_RS.RData")
        exp2_moduleLabels<-moduleLabels
        exp2_moduleColors<-moduleColors
        exp1_interesting_mods<-c("blue","magenta","darkorange","brown") #grey
	exp2_interesting_mods<-c("greenyellow","pink","skyblue","white","orange")
	myContrasts<-as.data.frame(cbind(exp1_mods=c("blue","blue","magenta","magenta","brown","darkorange","grey"), exp2_mods=c("greenyellow","pink","greenyellow","pink","orange","white","skyblue")))
} else{
        stop("Invalid cell type: must be 'both', 'LZ', or 'RS'")
}

#Create venn diagrams of genes in modules with similar patterns/associated with same traits
for(i in 1:nrow(myContrasts)){
	e1<-myContrasts$exp1_mods[i]
	e2<-myContrasts$exp2_mods[i]
	e1_genes<-names(exp1_moduleLabels[exp1_moduleColors==e1])
	e2_genes<-names(exp2_moduleLabels[exp2_moduleColors==e2])
	vd<-venn.diagram(x=list(e1_genes,e2_genes),category.names=c(paste("Experiment 1 -",e1),paste("Experiment 2 -",e2)),filename=paste("venn_diagram",e1,"vs",e2,ct,"png", sep="."),imagetype="png",output=TRUE)
}

#Also test for correlation (spearman?) of module membership eigenvalues between the two experiments (ex: if green associated w/ dom Y in exp1 and brown associated w/ dom Y in exp2, test for correlation between exp1 green membership and exp2 brown membership for each gene)
#Maybe identify genes w/ high correlation and high eigenvalues in both and test if these are enriched for certain GO terms or sex chromosomes?


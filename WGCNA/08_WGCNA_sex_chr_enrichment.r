#PURPOSE: Test if interesting modules are over or under-enriched for sex chromosomes

args<-commandArgs(TRUE)
if(length(args) != 1){
        stop("Missing command line arguments.\nArgument 1: dataset")
}

dataset<-args[1]

print(paste0("Running sex chromosome enrichment analysis for ", dataset))

library(WGCNA)
#library(anRichment)
library(biomaRt)

options(stringsAsFactors = FALSE)
#enableWGCNAThreads(nThreads=8)

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

#Get sex chromosome genes
XchrGenes_df<-myIDs[which(myIDs$chromosome_name=="X"), ]
XchrGenes_names<-unique(XchrGenes_df$ensembl_gene_id)
XchrGenes_inDataset<-XchrGenes_names[which(XchrGenes_names %in% names(moduleLabels))]
YchrGenes_df<-myIDs[which(myIDs$chromosome_name=="Y"), ]
YchrGenes_names<-unique(YchrGenes_df$ensembl_gene_id)
YchrGenes_inDataset<-YchrGenes_names[which(YchrGenes_names %in% names(moduleLabels))]

#Enrichment tests for each module and sex chromosome
#interesting_mods<-c("blue","turquoise","black","green","yellow","yellowgreen")
#Fisher's exact test for each module
#Xfisher_ps<-c()
#Yfisher_ps<-c()
#for(m in interesting_mods){
#	print(m)
#	mod_genes<-names(moduleLabels)[which(moduleColors==m)]
#	print("X chromosome:")
#	X_only<-XchrGenes_inDataset[which(!(XchrGenes_inDataset %in% mod_genes))]
#	mod_only<-mod_genes[which(!(mod_genes %in% XchrGenes_inDataset))]
#	both_genes<-intersect(mod_genes, XchrGenes_inDataset)
#	Xmod_union<-union(mod_genes, XchrGenes_inDataset)
#	none_genes<-names(moduleLabels)[which(!(names(moduleLabels) %in% Xmod_union))]
#	myTable<-cbind(rbind(length(both_genes), length(X_only)), rbind(length(mod_only), length(none_genes)))
#	rownames(myTable)<-c("in_module","NOT_in_module")
#	colnames(myTable)<-c("on_X","NOT_on_X")
#	print(myTable)
#	result<-fisher.test(myTable)
#	print(result)
#	Xfisher_ps<-c(Xfisher_ps,result$p.value)
#	print("Y chromosome:")
#	Y_only<-YchrGenes_inDataset[which(!(YchrGenes_inDataset %in% mod_genes))]
#        mod_only<-mod_genes[which(!(mod_genes %in% YchrGenes_inDataset))]
#        both_genes<-intersect(mod_genes, YchrGenes_inDataset)
#        Ymod_union<-union(mod_genes, YchrGenes_inDataset)
#        none_genes<-names(moduleLabels)[which(!(names(moduleLabels) %in% Ymod_union))]
#	myTable<-cbind(rbind(length(both_genes), length(Y_only)), rbind(length(mod_only), length(none_genes)))
#        rownames(myTable)<-c("in_module","NOT_in_module")
#        colnames(myTable)<-c("on_Y","NOT_on_Y")
#        print(myTable)
#        result<-fisher.test(myTable)
#        print(result)
#        Yfisher_ps<-c(Yfisher_ps,result$p.value)
#}
#
#Xcor_ps_f<-p.adjust(Xfisher_ps)
#names(Xcor_ps_f)<-interesting_mods
#print("FDR-corrected p-values for Fisher Exact Tests, X chromosome:")
#print(Xcor_ps_f)
#Ycor_ps_f<-p.adjust(Yfisher_ps)
#names(Ycor_ps_f)<-interesting_mods
#print("FDR-corrected p-values for Fisher Exact Tests, Y chromosome:")
#print(Ycor_ps_f)

#Binomial Test
Xbinom_ps<-c()
Ybinom_ps<-c()
for(m in interesting_mods){
        print(m)
        mod_genes<-names(moduleLabels)[which(moduleColors==m)]
        print("X chromosome:")
        X_only<-XchrGenes_inDataset[which(!(XchrGenes_inDataset %in% mod_genes))]
        mod_only<-mod_genes[which(!(mod_genes %in% XchrGenes_inDataset))]
        both_genes<-intersect(mod_genes, XchrGenes_inDataset)
        X_prop<-length(which(names(moduleLabels) %in% XchrGenes_inDataset)) / length(names(moduleLabels))
        result<-binom.test(length(both_genes), length(mod_genes), p=X_prop, alternative="two.sided", conf.level=0.95)
        print(result)
        Xbinom_ps<-c(Xbinom_ps, result$p.value)
	print("Y chromosome:")
        Y_only<-YchrGenes_inDataset[which(!(YchrGenes_inDataset %in% mod_genes))]
        mod_only<-mod_genes[which(!(mod_genes %in% YchrGenes_inDataset))]
        both_genes<-intersect(mod_genes, YchrGenes_inDataset)
	Y_prop<-length(which(names(moduleLabels) %in% YchrGenes_inDataset)) / length(names(moduleLabels))
        result<-binom.test(length(both_genes), length(mod_genes), p=Y_prop, alternative="two.sided", conf.level=0.95)
        print(result)
        Ybinom_ps<-c(Ybinom_ps, result$p.value)
}

Xcor_ps_b<-p.adjust(Xbinom_ps)
names(Xcor_ps_b)<-interesting_mods
print("FDR-corrected p-values for Binomial Tests, X chromosome:")
print(Xcor_ps_b)
Ycor_ps_b<-p.adjust(Ybinom_ps)
names(Ycor_ps_b)<-interesting_mods
print("FDR-corrected p-values for Binomial Tests, Y chromosome:")
print(Ycor_ps_b)

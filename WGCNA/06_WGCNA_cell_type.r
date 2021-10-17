#PURPOSE: Test if modules associated w/ early and late cell types also contain genes and/or have high module membership for genes expressed and induced in those cell types

args<-commandArgs(TRUE)
if(length(args) != 2){
        stop("Missing command line arguments.\nArgument 1: dataset (options: all, exp1, exp2)\nArgument 2: signed or unsigned? (options: TRUE, FALSE)")
}

dataset<-args[1] #all #exp1 #exp2
signed<-args[2] #Perform signed network analysis?

print(paste0("Running cell type association analysis for ", dataset, ", signed: ", signed))

library(VennDiagram)
library(gridExtra)
library(ggplot2)
library(biomaRt)
library(WGCNA)

load(paste("05-WGCNA",dataset,"RData", sep="."))
#load(paste0(dataset,"Data-block.1.RData"))

if(dataset=="exp1"){
	earlyMod<-"turquoise"
	lateMod<-"blue"
} else if(dataset=="exp2"){
	earlyMod<-"turquoise"
        lateMod<-"blue"
} else{
        stop("Invalid dataset: must be 'exp1' or 'exp2'")
}

earlyMod_genes<-names(moduleLabels[moduleColors==earlyMod])
lateMod_genes<-names(moduleLabels[moduleColors==lateMod])

#Read in expressed and induced gene lists from spermatogenesis molecular evolution paper
exp_early<-scan("/mnt/beegfs/ek112884/mus_expression_analysis/MULTI_MAP/gene_list_eve_LZ_edgeR_wholeGenome.ensemblOrthos.txt", what=character())
ind_early<-scan("/mnt/beegfs/ek112884/mus_expression_analysis/MULTI_MAP/gene_list_LZinduced_edgeR_wholeGenome.ensemblOrthos.txt", what=character())
exp_late<-scan("/mnt/beegfs/ek112884/mus_expression_analysis/MULTI_MAP/gene_list_eve_RS_edgeR_wholeGenome.ensemblOrthos.txt", what=character())
ind_late<-scan("/mnt/beegfs/ek112884/mus_expression_analysis/MULTI_MAP/gene_list_RSinduced_edgeR_wholeGenome.ensemblOrthos.txt", what=character())

#Some sanity checks to make sure early and late module have a lot of overlap with genes induced in these cell types in previous work
print("# genes in expressed early, induced early, and in early module:")
print(paste(length(exp_early), length(ind_early), length(earlyMod_genes)))
print("# genes expressed early and in early module:")
print(length(intersect(exp_early, earlyMod_genes)))
print("# genes induced early and in early module:")
print(length(intersect(ind_early, earlyMod_genes)))

print("# genes in expressed late, induced late, and in late module:")
print(paste(length(exp_late), length(ind_late), length(lateMod_genes)))
print("# genes expressed late and in late module:")
print(length(intersect(exp_late, lateMod_genes)))
print("# genes induced late and in late module:")
print(length(intersect(ind_late, lateMod_genes)))

#Venn diagrams
vd_early<-venn.diagram(x=list(exp_early,ind_early,earlyMod_genes),category.names=c("Expressed Early","Induced Early","Early Module"),filename=paste("venn_diagram",dataset,"LZ.png", sep="."),imagetype="png",output=TRUE)
vd_late<-venn.diagram(x=list(exp_late,ind_late,lateMod_genes),category.names=c("Expressed Late","Induced Late","Late Module"),filename=paste("venn_diagram",dataset,"RS.png", sep="."),imagetype="png",output=TRUE)

#Do genes induced early and late also have strongest module membership in these modules?
early_mm<-myKME[which(rownames(myKME) %in% earlyMod_genes), earlyMod]
earlyMod_ind<-c(earlyMod_genes %in% ind_early)
earlyMod_exp<-c(earlyMod_genes %in% exp_early)
early_df<-as.data.frame(cbind(MM=early_mm, is_ind=earlyMod_ind, is_exp=earlyMod_exp))
earlyInd_p<-ggplot(early_df, aes(x=as.factor(is_ind), y=as.numeric(as.character(MM)))) + geom_violin() + geom_boxplot(width=0.1)
earlyInd_p<-earlyInd_p + labs(title="Module membership in 'early' module for genes induced early", x="Induced Early?", y="Module Membership (eigengene-based connectivity)")
earlyInd_p<-earlyInd_p + theme(axis.text.y = element_text(size=20))
earlyInd_p<-earlyInd_p + theme_minimal()
earlyExp_p<-ggplot(early_df, aes(x=as.factor(is_exp), y=as.numeric(as.character(MM)))) + geom_violin() + geom_boxplot(width=0.1)
earlyExp_p<-earlyExp_p + labs(title="Module membership in 'early' module for genes expressed early", x="Expressed Early?", y="Module Membership (eigengene-based connectivity)")
earlyExp_p<-earlyExp_p + theme(axis.text.y = element_text(size=20))
earlyExp_p<-earlyExp_p + theme_minimal()

late_mm<-myKME[which(rownames(myKME) %in% lateMod_genes), lateMod]
lateMod_ind<-c(lateMod_genes %in% ind_late)
lateMod_exp<-c(lateMod_genes %in% exp_late)
late_df<-as.data.frame(cbind(MM=late_mm, is_ind=lateMod_ind, is_exp=lateMod_exp))
lateInd_p<-ggplot(late_df, aes(x=as.factor(is_ind), y=as.numeric(as.character(MM)))) + geom_violin() + geom_boxplot(width=0.1)
lateInd_p<-lateInd_p + labs(title="Module membership in 'late' module for genes induced late", x="Induced Late?", y="Module Membership (eigengene-based connectivity)")
lateInd_p<-lateInd_p + theme(axis.text.y = element_text(size=20))
lateInd_p<-lateInd_p + theme_minimal()
lateExp_p<-ggplot(late_df, aes(x=as.factor(is_exp), y=as.numeric(as.character(MM)))) + geom_violin() + geom_boxplot(width=0.1)
lateExp_p<-lateExp_p + labs(title="Module membership in 'late' module for genes expressed late", x="Expressed Late?", y="Module Membership (eigengene-based connectivity)")
lateExp_p<-lateExp_p + theme(axis.text.y = element_text(size=20))
lateExp_p<-lateExp_p + theme_minimal()

plots.list<-list(earlyInd_p, earlyExp_p, lateInd_p, lateExp_p)
plots<-marrangeGrob(plots.list,nrow=1,ncol=1)
ggsave(paste("moduleMembership_violin.byCellType",dataset,"pdf", sep="."),plots,width=11,height=8.5,units="in")

print("Here are the Wilcox test results:")
print("Expressed early:")
wilcox.test(as.numeric(as.character(early_df[which(early_df$is_ind==0),]$MM)),as.numeric(as.character(early_df[which(early_df$is_ind==1),]$MM)))
print("Induced early:")
wilcox.test(as.numeric(as.character(early_df[which(early_df$is_exp==0),]$MM)),as.numeric(as.character(early_df[which(early_df$is_exp==1),]$MM)))
print("Expressed late:")
wilcox.test(as.numeric(as.character(late_df[which(late_df$is_ind==0),]$MM)),as.numeric(as.character(late_df[which(late_df$is_ind==1),]$MM)))
print("Induced late:")
wilcox.test(as.numeric(as.character(late_df[which(late_df$is_exp==0),]$MM)),as.numeric(as.character(late_df[which(late_df$is_exp==1),]$MM)))

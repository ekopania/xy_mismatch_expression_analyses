#PURPOSE: Test if ampliconic gene families are associated with particular modules or have high module membership in interesting modules associated with cross types

args<-commandArgs(TRUE)
if(length(args) != 2){
        stop("Missing command line arguments.\nArgument 1: dataset (options: all, exp1, exp2)\nArgument 2: signed or unsigned? (options: TRUE, FALSE)")
}

dataset<-args[1] #all #exp1 #exp2
signed<-args[2] #FALSE #Perform signed network analysis?

print(paste0("Running ampliconic gene family association analysis for ", dataset, ", signed: ", signed))

library(VennDiagram)
library(gridExtra)
library(ggplot2)
library(biomaRt)
library(WGCNA)

load(paste("05-WGCNA",dataset,"RData", sep="."))
#load(paste0(dataset,"Data-block.1.RData"))

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
        stop("Invalid dataset: must be 'all', 'exp1', 'exp1_LZ', 'exp1_RS', 'exp2', 'exp2_LZ', or 'exp2_RS'")
}

#See which modules ampliconic genes fall into
ens_mus<-useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl")
my_genes<-getBM(attributes=c('ensembl_gene_id','external_gene_name'), mart=ens_mus)

#amplicon_list<-c("astx", "asty", "slx", "slxl1", "sly", "srsx", "srsy", "sstx", "ssty1", "ssty2", "atakusan", "speer")
amplicon_list<-c("slx", "slxl1", "sly", "sstx", "ssty1", "ssty2", "atakusan", "speer")
amp_geneIDs<-c()
amp_family<-c()
amp_modules<-c()
for(a in amplicon_list){
        print(a)
        this_table<-read.table(paste0("../GENE_FAMILY_FILES/",a,".pID97.gene_family_paralogs.bed"),header=FALSE, sep="\t")
        paralogs0<-unlist(sapply(this_table$V6, function(x) my_genes$ensembl_gene_id[which(my_genes$external_gene_name==x)]))
        paralogs<-paralogs0[which(paralogs0 %in% names(moduleLabels))]
        mod_nums<-moduleLabels[paralogs]
        assign(paste0(a, "_modules"), unlist(sapply(mod_nums, function(x) label_to_color$col[which(label_to_color$lab==x)])))
        print(get(paste0(a, "_modules")))
        amp_geneIDs<-c(amp_geneIDs, paralogs)
        amp_family<-c(amp_family, rep(a, length(paralogs)))
        amp_modules<-c(amp_modules, get(paste0(a, "_modules")))
}

print(length(amp_geneIDs))
print(length(amp_family))
print(length(amp_modules))
amplicons_df<-as.data.frame(cbind(geneID=amp_geneIDs, family=amp_family, module=amp_modules))

#Get module membership for ampliconic genes for interesting modules
#interesting_mods<-c("blue","turquoise","black","green","yellow","yellowgreen")
amp_MM<-c()
for(g in amplicons_df$geneID){
        my_MM<-myKME[g, interesting_mods]
        amp_MM<-rbind(amp_MM, my_MM)
}

amplicons_df_withMM<-as.data.frame(cbind(amplicons_df, amp_MM))
write.table(amplicons_df_withMM, paste("amplicon_modules",dataset,"txt", sep="."), row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)

#Violin plots: Median module membership for each "interesting module", separated by ampliconic family on X axis
plots.list<-list()
for(m in interesting_mods){
        print(m)
        p<-ggplot(amplicons_df_withMM, aes_string(x="family", y=m)) + geom_violin() + geom_boxplot(width=0.1)
        p<-p + labs(title=paste("Module membership in",m,"module for ampliconic families"), x="Gene Family", y="Module Membership (eigengene-based connectivity)")
        p<-p + theme(axis.text.y = element_text(size=20))
        p<-p + theme_minimal()
        plots.list[[length(plots.list)+1]]<-p
}
plots<-marrangeGrob(plots.list,nrow=1,ncol=1)
ggsave(paste("moduleMembership_violin.byGeneFamily",dataset,"pdf", sep="."),plots,width=11,height=8.5,units="in")
#To make sure nothing weird happened with ggplot:
warnings()

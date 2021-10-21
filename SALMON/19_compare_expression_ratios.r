#PURPOSE: Compare ratio of median X chromosome expression to median autosome expression for different cross types (repeat for Y)

library(biomaRt)
library(ggplot2)
library(ggbeeswarm)

dataset<-"Yintro_exp2"
cell_type<-"LZ"

#Read in data
fpkm<-read.table(paste("fpkm_filtered_table", dataset, cell_type, "97.txt", sep="."), header=TRUE)

#Get gene coordinates
ens_mus<-useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl")
all_genes<-getBM(attributes=c('ensembl_gene_id','chromosome_name','start_position','end_position'), mart=ens_mus)
colnames(all_genes)<-c("geneID", "seqnames", "start", "end")
print(head(all_genes))

auto_genes<-all_genes$geneID[which(all_genes$seqnames %in% c(1:19))]
x_genes<-all_genes$geneID[which(all_genes$seqnames=="X")]
y_genes<-all_genes$geneID[which(all_genes$seqnames=="Y")]

fpkm_auto<-fpkm[which(rownames(fpkm) %in% auto_genes),]
fpkm_x<-fpkm[which(rownames(fpkm) %in% x_genes),]
fpkm_y<-fpkm[which(rownames(fpkm) %in% y_genes),]

if(dataset=="Yintro_exp1"){
	cross_types<-c("CCPP","CPLY","WLPY","WWLL")
} else if(dataset=="Yintro_exp2"){
	cross_types<-c("PLPY","PPLL","LLPP","LPLY")
} else{
	stop("Invalid dataset; must be 'Yintro_exp1' or 'Yintro_exp2'")
}

auto_meds<-c()
x_meds<-c()
y_meds<-c()
for(c in cross_types){
	#Get median fpkm for each chromosome type for each sample
	temp_auto<-fpkm_auto[, grepl(c, colnames(fpkm_auto))]
	temp_x<-fpkm_x[, grepl(c, colnames(fpkm_x))]
	temp_y<-fpkm_y[, grepl(c, colnames(fpkm_y))]
	a<-apply(as.matrix(temp_auto), 2, median)
	print(a)
	x<-apply(as.matrix(temp_x), 2, median)
	print(x)
	y<-apply(as.matrix(temp_y), 2, median)
        print(y)
	#DF with rows as samples and columns as medians and cross type
	temp_auto_meds<-as.data.frame(cbind(medians=a, group=rep(c, length(a))))
	rownames(temp_auto_meds)<-names(a)
	temp_x_meds<-as.data.frame(cbind(medians=x, group=rep(c, length(x))))
        rownames(temp_x_meds)<-names(x)
	temp_y_meds<-as.data.frame(cbind(medians=y, group=rep(c, length(y))))
        rownames(temp_y_meds)<-names(y)
	#Append DF to include all cross types
	auto_meds<-as.data.frame(rbind(auto_meds, temp_auto_meds))
	x_meds<-as.data.frame(rbind(x_meds, temp_x_meds))
	y_meds<-as.data.frame(rbind(y_meds, temp_y_meds))
}
print("Medians:")
print(auto_meds)
print(x_meds)
print(y_meds)

#Get ratios
ratio_xa<-as.numeric(as.character(x_meds$medians))/as.numeric(as.character(auto_meds$medians))
ratio_xa_df<-as.data.frame(cbind(ratio=ratio_xa, group=auto_meds$group))
rownames(ratio_xa_df)<-rownames(auto_meds)
ratio_xa_df$ratio<-as.numeric(as.character(ratio_xa_df$ratio))
ratio_ya<-as.numeric(as.character(y_meds$medians))/as.numeric(as.character(auto_meds$medians))
ratio_ya_df<-as.data.frame(cbind(ratio=ratio_ya, group=auto_meds$group))
rownames(ratio_ya_df)<-rownames(auto_meds)
ratio_ya_df$ratio<-as.numeric(as.character(ratio_ya_df$ratio))

#Compare ratios
print("Wilcoxon test results, X:auto")
print(pairwise.wilcox.test(ratio_xa_df$ratio, ratio_xa_df$group, method="fdr"))
print("Wilcoxon test results, Y:auto")
print(pairwise.wilcox.test(ratio_ya_df$ratio, ratio_ya_df$group, method="fdr"))

print("t-test results, X:auto")
print(pairwise.t.test(ratio_xa_df$ratio, ratio_xa_df$group, method="fdr"))
print("t-test results, Y:auto")
print(pairwise.t.test(ratio_ya_df$ratio, ratio_ya_df$group, method="fdr"))

#Plot boxplots
print("Plotting...")
pdf(paste("chr_type_expression_ratios.boxplots", dataset, cell_type, "pdf", sep="."), onefile=TRUE)
px<-ggplot(ratio_xa_df, aes(x=group, y=ratio)) + geom_boxplot() + geom_quasirandom()
px<-px + labs(title="Median expression ratio X:auto", x="", y="X:auto") + theme_minimal()
px<-px + theme(axis.text=element_text(size=18), axis.title=element_text(size=21), plot.title=element_text(size=28))
px<-px + scale_x_discrete(limits=cross_types) + ylim(0,1)
py<-ggplot(ratio_ya_df, aes(x=group, y=ratio)) + geom_boxplot() + geom_quasirandom()
py<-py + labs(title="Median expression ratio Y:auto", x="", y="Y:auto") + theme_minimal()
py<-py + theme(axis.text=element_text(size=18), axis.title=element_text(size=21), plot.title=element_text(size=28))
py<-py + scale_x_discrete(limits=cross_types) + ylim(0,1)
print(px)
print(py)
dev.off()

print("Done with 19_compare_expression_ratios.r")

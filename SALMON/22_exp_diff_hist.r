#PURPOSE: Plot logFC values as histograms for X and Y on opposite sides

library(biomaRt)
library(vegan)
library(ggplot2)

dataset<-"Yintro_exp2"
cell_type<-"RS"
procoOnly<-FALSE

#Read in data
fpkm_raw<-read.table(paste("fpkm_filtered_table",dataset,cell_type,"97.txt", sep="."), header=TRUE)
#Normalize expression data
fpkm_data<-decostand(fpkm_raw, "normalize")

if(dataset=="Yintro_exp1"){
        #Set up with sex chr mismatch always first
        cross_types<-list(c("CPLY","CCPP"), c("WLPY","CCPP"), c("WLPY","WWLL"), c("CPLY","WWLL"))
        names(cross_types)<-c("musX","musY","domX","domY")
} else if(dataset=="Yintro_exp2"){
        cross_types<-list(c("PPLL","PLPY"), c("LLPP","PLPY"), c("LLPP","LPLY"), c("PPLL","LPLY"))
        names(cross_types)<-c("musX","musY","domX","domY")
} else{
        stop("Invalid dataset: must be 'Yintro_exp1' or 'Yintro_exp2'")
}

print("Here are the cross type comparisons:")
print(cross_types)

#Separate expressed genes by chr type
ens_mus<-useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl")
all_genes<-getBM(attributes=c('ensembl_gene_id','chromosome_name','start_position','end_position'), mart=ens_mus)
colnames(all_genes)<-c("geneID", "chr", "start", "end")
print(head(all_genes))
auto_genes<-all_genes$geneID[which(all_genes$chr %in% c(1:19))]
x_genes<-all_genes$geneID[which(all_genes$chr=="X")]
y_genes<-all_genes$geneID[which(all_genes$chr=="Y")]

fpkm_auto<-fpkm_data[which(rownames(fpkm_data) %in% auto_genes),]
fpkm_x<-fpkm_data[which(rownames(fpkm_data) %in% x_genes),]
fpkm_y<-fpkm_data[which(rownames(fpkm_data) %in% y_genes),]

pdf(paste("expression_difference_hist",dataset,cell_type,"pdf", sep="."), onefile=TRUE)
for(c in names(cross_types)){
	print(paste("Working on cross type:", c))
	#Get fpkms for each cross type in comparison
	fpkm1<-fpkm_data[, grepl(cross_types[[c]][1],colnames(fpkm_data))]
	fpkm2<-fpkm_data[, grepl(cross_types[[c]][2],colnames(fpkm_data))]
	#Separate by chromosome
	fpkm1_auto<-fpkm1[which(rownames(fpkm1) %in% auto_genes),]
	fpkm1_x<-fpkm1[which(rownames(fpkm1) %in% x_genes),]
	fpkm1_y<-fpkm1[which(rownames(fpkm1) %in% y_genes),]
	fpkm2_auto<-fpkm2[which(rownames(fpkm2) %in% auto_genes),]
        fpkm2_x<-fpkm2[which(rownames(fpkm2) %in% x_genes),]
        fpkm2_y<-fpkm2[which(rownames(fpkm2) %in% y_genes),]
	#get difference in mean expression
	#auto
	stopifnot(all.equal(rownames(fpkm1_auto), rownames(fpkm2_auto)))
	fpkm1_auto_means<-rowMeans(fpkm1_auto)
	fpkm2_auto_means<-rowMeans(fpkm2_auto)
	auto_diff<-fpkm1_auto_means - fpkm2_auto_means
	print(length(auto_diff))
	print(head(auto_diff))
	#X
	stopifnot(all.equal(rownames(fpkm1_x), rownames(fpkm2_x)))
        fpkm1_x_means<-rowMeans(fpkm1_x)
        fpkm2_x_means<-rowMeans(fpkm2_x)
        x_diff<-fpkm1_x_means - fpkm2_x_means
        print(length(x_diff))
        print(head(x_diff))
	#Y
	stopifnot(all.equal(rownames(fpkm1_y), rownames(fpkm2_y)))
        fpkm1_y_means<-rowMeans(fpkm1_y)
        fpkm2_y_means<-rowMeans(fpkm2_y)
        y_diff<-fpkm1_y_means - fpkm2_y_means
        print(length(y_diff))
        print(head(y_diff))
	#plot
	df_auto<-as.data.frame(cbind(auto=auto_diff))
	df_x<-as.data.frame(cbind(X=x_diff))
	df_y<-as.data.frame(cbind(Y=y_diff))
	#print(head(df_auto))
	#Setting y=..count../sum(..count..)) should make the y=axis of histogram proportion of total
	p<-ggplot(df_auto, aes(x)) + geom_histogram(aes(x=auto, y=..count../sum(..count..)), fill="black", binwidth=diff(range(df_auto$auto))/30)
	p<-p + geom_histogram(data=df_x, aes(x=X, y=-..count../sum(..count..)), fill="lightgrey", alpha=0.5, binwidth=diff(range(df_auto$auto))/30)
	p<-p + geom_histogram(data=df_y, aes(x=Y, y=-..count../sum(..count..)), fill="darkgrey", alpha=0.5, binwidth=diff(range(df_auto$auto))/30)
	p<-p + labs(title=paste("Histogram of normalized gene expression:\n",cross_types[[c]][1],"vs",cross_types[[c]][2]), x="Normalized FPKM Difference", y="Proportion of Genes")
	p<-p + theme_minimal() + xlim(-0.6, 0.6) + ylim(-0.3, 0.3)
	p<-p + theme(axis.text=element_text(size=18), axis.title=element_text(size=21), plot.title=element_text(size=22))
	print(p + coord_flip())
	#Kolmogorov-Smirnov tests
	print("autos vs X:")
	print(ks.test(df_auto$auto, df_x$X))
	print("autos vs Y:")
	print(ks.test(df_auto$auto, df_y$Y))
	print("X vs Y:")
	print(ks.test(df_x$X, df_y$Y))
}
dev.off()

warnings()

print("Done with 22_exp_diff_hist.r")

#PURPOSE: Get median expression levels on the sex chromosomes for different genotypes and cell types and compare

#change these to test different parameters
#Also, some of these gene families (srsx I believe) have FPKM <1 for every paralog across all samples
min_rpkm<-1
min_samples<-4
#min_logFC<-0
#induced_cutoff<-2

#Choose which dataset to look at
dataset<-"LarsonEtal"
#cell_type<-"LZandRS"
#cell_type<-"LZ"
cell_type<-"RS"
#Percent ID threshold for considering genes paralogs/part of the same gene family
#paralog_thresh<-95
induced<-TRUE
testis_specific<-FALSE

print(paste("Comparing sex chromosome gene expression for:", dataset, cell_type))
print(paste("Induced?", induced))
print(paste("Testis-specific?", testis_specific))

print("Loading R libraries...")
library(edgeR)
#library(EnsDb.Mmusculus.v79)
library(ggplot2)
library(gridExtra)
library(biomaRt)
library(matrixStats)

print("Setting up DGE object...")
mydata<-read.table("salmon_output.counts.txt")
mylengths<-read.table("salmon_output.lengths.txt")
#Groups: 1=CCPP_LZ, 2=CCPP_RS, 3=WWLL_LZ 4=WWLL_RS, 5=LLPP_LZ, 6=LLPP_RS, 7=PPLL_LZ, 8=PPLL_RS
if(dataset=="LarsonEtal"){
        if(cell_type=="LZandRS"){
                these_data<-mydata
                these_lengths<-mylengths
                mygroups<-c(1,2,1,2,1,2,6,5,6,6,5,5,8,7,8,7,8,7,4,3,4,3,4,3)
                geno_order<-c("CCPP_LZ","CCPP_RS","WWLL_LZ","WWLL_RS","LLPP_LZ","LLPP_RS","PPLL_LZ","PPLL_RS")
        } else if(cell_type=="LZ"){
                these_data<-mydata[,grepl("LZ",colnames(mydata))]
                these_lengths<-mylengths[,grepl("LZ",colnames(mylengths))]
                mygroups<-c(1,1,1,5,5,5,7,7,7,3,3,3)
                geno_order<-c("CCPP_LZ","WWLL_LZ","LLPP_LZ","PPLL_LZ")
        } else if(cell_type=="RS"){
                these_data<-mydata[,grepl("RS",colnames(mydata))]
                these_lengths<-mylengths[,grepl("RS",colnames(mylengths))]
                mygroups<-c(2,2,2,6,6,6,8,8,8,4,4,4)
                geno_order<-c("CCPP_RS","WWLL_RS","LLPP_RS","PPLL_RS")
        } else{
                stop("Invalid cell type")
        }
} else{
        stop("Invalid dataset")
}

#NEED TO ACCOUNT FOR GENE LENGTHS - they may be different among samples because we originally came from txpt level counts
#SEE TXIMPORT DOCUMENTATION: http://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html#Downstream_DGE_in_Bioconductor
normMat<-these_lengths/exp(rowMeans(log(these_lengths)))
normCts<-these_data/normMat
eff.lib <- calcNormFactors(normCts) * colSums(normCts)
normMat <- sweep(normMat, 2, eff.lib, "*")
normMat <- log(normMat)
num_normMat<-data.matrix(normMat) #Convert data frame to numeric matrix to keep EdgeR happy
#print(dim(these_data))
#print(length(mygroups))
#print(head(these_data))
#print(is.numeric(these_data))
myDGE<-DGEList(counts=these_data, group=mygroups)
myDGE<-scaleOffset(myDGE, num_normMat)
print("Dimensions of myDGE counts:")
print(dim(myDGE$counts))
print(myDGE$counts[1:5,1:5])
print("myDGE samples - CHECK GROUP LABELS:")
print(myDGE$samples)
#This will be useful later to put meaningful labels on graphs, instead of the numeric group labels EdgeR uses
group_to_cross<-c(1, 2, 3, 4, 5, 6, 7, 8)
names(group_to_cross)<-c("CCPP_LZ", "CCPP_RS", "WWLL_LZ", "WWLL_RS", "LLPP_LZ", "LLPP_RS", "PPLL_LZ", "PPLL_RS")

print("Filtering by expression level...")
#EdgeR rpkm function seems to require same gene length across all samples; taking median of tximport output gene lengths and using that
new_lengths<-rowMedians(as.matrix(these_lengths))
names(new_lengths)<-rownames(these_lengths)
print(head(new_lengths))
print(dim(myDGE$offset))
print(myDGE$samples$lib.size)
my_rpkm<-rpkm(myDGE, gene.length=new_lengths)
print(dim(my_rpkm))

keep<-rowSums(my_rpkm > min_rpkm) >= min_samples
keep[is.na(keep)]<-FALSE #NAs show up for gene IDs that aren't in myDGE; replace then with false otherwise myDGE gets confused

myDGE<-myDGE[keep, , keep.lib.sizes=FALSE]
myDGE<-calcNormFactors(myDGE)
print(myDGE$samples)

print("Dimensions of filtered myDGE counts:")
print(dim(myDGE$counts))

#Calculate tpm for ALL GENES and then extract sex chromosome genes
print("Calculating tpm...")
newOrdered_lengths<-new_lengths[match(rownames(myDGE$counts),names(new_lengths))] #filtered myDGE order
new_rpkm<-rpkm(myDGE,gene.length=newOrdered_lengths) #Based on filtered gene set
print(dim(new_rpkm))
print(new_rpkm[1:5,1:5])
divLen<-apply(myDGE$counts, 2, function(x) x/as.numeric(newOrdered_lengths)) #normalize counts by gene length
#Some checks:
#print(dim(myDGE$counts))
#print(dim(divLen))
#print(myDGE$counts[1:5,1:5])
#print(head(newOrdered_lengths))
#print(divLen[1:5,1:5])
scaling_factors<-colSums(divLen)/1000000
print(scaling_factors)
tpm<-apply(divLen, 1, function(x) x/scaling_factors)
tpm_transformed<-t(tpm)
print(dim(tpm_transformed))
print(head(tpm_transformed))
#Get tpms for each sex chromosome
print("Getting sex chromosome genes...")
#edb<-EnsDb.Mmusculus.v79
#edb_y<-addFilter(edb, SeqNameFilter("Y"))
#y_genes<-genes(edb_y)
#edb_x<-addFilter(edb, SeqNameFilter("X"))
#x_genes<-genes(edb_x)
#edb_auto<-addFilter(edb, SeqNameFilter(c(1:19)))
#auto_genes<-genes(edb_auto)
ens_mus<-useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl")
all_genes<-getBM(attributes=c('ensembl_gene_id','external_gene_name','chromosome_name','start_position','end_position','strand','gene_biotype'), mart=ens_mus)
colnames(all_genes)<-c("gene_id","gene_name","seqnames","start","end","strand","gene_biotype")
#Subset to induced genes AFTER all normalization steps
if(induced){
	ind_genes<-scan(paste0("gene_list_",cell_type,"induced_edgeR.lessStrict.",dataset,".txt"), what=character())
	all_genes<-all_genes[which(all_genes$gene_id %in% ind_genes), ]
}
if(testis_specific){
	ts_genes<-scan("/mnt/beegfs/ek112884/mus_expression_analysis/chalmel_testis_specific.txt.gz", what=character())
	all_genes<-all_genes[which(all_genes$gene_id %in% ts_genes), ]
}
auto_genes<-all_genes[which(all_genes$seqnames %in% c(1:19)),]
x_genes<-all_genes[which(all_genes$seqnames=="X"),]
y_genes<-all_genes[which(all_genes$seqnames=="Y"),]

tpm_x<-tpm_transformed[which(rownames(tpm_transformed) %in% x_genes$gene_id),] #GRCm38, Ensembl ID
tpm_y<-tpm_transformed[which(rownames(tpm_transformed) %in% y_genes$gene_id),] #GRCm38, Ensembl ID
#tpm_x<-tpm_transformed[which(rownames(tpm_transformed) %in% x_genes$gene_name),] #GRCm39, gene name
#tpm_y<-tpm_transformed[which(rownames(tpm_transformed) %in% y_genes$gene_name),] #GRCm39, gene name

#Separate by genotypes - FIGURE OUT THE "MEDIAN OF MEDIANS" PROBLEM!!! Pretty sure I did this before; check Erica's papers

#Make separate dfs of median tpm values for each sex chr, where one column is sample, one is genotype, one is XY mismatch or not, and one is median tpm value

#OK edgeR has an rpkm by group function; going to try this first, then figure out how it's done, then do the same thing for tpm
print("Getting group FPKMs...")
group_rpkm<-rpkmByGroup(myDGE, group=myDGE$sample$group, gene.length=as.numeric(newOrdered_lengths))
print(head(group_rpkm))
colnames(group_rpkm)<-names(group_to_cross)[which(group_to_cross %in% colnames(group_rpkm))]
print(head(group_rpkm))

#print(dim(group_rpkm))
#print(head(group_rpkm))
group_rpkm_x<-group_rpkm[which(rownames(group_rpkm) %in% x_genes$gene_id),] #GRCm38, Ensembl ID
group_rpkm_y<-group_rpkm[which(rownames(group_rpkm) %in% y_genes$gene_id),] #GRCm38, Ensembl ID
#group_rpkm_x<-group_rpkm[which(rownames(group_rpkm) %in% x_genes$gene_name),] #GRCm39, gene name
#group_rpkm_y<-group_rpkm[which(rownames(group_rpkm) %in% y_genes$gene_name),] #GRCm39, gene name
print(head(group_rpkm_x))
print(head(group_rpkm_y))

#set up X only df
rpkm_value<-as.vector(group_rpkm_x)
geno<-c(rep(colnames(group_rpkm_x)[1], nrow(group_rpkm_x)), rep(colnames(group_rpkm_x)[2], nrow(group_rpkm_x)), rep(colnames(group_rpkm_x)[3], nrow(group_rpkm_x)), rep(colnames(group_rpkm_x)[4], nrow(group_rpkm_x)))
#mismatch<-grepl("CCPPLY|WWLLPY|LLPP_|PPLL_",geno)
rpkm_x_df<-as.data.frame(cbind(rpkm_value, geno))
rpkm_x_df$rpkm_value<-as.numeric(as.character(rpkm_x_df$rpkm_value))
print(dim(rpkm_x_df))
print(head(rpkm_x_df))
#set up Y only df
rpkm_value<-as.vector(group_rpkm_y)
geno<-c(rep(colnames(group_rpkm_y)[1], nrow(group_rpkm_y)), rep(colnames(group_rpkm_y)[2], nrow(group_rpkm_y)), rep(colnames(group_rpkm_y)[3], nrow(group_rpkm_y)), rep(colnames(group_rpkm_y)[4], nrow(group_rpkm_y)))
#mismatch<-grepl("CCPPLY|WWLLPY|LLPP_|PPLL_",geno)
rpkm_y_df<-as.data.frame(cbind(rpkm_value, geno))
rpkm_y_df$rpkm_value<-as.numeric(as.character(rpkm_y_df$rpkm_value))
print(dim(rpkm_y_df))
print(head(rpkm_y_df))

print("Comparing overall expression levels on the sex chromosomes...")
#pairwise wilcoxon test
pwt_x<-pairwise.wilcox.test(rpkm_x_df$rpkm_value, rpkm_x_df$geno, p.adjust.method="fdr")
pwt_y<-pairwise.wilcox.test(rpkm_y_df$rpkm_value, rpkm_y_df$geno, p.adjust.method="fdr")
print(pwt_x)
print(pwt_y)

print("Plotting...")
p_x<-ggplot(rpkm_x_df, aes(x=geno, y=log(rpkm_value))) + geom_violin() + geom_boxplot(width=0.1)
p_x<-p_x + labs(title="FPKM on the X chromosome by cross type", x="Cross type", y="log(Group FPKM)")
p_x<-p_x + theme(axis.text.y = element_text(size=20))
p_x<-p_x + theme_minimal() #+ ylim(0,100)
#p_x<-p_x + scale_fill_manual(values=c("darkolivegreen","darkorchid4"))

p_y<-ggplot(rpkm_y_df, aes(x=geno, y=log(rpkm_value))) + geom_violin() + geom_boxplot(width=0.1)
p_y<-p_y + labs(title="FPKM on the Y chromosome by cross type", x="Cross type", y="log(Group FPKM)")
p_y<-p_y + theme(axis.text.y = element_text(size=20))
p_y<-p_y + theme_minimal()
#p_y<-p_y + scale_fill_manual(values=c("darkolivegreen","darkorchid4"))

plots.list<-list(p_x,p_y)
plots<-marrangeGrob(plots.list,nrow=1,ncol=1)
if(induced && testis_specific){
	ggsave(paste("fpkm_by_sexChr",dataset,cell_type,"induced_testisSpecific.pdf",sep="."), plots, width = 11, height = 8.5, units = "in")
} else if(induced){
	ggsave(paste("fpkm_by_sexChr",dataset,cell_type,"induced.pdf",sep="."), plots, width = 11, height = 8.5, units = "in")
} else if(testis_specific){
	ggsave(paste("fpkm_by_sexChr",dataset,cell_type,"testisSpecific.pdf",sep="."), plots, width = 11, height = 8.5, units = "in")
} else{
	ggsave(paste("fpkm_by_sexChr",dataset,cell_type,"pdf",sep="."), plots, width = 11, height = 8.5, units = "in")
}

##Repeat with normalization instead of log
#normalize<-function(x) (x- min(x))/(max(x) - min(x))
#p_x<-ggplot(rpkm_x_df, aes(x=geno, y=normalize(rpkm_value))) + geom_violin(aes(fill=mismatch)) + geom_boxplot(width=0.1)
#p_x<-p_x + labs(title="FPKM on the X chromosome by cross type", x="Cross type", y="Normalized group FPKM)")
#p_x<-p_x + theme(axis.text.y = element_text(size=20))
#p_x<-p_x + theme_minimal() #+ ylim(0,100)
#p_x<-p_x + scale_fill_manual(values=c("darkolivegreen","darkorchid4"))
#
#p_y<-ggplot(rpkm_y_df, aes(x=geno, y=normalize(rpkm_value))) + geom_violin(aes(fill=mismatch)) + geom_boxplot(width=0.1)
#p_y<-p_y + labs(title="FPKM on the Y chromosome by cross type", x="Cross type", y="Normalized group FPKM)")
#p_y<-p_y + theme(axis.text.y = element_text(size=20))
#p_y<-p_y + theme_minimal()
#p_y<-p_y + scale_fill_manual(values=c("darkolivegreen","darkorchid4"))
#
#plots.list<-list(p_x,p_y)
#plots<-marrangeGrob(plots.list,nrow=1,ncol=1)
#ggsave(paste("fpkm_by_sexChr.normalized",dataset,cell_type,"pdf",sep="."), plots, width = 11, height = 8.5, units = "in")

#Also plot tpm values against each other for Y intro vs control (to test if overall genes are more highly expresed in one vs other or if there are some over exp and other under exp)
if(induced && testis_specific){
	pdf(paste("fpkm_by_gene",dataset,cell_type,"induced_testisSpecific.pdf", sep="."), onefile=TRUE)
} else if(induced){
	pdf(paste("fpkm_by_gene",dataset,cell_type,"induced.pdf", sep="."), onefile=TRUE)
} else if(testis_specific){
	pdf(paste("fpkm_by_gene",dataset,cell_type,"testisSpecific.pdf", sep="."), onefile=TRUE)
} else{
	pdf(paste("fpkm_by_gene",dataset,cell_type,"pdf", sep="."), onefile=TRUE)
}
plot(x=as.numeric(as.character(group_rpkm_x[,1])), y=as.numeric(as.character(group_rpkm_x[,2])), main="Y-intro vs control - X chromosome", xlab=colnames(group_rpkm_x)[1], ylab=colnames(group_rpkm_x)[2])
abline(a=0, b=1)
plot(x=as.numeric(as.character(group_rpkm_x[,4])), y=as.numeric(as.character(group_rpkm_x[,3])), main="Y-intro vs control - X chromosome", xlab=colnames(group_rpkm_x)[4], ylab=colnames(group_rpkm_x)[3])
abline(a=0, b=1)
plot(x=as.numeric(as.character(group_rpkm_y[,1])), y=as.numeric(as.character(group_rpkm_y[,2])), main="Y-intro vs control - Y chromosome", xlab=colnames(group_rpkm_y)[1], ylab=colnames(group_rpkm_y)[2])
abline(a=0, b=1)
plot(x=as.numeric(as.character(group_rpkm_y[,4])), y=as.numeric(as.character(group_rpkm_y[,3])), main="Y-intro vs control - Y chromosome", xlab=colnames(group_rpkm_y)[4], ylab=colnames(group_rpkm_y)[3])
abline(a=0, b=1)
dev.off()
print(head(group_rpkm_y))

#print(head(which(group_rpkm_x[,1]==0)))
if(induced && testis_specific){
	write.table(group_rpkm, file=paste("group_fpkm",dataset,cell_type,"induced_testisSpecific.txt", sep="."), sep="\t", quote=FALSE)
} else if(induced){
	 write.table(group_rpkm, file=paste("group_fpkm",dataset,cell_type,"induced.txt", sep="."), sep="\t", quote=FALSE)
} else if(testis_specific){
	write.table(group_rpkm, file=paste("group_fpkm",dataset,cell_type,"testisSpecific.txt", sep="."), sep="\t", quote=FALSE)
} else{
	write.table(group_rpkm, file=paste("group_fpkm",dataset,cell_type,"txt", sep="."), sep="\t", quote=FALSE)
}

print("Protein coding only...")
#all_genes<-genes(edb)
#proco_genes<-all_genes$gene_name[which(all_genes$gene_biotype=="protein_coding")] #GRCm39
proco_genes<-all_genes$gene_id[which(all_genes$gene_biotype=="protein_coding")] #GRCm38
group_rpkm_x_proco<-group_rpkm_x[which(rownames(group_rpkm_x) %in% proco_genes),]
group_rpkm_y_proco<-group_rpkm_y[which(rownames(group_rpkm_y) %in% proco_genes),]
print(head(group_rpkm_x_proco))
print(head(group_rpkm_y_proco))

#set up X only df
rpkm_value<-as.vector(group_rpkm_x_proco)
geno<-c(rep(colnames(group_rpkm_x_proco)[1], nrow(group_rpkm_x_proco)), rep(colnames(group_rpkm_x_proco)[2], nrow(group_rpkm_x_proco)), rep(colnames(group_rpkm_x_proco)[3], nrow(group_rpkm_x_proco)), rep(colnames(group_rpkm_x_proco)[4], nrow(group_rpkm_x_proco)))
#mismatch<-grepl("CCPPLY|WWLLPY|LLPP_|PPLL_",geno)
rpkm_x_df_proco<-as.data.frame(cbind(rpkm_value, geno))
rpkm_x_df_proco$rpkm_value<-as.numeric(as.character(rpkm_x_df_proco$rpkm_value))
print(dim(rpkm_x_df_proco))
print(head(rpkm_x_df_proco))
#set up Y only df
rpkm_value<-as.vector(group_rpkm_y_proco)
geno<-c(rep(colnames(group_rpkm_y_proco)[1], nrow(group_rpkm_y_proco)), rep(colnames(group_rpkm_y_proco)[2], nrow(group_rpkm_y_proco)), rep(colnames(group_rpkm_y_proco)[3], nrow(group_rpkm_y_proco)), rep(colnames(group_rpkm_y_proco)[4], nrow(group_rpkm_y_proco)))
#mismatch<-grepl("CCPPLY|WWLLPY|LLPP_|PPLL_",geno)
rpkm_y_df_proco<-as.data.frame(cbind(rpkm_value, geno))
rpkm_y_df_proco$rpkm_value<-as.numeric(as.character(rpkm_y_df_proco$rpkm_value))
print(dim(rpkm_y_df_proco))
print(head(rpkm_y_df_proco))

print("Comparing protein coding gene expression levels on the sex chromosomes...")
#pairwise wilcoxon test
pwt_x<-pairwise.wilcox.test(rpkm_x_df_proco$rpkm_value, rpkm_x_df_proco$geno, p.adjust.method="fdr")
pwt_y<-pairwise.wilcox.test(rpkm_y_df_proco$rpkm_value, rpkm_y_df_proco$geno, p.adjust.method="fdr")
print(pwt_x)
print(pwt_y)

print("Plotting...")
p_x<-ggplot(rpkm_x_df_proco, aes(x=geno, y=log(rpkm_value))) + geom_violin() + geom_boxplot(width=0.1)
p_x<-p_x + labs(title="FPKM on the X chromosome by cross type - protein coding genes only", x="Cross type", y="log(Group FPKM)")
p_x<-p_x + theme(axis.text.y = element_text(size=20))
p_x<-p_x + theme_minimal() #+ ylim(0,100)
#p_x<-p_x + scale_fill_manual(values=c("darkolivegreen","darkorchid4"))

p_y<-ggplot(rpkm_y_df_proco, aes(x=geno, y=log(rpkm_value))) + geom_violin() + geom_boxplot(width=0.1)
p_y<-p_y + labs(title="FPKM on the Y chromosome by cross type - protein coding genes only", x="Cross type", y="log(Group FPKM)")
p_y<-p_y + theme(axis.text.y = element_text(size=20))
p_y<-p_y + theme_minimal()
#p_y<-p_y + scale_fill_manual(values=c("darkolivegreen","darkorchid4"))

plots.list<-list(p_x,p_y)
plots<-marrangeGrob(plots.list,nrow=1,ncol=1)
if(induced && testis_specific){
	ggsave(paste("fpkm_by_sexChr.protein_coding_only",dataset,cell_type,"induced_testisSpecific.pdf",sep="."), plots, width = 11, height = 8.5, units = "in")
} else if(induced){
	ggsave(paste("fpkm_by_sexChr.protein_coding_only",dataset,cell_type,"induced.pdf",sep="."), plots, width = 11, height = 8.5, units = "in")
} else if(testis_specific){
	ggsave(paste("fpkm_by_sexChr.protein_coding_only",dataset,cell_type,"testisSpecific.pdf",sep="."), plots, width = 11, height = 8.5, units = "in")
} else{
	ggsave(paste("fpkm_by_sexChr.protein_coding_only",dataset,cell_type,"pdf",sep="."), plots, width = 11, height = 8.5, units = "in")
}

#PROTEIN CODING ONLY - plot tpm values against each other for Y intro vs control (to test if overall genes are more highly expresed in one vs other or if there are some over exp and other under exp)
if(induced && testis_specific){
	pdf(paste("fpkm_by_gene.protein_coding_only",dataset,cell_type,"induced_testisSpecific.pdf", sep="."), onefile=TRUE)
} else if(induced){
	pdf(paste("fpkm_by_gene.protein_coding_only",dataset,cell_type,"induced.pdf", sep="."), onefile=TRUE)
} else if(testis_specific){
	pdf(paste("fpkm_by_gene.protein_coding_only",dataset,cell_type,"testisSpecific.pdf", sep="."), onefile=TRUE)
} else{
	pdf(paste("fpkm_by_gene.protein_coding_only",dataset,cell_type,"pdf", sep="."), onefile=TRUE)
}
plot(x=as.numeric(as.character(group_rpkm_x_proco[,1])), y=as.numeric(as.character(group_rpkm_x_proco[,2])), main="Y-intro vs control - X chromosome", xlab=colnames(group_rpkm_x_proco)[1], ylab=colnames(group_rpkm_x_proco)[2])
abline(a=0, b=1)
plot(x=as.numeric(as.character(group_rpkm_x_proco[,4])), y=as.numeric(as.character(group_rpkm_x_proco[,3])), main="Y-intro vs control - X chromosome", xlab=colnames(group_rpkm_x_proco)[4], ylab=colnames(group_rpkm_x_proco)[3])
abline(a=0, b=1)
plot(x=as.numeric(as.character(group_rpkm_y_proco[,1])), y=as.numeric(as.character(group_rpkm_y_proco[,2])), main="Y-intro vs control - Y chromosome", xlab=colnames(group_rpkm_y_proco)[1], ylab=colnames(group_rpkm_y_proco)[2])
abline(a=0, b=1)
plot(x=as.numeric(as.character(group_rpkm_y_proco[,4])), y=as.numeric(as.character(group_rpkm_y_proco[,3])), main="Y-intro vs control - Y chromosome", xlab=colnames(group_rpkm_y_proco)[4], ylab=colnames(group_rpkm_y_proco)[3])
abline(a=0, b=1)
dev.off()

print(head(which(group_rpkm_x[,1]==0)))
if(induced && testis_specific){
	write.table(group_rpkm, file=paste("group_fpkm",dataset,cell_type,"induced_testisSpecific.pdf", sep="."), sep="\t", quote=FALSE)
} else if(induced){
	write.table(group_rpkm, file=paste("group_fpkm",dataset,cell_type,"induced.pdf", sep="."), sep="\t", quote=FALSE)
} else if(testis_specific){
	write.table(group_rpkm, file=paste("group_fpkm",dataset,cell_type,"testisSpecific.pdf", sep="."), sep="\t", quote=FALSE)
} else{
	write.table(group_rpkm, file=paste("group_fpkm",dataset,cell_type,"pdf", sep="."), sep="\t", quote=FALSE)
}

##Repeat with normalization instead of log
#p_x<-ggplot(rpkm_x_df_proco, aes(x=geno, y=normalize(rpkm_value))) + geom_violin(aes(fill=mismatch)) + geom_boxplot(width=0.1)
#p_x<-p_x + labs(title="FPKM on the X chromosome by cross type - protein coding genes only", x="Cross type", y="Normalized group FPKM)")
#p_x<-p_x + theme(axis.text.y = element_text(size=20))
#p_x<-p_x + theme_minimal() #+ ylim(0,100)
#p_x<-p_x + scale_fill_manual(values=c("darkolivegreen","darkorchid4"))
#
#p_y<-ggplot(rpkm_y_df_proco, aes(x=geno, y=normalize(rpkm_value))) + geom_violin(aes(fill=mismatch)) + geom_boxplot(width=0.1)
#p_y<-p_y + labs(title="FPKM on the Y chromosome by cross type - protein coding genes only", x="Cross type", y="normalized group FPKM)")
#p_y<-p_y + theme(axis.text.y = element_text(size=20))
#p_y<-p_y + theme_minimal()
#p_y<-p_y + scale_fill_manual(values=c("darkolivegreen","darkorchid4"))
#
#plots.list<-list(p_x,p_y)
#plots<-marrangeGrob(plots.list,nrow=1,ncol=1)
#ggsave(paste("fpkm_by_sexChr.protein_coding_only.normalized",dataset,cell_type,"pdf",sep="."), plots, width = 11, height = 8.5, units = "in")


print("Done with 07_compare_sex_chromosome_expression.r")

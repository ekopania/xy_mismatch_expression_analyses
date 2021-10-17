#PURPOSE: Read in salmon and tximport output (raw expression readcounts), normalize in EdgeR, determine gene family expression levels with four methods
#1) Sum normalized counts from EdgeR object across all paralogs for each gene family
#2) Calculate FPKM for each gene, then sum FPKM values across all paralogs for each gene family
#3) Calculate FPKM treating all paralogs for a gene family as a single "gene" - sum normalized read counts across all paralogs then sum lengths of all paralogs to get total counts and lengths to use for calculating FPKM
#4) Calculated TPM treating all paralogs for a gene family as a single "gene" - same approach as method #3

#change these to test different parameters
#Normally I'd filter out low exp things but given that we want to catch everything that's part of these gene families it might be best to leave everything zero for now; one individual paralog might have a low expression level that looks like noise but it still might be contributing to the overall expression of the gene family
#Also, some of these gene families (srsx I believe) have FPKM <1 for every paralog across all samples
min_rpkm<-1 #1
min_samples<-4 #4 #16
#min_logFC<-0
#induced_cutoff<-2

#Using Larson et al. 2017 cell sorted F1s as a test at first; will eventually use the Y introgression data
dataset<-"LarsonEtal"
#cell_type<-"LZandRS"
#cell_type<-"LZ"
cell_type<-"RS"
#Percent ID threshold for considering genes paralogs/part of the same gene family
paralog_thresh<-90

print("Loading R libraries...")
library(edgeR)
library(ggplot2)
library(gridExtra)
library(EnsDb.Mmusculus.v79)
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
print(dim(these_data))
print(length(mygroups))
myDGE<-DGEList(counts=these_data, group=mygroups)
myDGE<-scaleOffset(myDGE, num_normMat)
print("Dimensions of myDGE counts:")
print(dim(myDGE$counts))
print("myDGE samples - CHECK GROUP LABELS:")
print(myDGE$samples)

print("Filtering by expression level...")
#EdgeR rpkm function seems to require same gene length across all samples; taking median of tximport output gene lengths and using that
new_lengths<-rowMedians(as.matrix(these_lengths))
names(new_lengths)<-rownames(these_lengths)
print(head(new_lengths))
my_rpkm<-rpkm(myDGE, gene.length=new_lengths)
write.table(my_rpkm,paste("fpkm_full_table",dataset,cell_type,paralog_thresh,"txt",sep="."),quote=FALSE,sep="\t")

keep<-rowSums(my_rpkm > min_rpkm) >= min_samples
keep[is.na(keep)]<-FALSE #NAs show up for gene IDs that aren't in myDGE; replace then with false otherwise myDGE gets confused

myDGE<-myDGE[keep, , keep.lib.sizes=FALSE]
myDGE<-calcNormFactors(myDGE)
print(myDGE$samples)

print("Dimensions of filtered myDGE counts:")
print(dim(myDGE$counts))
#pdf(paste("mds_plot",dataset,cell_type,"pdf",sep="."), height=8.5,width=11)
pdf(paste("mds_plot",dataset,cell_type,"withLabels.pdf",sep="."), height=8.5,width=11)
#Assign labels based on file names
labs<-sapply(colnames(myDGE$counts), function(x) gsub("\\..*_counts","",gsub("^.*/","",x)))
print(labs)
#Assign group number to group name
group_names_all<-c()
for(i in mygroups){
        if(i == 1){
		group_names_all<-c(group_names_all, "CCPP_LZ")
	} else if(i==2){
                group_names_all<-c(group_names_all, "CCPP_RS")
        } else if(i == 3){
		group_names_all<-c(group_names_all, "WWLL_LZ")
	} else if(i==4){
                group_names_all<-c(group_names_all, "WWLL_RS")
        } else if(i == 5){
		group_names_all<-c(group_names_all, "LLPP_LZ")
	} else if(i==6){
                group_names_all<-c(group_names_all, "LLPP_RS")
        } else if(i == 7){
		group_names_all<-c(group_names_all, "PPLL_LZ")
	} else if(i==8){
                group_names_all<-c(group_names_all, "PPLL_RS")
        } else{
                stop("Invalid group number")
        }
}
group_names<-unique(group_names_all)
#print(group_names)
#Assign color by cell type
cts_short<-c()
for(i in group_names){
	if(grepl("LZ", i)){
		cts_short<-c(cts_short, "chocolate1")
	} else if(grepl("RS", i)){
		cts_short<-c(cts_short, "lightsteelblue")
	} else{
		stop("Invalid cell type in group name")
	}
}
cts<-c()
for(i in group_names_all){
        if(grepl("LZ", i)){
                cts<-c(cts, "chocolate1")
        } else if(grepl("RS", i)){
                cts<-c(cts, "lightsteelblue")
        } else{
                stop("Invalid cell type in group name")
        }
}
myMDS<-plotMDS(myDGE, labels=labs, col=cts) #plot with labels
myMDS
#plotMDS(myDGE, top=nrow(myDGE$counts), pch=as.numeric(as.character(myDGE$samples$group)), col=cts, cex=3) #plot with symbols
#legend("top", legend=group_names, pch=unique(mygroups), col=cts_short, bty="n") #legend for symbols plot
dev.off()
#print(dim(myMDS$distance.matrix))
#print(myMDS$distance.matrix)
#print(dim(myMDS$cmdscale.out))
#print(myMDS$cmdscale.out)
#Save cmd scale values (what is actually plotted, seems analagous to principal component values)
#write.table(myMDS$cmdscale.out, file=paste("mds_values",dataset,"txt",sep="."), quote=FALSE, sep="\t")

#List of gene families
#gene_fams<-c("atakusan","speer","slx","slxl1","slx-slxl1","srsx","sstx","sly","srsy","ssty1","ssty2","ssty1-ssty2")
#gene_fams<-c("atakusan","speer","slx","slxl1","slx-slxl1","sly","ssty1","ssty2","ssty1-ssty2")
gene_fams<-c("atakusan","speer","astx","slx","slxl1","asty","sly","ssty1","ssty2","eif2s3x","eif2s3y")
print("Calculating gene family expression levels...")
newOrdered_lengths<-new_lengths[match(rownames(myDGE$counts),names(new_lengths))] #filtered myDGE order
new_rpkm<-rpkm(myDGE,gene.length=newOrdered_lengths) #Based on filtered gene set
print(dim(new_rpkm))
#print(new_rpkm[1:5,1:5])
write.table(new_rpkm,paste("fpkm_filtered_table",dataset,cell_type,paralog_thresh,"txt",sep="."),quote=FALSE,sep="\t")
edb<-EnsDb.Mmusculus.v79
ens_gene_table<-genes(edb, return.type="DataFrame")
head(ens_gene_table)
method1_sums<-c()
method2_sums<-c()
method3_lengths<-c()
for(g in gene_fams){
	print(paste("Working on gene family",g))
	#Get gene IDs for all paralogs in gene family ../GENE_FAMILY_FILES/atakusan.pID95.gene_family_paralogs.bed
	my_bed<-read.table(paste0("../GENE_FAMILY_FILES/",g,".pID",paralog_thresh,".gene_family_paralogs.bed"), header=FALSE, sep="\t")
	gene_names<-as.character(my_bed[,6])
	print(head(gene_names))
	gene_ids<-ens_gene_table$gene_id[which(ens_gene_table$gene_name %in% gene_names)]
	print(head(gene_ids))
	#Method 1 - sum normalized counts
	my_counts<-myDGE$counts[which(rownames(myDGE$counts) %in% gene_ids),]
	if(is.vector(my_counts)){ #Dealing with weird quirk if exactly one gene matches myDGE and gene_ids
		my_counts<-as.data.frame(rbind(my_counts, rep(0,length(my_counts))))
	}
	print(dim(my_counts))
	#if(nrow(my_counts)==1){
	#	my_sums1<-my_counts[1,]
	#} else{
	my_sums1<-colSums(my_counts)
	#}
	method1_sums<-rbind(method1_sums,my_sums1)
	#Method 2 - sum FPKMs
	my_fpkms<-new_rpkm[which(rownames(new_rpkm) %in% gene_ids),]
	if(is.vector(my_fpkms)){ #Dealing with weird quirk if exactly one gene matches myDGE and gene_ids
                my_fpkms<-as.data.frame(rbind(my_fpkms, rep(0,length(my_fpkms))))
        }
	print(dim(my_fpkms))
	if(nrow(my_fpkms)==1){
		my_sums2<-my_fpkms[1,]
	} else{
		my_sums2<-colSums(my_fpkms)
	}
	method2_sums<-rbind(method2_sums,my_sums2)
	#Method 3 - calculate FPKM based on sums and total gene length
	method3_lengths<-c(method3_lengths, sum(my_bed[,3] - my_bed[,2]))
}

print("Writing output...")
method1_sums<-as.data.frame(method1_sums)
colnames(method1_sums)<-colnames(myDGE$counts)
rownames(method1_sums)<-gene_fams
write.table(method1_sums, paste("normalized_count_sums",dataset,cell_type,paralog_thresh,"txt",sep="."), quote=FALSE, sep="\t")

method2_sums<-as.data.frame(method2_sums)
colnames(method2_sums)<-colnames(my_rpkm)
rownames(method2_sums)<-gene_fams
write.table(method2_sums, paste("fpkm_sums",dataset,cell_type,paralog_thresh,"txt",sep="."), quote=FALSE, sep="\t")

print(dim(method1_sums))
print(length(method3_lengths))
print(head(method1_sums))
print(head(method3_lengths))
method3_dge<-DGEList(counts=method1_sums)
method3_fpkms<-rpkm(method3_dge, method3_lengths)
rownames(method3_fpkms)<-gene_fams
print(head(method3_fpkms))
write.table(method3_fpkms, paste("gene_family_fpkms",dataset,cell_type,paralog_thresh,"txt",sep="."), quote=FALSE, sep="\t")

#Calculate tpm for ALL GENES and then sum tpms for genes in each gene family
divLen<-apply(myDGE$counts, 2, function(x) x/as.numeric(newOrdered_lengths)) #normalize counts by gene length
#Some checks:
#print(dim(myDGE$counts))
#print(dim(divLen))
#print(myDGE$counts[1:5,1:5])
#print(head(newOrdered_fcLengths))
#print(divLen[1:5,1:5])
scaling_factors<-colSums(divLen)/1000000
print(scaling_factors)
tpm<-apply(divLen, 1, function(x) x/scaling_factors)
tpm_transformed<-t(tpm)
print(dim(tpm_transformed))
print(head(tpm_transformed))
#Sum tpms for each gene family
gene_fam_tpms<-c()
for(g in gene_fams){
        print(paste("Working on gene family",g))
        #Get gene IDs for all paralogs in gene family
        my_bed<-read.table(paste0("../GENE_FAMILY_FILES/",g,".pID",paralog_thresh,".gene_family_paralogs.bed"), header=FALSE, sep="\t")
	gene_names<-as.character(my_bed[,6])
        print(head(gene_names))
        gene_ids<-ens_gene_table$gene_id[which(ens_gene_table$gene_name %in% gene_names)]
        print(head(gene_ids))
        #Method 4 - sum tpms
        my_tpms<-tpm_transformed[which(rownames(tpm_transformed) %in% gene_ids),]
        if(is.vector(my_tpms)){ #Dealing with weird quirk if exactly one gene matches myDGE and gene_ids
                my_tpms<-as.data.frame(rbind(my_tpms, rep(0,length(my_tpms))))
        }
	print(dim(my_tpms))
        my_sums<-colSums(my_tpms)
        gene_fam_tpms<-rbind(gene_fam_tpms,my_sums)
}
rownames(gene_fam_tpms)<-gene_fams
print(dim(gene_fam_tpms))
print(gene_fam_tpms)
write.table(tpm_transformed, paste("gene_family_tpms",dataset,cell_type,paralog_thresh,"txt",sep="."), quote=FALSE, sep="\t")

#print("Plotting tpm PCA...")
#new_tpm<-tpm[,colSums(tpm[])>0]
#my_pca<-prcomp(new_tpm,center=TRUE,scale.=TRUE) 
#summary(my_pca)
#pc1_var<-round(100*summary(my_pca)$importance[2,1], 2)
#pc2_var<-round(100*summary(my_pca)$importance[2,2], 2)
#genos<-sapply(rownames(tpm), function(x) substring(x,15,18))
#print(genos)
#mylabs<-sapply(rownames(tpm), function(x) gsub("\\.Yintro_counts","",gsub("COUNTS_Yintro\\/","",x)))
##tpm_pca<-ggplot(as.data.frame(my_pca$x),aes(x=PC1,y=PC2,shape=genos)) + geom_point(size=3)
#tpm_pca<-ggplot(as.data.frame(my_pca$x),aes(x=PC1,y=PC2)) + geom_text(label=mylabs)
#tpm_pca<-tpm_pca + labs(title=dataset, x=paste0("PC1, ", pc1_var, "% of var explaied"), y=paste0("PC2, ", pc2_var, "% of var explaied"))
#tpm_pca<-tpm_pca + theme(plot.title = element_text(size = rel(2)), axis.text=element_text(size=18), axis.title=element_text(size=21)) + theme_minimal()
#ggsave(paste("tpm_pca",dataset,paralog_thresh,"pdf",sep="."), tpm_pca, width = 11, height = 8.5, units = "in")

#Get expression levels for cell-type-specific genes from Green et al. 2018
#print("Getting cell type-specific genes...")
##edb<-EnsDb.Mmusculus.v79
#ens_mus<-useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl")
#my_genes<-getBM(attributes=c('ensembl_gene_id','external_gene_name','chromosome_name','start_position','end_position','strand','gene_biotype'), mart=ens_mus)
##edb_genes<-addFilter(edb, SeqNameFilter(c(1:19,"X","Y")))
##my_genes<-genes(edb_genes)
#
#ct_gene_table<-read.csv("../green_2018_marker_genes_byCellType.csv", header=TRUE, stringsAsFactors=FALSE)
#print(ct_gene_table)
#ct_gene_ids<-c()
#for(g in ct_gene_table$gene){
#	ct_gene_ids<-c(ct_gene_ids, my_genes$gene_id[which(my_genes$symbol==g)])
#}
#print(length(ct_gene_ids))
#ct_specific_genes<-as.data.frame(cbind(ct_gene_ids, ct_gene_table))
#print(dim(ct_specific_genes))
#print(ct_specific_genes)
#ct_specific_rpkms<-new_rpkm[which(rownames(new_rpkm) %in% ct_gene_ids),] #ENSMUSG..., use with GRCm38
##ct_specific_rpkms<-new_rpkm[which(rownames(new_rpkm) %in% ct_gene_table$gene),] #Gene names (ex: Slx), GRCm39
#print(dim(ct_specific_rpkms))
#print(ct_specific_rpkms)
#ct_specific_genes$ct_gene_ids<-as.character(ct_specific_genes$ct_gene_ids)
#print("Plotting rpkm for cell type-specific genes...")
#ordered_ctRpkms<-ct_specific_rpkms[order(rownames(ct_specific_rpkms)),]
#print(ordered_ctRpkms)
#ordered_ctSpecGenes<-ct_specific_genes[order(ct_specific_genes$ct_gene_ids),] #ENSMUSG..., use with GRCm38
##ordered_ctSpecGenes<-ct_specific_genes[order(ct_specific_genes$gene), ] #Gene names (ex: Slx), use w/ GRCm39
#print(ordered_ctSpecGenes)
##ct_values<-as.vector(ct_specific_rpkms)
#stopifnot(all.equal(rownames(ordered_ctRpkms), as.character(ordered_ctSpecGenes$ct_gene_ids))) #ENSMUSG..., GRCm38
##stopifnot(all.equal(rownames(ordered_ctRpkms), as.character(ordered_ctSpecGenes$gene))) #Gene names (ex: Slx), GRCm39
##cell_types<-rep(ct_specific_genes$CellType, ncol(ordered_ctRpkms))
##cross_types<-sapply(rownames(tpm), function(x) substring(x,15,18))
##cross_types_full<-rep(cross_types, nrow(ordered_ctRpkms))
##ct_df<-as.data.frame(cbind(ct_values, cell_types, cross_types_full))
##ct_df$ct_values<-as.numeric(as.character(ct_df$ct_values))
#ct_df<-as.data.frame(cbind(ordered_ctRpkms, CellType=as.character(ordered_ctSpecGenes$CellType), gene=as.character(ordered_ctSpecGenes$gene)))
#plots.list<-list()
#for(i in 1:(ncol(ct_df)-2)){ #loop through all RPKM value columns
#	this_samp<-gsub("\\.Yintro_counts","",gsub("COUNTS_Yintro_GRCm38\\/","",colnames(ct_df)[i]))
#	this_df<-as.data.frame(cbind(rpkm_val=as.character(ct_df[,i]), CellType=as.character(ct_df$CellType), gene=as.character(ct_df$gene)))
#	this_df$rpkm_val<-as.numeric(as.character(this_df$rpkm_val))
#	print(this_df)
#	p<-ggplot(this_df, aes(x=gene, y=rpkm_val, fill=CellType)) + geom_col()
#	p<-p + labs(title=this_samp, x="Gene", y="RPKM") + theme_minimal() + theme(axis.text.x = element_text(angle = 90))
#	p<-p + scale_x_discrete(limits=as.character(ct_gene_table$gene)) #Order by cell type
#	plots.list[[length(plots.list)+1]]<-p
#	
#}
#plots<-marrangeGrob(plots.list,nrow=2,ncol=2)
#ggsave(paste("cell_type_marker_rpkm",dataset,cell_type,"pdf",sep="."), plots, height=8.5, width=11, units="in")

print("Plotting...")
#Reformat df as three columns: gene family, type (genotype and cell type), expression value
gene_fam_vec<-c()
for(g in gene_fams){
	temp_vec<-rep(g,ncol(method1_sums))
	gene_fam_vec<-c(gene_fam_vec,temp_vec)
}
sample_vec<-rep(colnames(method1_sums),nrow(method1_sums))
print(head(gene_fam_vec))
print(head(sample_vec))
my_df<-as.data.frame(cbind(fam=gene_fam_vec,sample=sample_vec,exp=as.vector(as.matrix(t(method1_sums)))))
my_groups<-c()
my_labs_all<-c()
for(i in 1:nrow(my_df)){
	geno<-substr(sapply(my_df$sample[i], function(x) gsub("^.*/","",x)),1,4)
	ct<-gsub("\\..*","",gsub(".*M","",my_df$sample[i]))
	gf<-my_df$fam[i]
	temp_group<-paste(geno,ct,gf,sep="_")
	my_groups<-c(my_groups,temp_group)
	temp_lab<-paste(geno,ct,sep="_")
	#print(temp_group)
	#print(temp_lab)
	my_labs_all<-c(my_labs_all,temp_lab) #Making a vector of just genotype and cell type to label x-axes
}
my_labs<-unique(sort(my_labs_all))
print(my_labs)
df_counts<-as.data.frame(cbind(group=my_groups,my_df))
print(dim(df_counts))
print(head(df_counts))
df_counts_LZ<-df_counts[grepl("LZ",df_counts$group),]
df_counts_LZ_wCT<-as.data.frame(cbind(df_counts_LZ, ct=rep("LZ",nrow(df_counts_LZ))))
df_counts_RS<-df_counts[grepl("RS",df_counts$group),]
df_counts_RS_wCT<-as.data.frame(cbind(df_counts_RS, ct=rep("RS",nrow(df_counts_RS))))
df_counts_final<-as.data.frame(rbind(df_counts_LZ_wCT, df_counts_RS_wCT))
#Repeat for fpkm (instead of normalized counts) as the measure of expression level
df_fpkm<-as.data.frame(cbind(group=my_groups,fam=gene_fam_vec,sample=sample_vec,exp=as.vector(as.matrix(t(method2_sums)))))
print(dim(df_fpkm))
print(head(df_fpkm))
df_fpkm_LZ<-df_fpkm[grepl("LZ",df_fpkm$group),]
df_fpkm_LZ_wCT<-as.data.frame(cbind(df_fpkm_LZ, ct=rep("LZ",nrow(df_fpkm_LZ))))
df_fpkm_RS<-df_fpkm[grepl("RS",df_fpkm$group),]
df_fpkm_RS_wCT<-as.data.frame(cbind(df_fpkm_RS, ct=rep("RS",nrow(df_fpkm_RS))))
df_fpkm_final<-as.data.frame(rbind(df_fpkm_LZ_wCT, df_fpkm_RS_wCT))
#Repeat for gene family fpkm
df_gf_fpkm<-as.data.frame(cbind(group=my_groups,fam=gene_fam_vec,sample=sample_vec,exp=as.vector(as.matrix(t(method3_fpkms)))))
print(dim(df_gf_fpkm))
print(head(df_gf_fpkm))
df_gf_fpkm_LZ<-df_gf_fpkm[grepl("LZ",df_gf_fpkm$group),]
df_gf_fpkm_LZ_wCT<-as.data.frame(cbind(df_gf_fpkm_LZ, ct=rep("LZ",nrow(df_gf_fpkm_LZ))))
df_gf_fpkm_RS<-df_gf_fpkm[grepl("RS",df_gf_fpkm$group),]
df_gf_fpkm_RS_wCT<-as.data.frame(cbind(df_gf_fpkm_RS, ct=rep("RS",nrow(df_gf_fpkm_RS))))
df_gf_fpkm_final<-as.data.frame(rbind(df_gf_fpkm_LZ_wCT, df_gf_fpkm_RS_wCT))
#Repeat for gene family tpm
df_gf_tpm<-as.data.frame(cbind(group=my_groups,fam=gene_fam_vec,sample=sample_vec,exp=as.vector(as.matrix(t(gene_fam_tpms)))))
print(dim(df_gf_tpm))
print(head(df_gf_tpm))
df_gf_tpm_LZ<-df_gf_tpm[grepl("LZ",df_gf_tpm$group),]
df_gf_tpm_LZ_wCT<-as.data.frame(cbind(df_gf_tpm_LZ, ct=rep("LZ",nrow(df_gf_tpm_LZ))))
df_gf_tpm_RS<-df_gf_tpm[grepl("RS",df_gf_tpm$group),]
df_gf_tpm_RS_wCT<-as.data.frame(cbind(df_gf_tpm_RS, ct=rep("RS",nrow(df_gf_tpm_RS))))
df_gf_tpm_final<-as.data.frame(rbind(df_gf_tpm_LZ_wCT, df_gf_tpm_RS_wCT))
print(head(df_gf_tpm_final))

#plot by gene fam
plots.list.counts<-list()
plots.list.fpkm<-list()
plots.list.gf_fpkm<-list()
plots.list.gf_tpm<-list()
tpm_ps<-c()
comparisons_order<-c()
#geno_order<-c("CCPP_LZ","CCPP_RS","CPLY_LZ","CPLY_RS","WWLL_LZ","WWLL_RS","WLPY_LZ","WLPY_RS","LLPP_LZ","LLPP_RS","PPLL_LZ","PPLL_RS","LPLY_LZ","LPLY_RS","PLPY_LZ","PLPY_RS")
for(g in gene_fams){
	print(paste("plotting",g))
	g2<-paste0("^",g,"$") #Appends "^" and "$" to gene name to get exact matches only with grepl
	#Use this to plot LZ and RS together
	count_df<-df_counts_final[grepl(g2,df_counts_final$fam),]
	print(head(count_df))
	this_geno_order<-sapply(geno_order, function(x) paste(x, g, sep="_"))
	p_counts<-ggplot(count_df, aes(x=group, y=as.numeric(as.character(exp)), fill=ct)) + geom_boxplot(width=0.1) + labs(title=paste("Gene expression normalized counts:",g), x="Genotype", y="Normalized Expression Counts")
	p_counts<-p_counts + scale_x_discrete(limits=this_geno_order)
	p_counts<-p_counts + theme(axis.text.y = element_text(size=20))
	p_counts<-p_counts + theme_minimal() #+ scale_fill_manual(values=c("chocolate1","lightsteelblue"))
	#p_counts<-p_counts + scale_x_discrete(labels=my_labs)
	myGroup1<-this_geno_order[1] #paste(geno_order[1],g,sep="_")
        myGroup2<-this_geno_order[2] #paste(geno_order[2],g,sep="_")
        counts_result_mus<-wilcox.test(as.numeric(as.character(count_df$exp[which(count_df$group==myGroup1)])), as.numeric(as.character(count_df$exp[which(count_df$group== myGroup2)])))
        print("Wilcoxon test for musculus:")
        print(paste(myGroup1, "VS", myGroup2))
	print(counts_result_mus)
        myGroup1<-this_geno_order[3] #paste(geno_order[3],g,sep="_")
        myGroup2<-this_geno_order[4] #paste(geno_order[4],g,sep="_")
        counts_result_dom<-wilcox.test(as.numeric(as.character(count_df$exp[which(count_df$group==myGroup1)])), as.numeric(as.character(count_df$exp[which(count_df$group==myGroup2)])))
        print("Wilcoxon test for domesticus:")
        print(paste(myGroup1, "VS", myGroup2))
	print(counts_result_dom)
	#FPKMS
	fpkm_df<-df_fpkm_final[grepl(g2,df_fpkm_final$fam),]
	p_fpkm<-ggplot(fpkm_df, aes(x=group, y=as.numeric(as.character(exp)), fill=ct)) + geom_boxplot(width=0.1) + labs(title=paste("Gene expression FPKM:",g), x="Genotype", y="FPKM")
	p_fpkm<-p_fpkm + scale_x_discrete(limits=this_geno_order)
	p_fpkm<-p_fpkm + theme(axis.text.y = element_text(size=20))
	p_fpkm<-p_fpkm + theme_minimal() #+ scale_fill_manual(values=c("chocolate1","lightsteelblue"))
	#p_fpkm<-p_fpkm + scale_x_discrete(labels=my_labs)
	#GENE FAMILY FPKM
	gf_fpkm_df<-df_gf_fpkm_final[grepl(g2,df_gf_fpkm_final$fam),]
        p_gf_fpkm<-ggplot(gf_fpkm_df, aes(x=group, y=as.numeric(as.character(exp)), fill=ct)) + geom_boxplot(width=0.1) + labs(title=paste("Gene expression gene family FPKM:",g), x="Genotype", y="Gene Family FPKM")
        p_gf_fpkm<-p_gf_fpkm + scale_x_discrete(limits=this_geno_order)
	p_gf_fpkm<-p_gf_fpkm + theme(axis.text.y = element_text(size=20))
        p_gf_fpkm<-p_gf_fpkm + theme_minimal() #+ scale_fill_manual(values=c("chocolate1","lightsteelblue"))
	#p_gf_fpkm<-p_gf_fpkm + scale_x_discrete(labels=my_labs)
	#GENE FAMILY TPM
        gf_tpm_df<-df_gf_tpm_final[grepl(g2,df_gf_tpm_final$fam),]
	print(head(gf_tpm_df))
        p_gf_tpm<-ggplot(gf_tpm_df, aes(x=group, y=as.numeric(as.character(exp)))) + geom_boxplot(width=0.1) + labs(title=paste("Gene expression gene family TPM:",g), x="Genotype", y="Gene Family TPM") #,fill=ct))
        p_gf_tpm<-p_gf_tpm + scale_x_discrete(limits=this_geno_order)
	p_gf_tpm<-p_gf_tpm + theme(axis.text.y = element_text(size=20))
        p_gf_tpm<-p_gf_tpm + theme_minimal() #+ scale_fill_manual(values=c("chocolate1","lightsteelblue"))
	#p_gf_tpm<-p_gf_tpm + scale_x_discrete(labels=my_labs)
	myGroup1<-this_geno_order[1]
        myGroup2<-this_geno_order[2]
	tpm_result_mus<-wilcox.test(as.numeric(as.character(gf_tpm_df$exp[which(gf_tpm_df$group==myGroup1)])), as.numeric(as.character(gf_tpm_df$exp[which(gf_tpm_df$group== myGroup2)])))
        print("Wilcoxon test for musculus TPM:")
        print(paste(myGroup1, "VS", myGroup2))
        print(tpm_result_mus)
        myGroup1<-this_geno_order[3]
        myGroup2<-this_geno_order[4]
        tpm_result_dom<-wilcox.test(as.numeric(as.character(gf_tpm_df$exp[which(gf_tpm_df$group==myGroup1)])), as.numeric(as.character(gf_tpm_df$exp[which(gf_tpm_df$group==myGroup2)])))
        print("Wilcoxon test for domesticus TPM:")
        print(paste(myGroup1, "VS", myGroup2))
        print(tpm_result_dom)
	tpm_ps<-c(tpm_ps, tpm_result_mus$p.value, tpm_result_dom$p.value)
	comparisons_order<-c(comparisons_order, this_geno_order[1], this_geno_order[3])
	print("Pairwise Wilcoxon test (TPM):")
	tpm_result_pairwise<-pairwise.wilcox.test(as.numeric(as.character(gf_tpm_df$exp)), gf_tpm_df$group, p.adjust.method="fdr")
	print(tpm_result_pairwise)
	print("Pairwise t test (TPM)")
	tpm_result_pairwiseT<-pairwise.t.test(as.numeric(as.character(gf_tpm_df$exp)), gf_tpm_df$group, p.adjust.method="fdr", var.equal=FALSE)
        print(tpm_result_pairwiseT)
	for(i in this_geno_order){
		print(paste("Median expression levels (TPM):", i))
		print(median(as.numeric(as.character(gf_tpm_df$exp[which(gf_tpm_df$group==i)]))))
	}
	#pdf(paste("qqPlot",dataset,cell_type,paralog_thresh,"pdf",sep="."))
	#qqnorm(as.numeric(as.character(gf_tpm_df$exp)))
	#qqline(as.numeric(as.character(gf_tpm_df$exp)))
	#dev.off()
	#APPEND TO LIST
	#plots.list.counts[[length(plots.list.counts)+1]]<-pLZ_counts
	#plots.list.counts[[length(plots.list.counts)+1]]<-pRS_counts
	plots.list.counts[[length(plots.list.counts)+1]]<-p_counts
        plots.list.fpkm[[length(plots.list.fpkm)+1]]<-p_fpkm
        plots.list.gf_fpkm[[length(plots.list.gf_fpkm)+1]]<-p_gf_fpkm
        plots.list.gf_tpm[[length(plots.list.gf_tpm)+1]]<-p_gf_tpm
}

#save
plots.counts<-marrangeGrob(plots.list.counts,nrow=1,ncol=1)
plots.fpkm<-marrangeGrob(plots.list.fpkm,nrow=1,ncol=1)
plots.gf_fpkm<-marrangeGrob(plots.list.gf_fpkm,nrow=1,ncol=1)
plots.gf_tpm<-marrangeGrob(plots.list.gf_tpm,nrow=1,ncol=1)
ggsave(paste("gene_family_expression_counts",dataset,cell_type,paralog_thresh,"pdf",sep="."), plots.counts, width = 11, height = 8.5, units = "in")
ggsave(paste("gene_family_expression_fpkms",dataset,cell_type,paralog_thresh,"pdf",sep="."), plots.fpkm, width = 11, height = 8.5, units = "in")
ggsave(paste("gene_family_expression_geneFamFpkms",dataset,cell_type,paralog_thresh,"pdf",sep="."), plots.gf_fpkm, width = 11, height = 8.5, units = "in")
ggsave(paste("gene_family_expression_geneFamTpms",dataset,cell_type,paralog_thresh,"pdf",sep="."), plots.gf_tpm, width = 11, height = 8.5, units = "in")

tpm_ps_adjusted<-p.adjust(tpm_ps, method="fdr")
names(tpm_ps_adjusted)<-comparisons_order
print("Adjusted p-values for gene family TPM comparisons:")
print(tpm_ps_adjusted)

print("Done with 05_edgeR_script.r")

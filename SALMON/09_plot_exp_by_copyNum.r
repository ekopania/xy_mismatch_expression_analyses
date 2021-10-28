#PURPOSE: Plot TPM by copy number for each gene family
#1) Plot TPM by copy number for same gene (ex: Does Slx copy number correlate w/ Slx expression level?)
#2) Plot TPM by copy number of other genes (ex: Does Slx copy number correlated w/ Sly, Sst1, etc expression levels?)
#3) Plot TPM by relative copy number ratios (ex: Does Slx:Sly copy number ratio correlated w/ expression levels of each gene family?)
	#NOTE: Doing correlation on a ratio might be problematic - maybe bin these and do a Wilcox test for Fisher's exact test? Should be easy to bin b/c Y-intro lines should be way off but pure strains should be the same; I guess this is kinda the same as doing a Wilcox test between the various cross types so maybe unnecessary

library(biomaRt)
library(ggplot2)

#Read in tpm data
exp1_tpms<-read.table("gene_family_tpms.Yintro_exp1.RS.97.txt", header=TRUE)
exp2_tpms<-read.table("gene_family_tpms.Yintro_exp2.RS.97.txt", header=TRUE)
larson_tpms<-read.table("../SALMON_LarsonEtAl2017/gene_family_tpms.LarsonEtal.RS.97.txt", header=TRUE)

#Read in copy number data
cn<-read.table("/mnt/beegfs/ek112884/cnvs/AMPLICONE_ANALYSES/OWN_INFOSITE_METHOD/cnv_estimates.amplicone.combined.txt", header=FALSE)

#Get geneIDs
ens_mus<-useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl")
all_genes<-getBM(attributes=c('ensembl_gene_id','external_gene_name','chromosome_name','start_position','end_position'), mart=ens_mus)
colnames(all_genes)<-c("gene_id", "gene_name", "chr", "start", "end")
print(head(all_genes))

#Vector of gene families to include in analysis
#gene_fams<-c("atakusan","speer","astx","slx","slxl1","asty","sly","ssty1","ssty2","eif2s3x","eif2s3y")
gene_fams<-c("slx","sly")

#Estimate gene family tpms
exp1_gene_fam_tpms<-c()
exp2_gene_fam_tpms<-c()
larson_gene_fam_tpms<-c()
for(g in gene_fams){
        print(paste("Working on gene family",g))
        #Get gene IDs for all paralogs in gene family
        my_bed<-read.table(paste0("../GENE_FAMILY_FILES/",g,".pID",paralog_thresh,".gene_family_paralogs.bed"), header=FALSE, sep="\t")
        gene_names<-as.character(my_bed[,6])
        print(head(gene_names))
        gene_ids<-all_genes$gene_id[which(all_genes$gene_name %in% gene_names)]
        print(head(gene_ids))
        #Sum tpms - exp1
        my_tpms<-exp1_tpms[which(rownames(exp1_tpms) %in% gene_ids),]
        if(is.vector(my_tpms)){ #Dealing with weird quirk if exactly one gene matches myDGE and gene_ids
                my_tpms<-as.data.frame(rbind(my_tpms, rep(0,length(my_tpms))))
        }
        print(dim(my_tpms))
        my_sums<-colSums(my_tpms)
        exp1_gene_fam_tpms<-rbind(exp1_gene_fam_tpms,my_sums)
	#Sum tpms - exp2
        my_tpms<-exp2_tpms[which(rownames(exp2_tpms) %in% gene_ids),]
        if(is.vector(my_tpms)){ #Dealing with weird quirk if exactly one gene matches myDGE and gene_ids
                my_tpms<-as.data.frame(rbind(my_tpms, rep(0,length(my_tpms))))
        }
        print(dim(my_tpms))
        my_sums<-colSums(my_tpms)
	exp2_gene_fam_tpms<-rbind(exp2_gene_fam_tpms,my_sums)
	#Sum tpms - larson
        my_tpms<-larson_tpms[which(rownames(larson_tpms) %in% gene_ids),]
        if(is.vector(my_tpms)){ #Dealing with weird quirk if exactly one gene matches myDGE and gene_ids
                my_tpms<-as.data.frame(rbind(my_tpms, rep(0,length(my_tpms))))
        }
        print(dim(my_tpms))
        my_sums<-colSums(my_tpms)
	larson_gene_fam_tpms<-rbind(larson_gene_fam_tpms,my_sums)
}
rownames(exp1_gene_fam_tpms)<-gene_fams
print(dim(exp1_gene_fam_tpms))
print(exp1_gene_fam_tpms)
rownames(exp2_gene_fam_tpms)<-gene_fams
print(dim(exp2_gene_fam_tpms))
print(exp2_gene_fam_tpms)
rownames(larson_gene_fam_tpms)<-gene_fams
print(dim(larson_gene_fam_tpms))
print(larson_gene_fam_tpms)

print("Done with 09_plot_exp_by_copyNum.r")

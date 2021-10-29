#PURPOSE: Plot TPM by copy number for each gene family
#1) Plot TPM by copy number for same gene (ex: Does Slx copy number correlate w/ Slx expression level?)
#2) Plot TPM by copy number of other genes (ex: Does Slx copy number correlated w/ Sly, Sst1, etc expression levels?)
#3) Plot TPM by relative copy number ratios (ex: Does Slx:Sly copy number ratio correlated w/ expression levels of each gene family?)
	#NOTE: Doing correlation on a ratio might be problematic - maybe bin these and do a Wilcox test for Fisher's exact test? Should be easy to bin b/c Y-intro lines should be way off but pure strains should be the same; I guess this is kinda the same as doing a Wilcox test between the various cross types so maybe unnecessary

library(biomaRt)
library(ggplot2)
library(ggbeeswarm)
library(lme4)
library(lmerTest)

#Vector of gene families to include in analysis
#gene_fams<-c("atakusan","speer","astx","slx","slxl1","asty","sly","ssty1","ssty2","eif2s3x","eif2s3y")
gene_fams<-c("slx","sly","atakusan")
paralog_thresh<-97

#Read in tpm data
exp1_tpms<-read.table("gene_family_tpms.Yintro_exp1.RS.97.txt", header=TRUE)
exp2_tpms<-read.table("gene_family_tpms.Yintro_exp2.RS.97.txt", header=TRUE)
larson_tpms<-read.table("../SALMON_LarsonEtAl2017/gene_family_tpms.LarsonEtal.RS.97.txt", header=TRUE)

#Read in copy number data
cn<-read.table("/mnt/beegfs/ek112884/cnvs/AMPLICONE_ANALYSES/OWN_INFOSITE_METHOD/cnv_estimates.amplicone.combined.txt", header=FALSE, col.names=c("sample","gene_fam","pID","copy_number"))

#Get geneIDs
ens_mus<-useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl")
all_genes<-getBM(attributes=c('ensembl_gene_id','external_gene_name','chromosome_name','start_position','end_position'), mart=ens_mus)
colnames(all_genes)<-c("gene_id", "gene_name", "chr", "start", "end")
print(head(all_genes))

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
	temp_df<-as.data.frame(cbind(rep(g,length(my_sums)), rep("exp1",length(my_sums)), my_sums))
        exp1_gene_fam_tpms<-as.data.frame(rbind(exp1_gene_fam_tpms,temp_df))
	#Sum tpms - exp2
        my_tpms<-exp2_tpms[which(rownames(exp2_tpms) %in% gene_ids),]
        if(is.vector(my_tpms)){ #Dealing with weird quirk if exactly one gene matches myDGE and gene_ids
                my_tpms<-as.data.frame(rbind(my_tpms, rep(0,length(my_tpms))))
        }
        print(dim(my_tpms))
        my_sums<-colSums(my_tpms)
	temp_df<-as.data.frame(cbind(rep(g,length(my_sums)), rep("exp2",length(my_sums)), my_sums))
	exp2_gene_fam_tpms<-as.data.frame(rbind(exp2_gene_fam_tpms,temp_df))
	#Sum tpms - larson
        my_tpms<-larson_tpms[which(rownames(larson_tpms) %in% gene_ids),]
        if(is.vector(my_tpms)){ #Dealing with weird quirk if exactly one gene matches myDGE and gene_ids
                my_tpms<-as.data.frame(rbind(my_tpms, rep(0,length(my_tpms))))
        }
        print(dim(my_tpms))
        my_sums<-colSums(my_tpms)
	temp_df<-as.data.frame(cbind(rep(g, length(my_sums)), rep("larson",length(my_sums)), my_sums))
	larson_gene_fam_tpms<-as.data.frame(rbind(larson_gene_fam_tpms,temp_df))
}
colnames(exp1_gene_fam_tpms)<-c("gene_family","dataset","expression")
print(dim(exp1_gene_fam_tpms))
print(exp1_gene_fam_tpms)
colnames(exp2_gene_fam_tpms)<-c("gene_family","dataset","expression")
print(dim(exp2_gene_fam_tpms))
print(exp2_gene_fam_tpms)
colnames(larson_gene_fam_tpms)<-c("gene_family","dataset","expression")
print(dim(larson_gene_fam_tpms))
print(larson_gene_fam_tpms)

#Get copy number for each cross type
#Each one is a little different so I think it will be easier to just have a lot of code instead of writing a loop unfortunately...
#EXP1
CP_cn_slx<-cn[Reduce(intersect, list(which(cn$sample=="CCCC"), which(cn$gene_fam=="slx"), which(cn$pID==paralog_thresh))), "copy_number"]
CP_cn_sly<-cn[Reduce(intersect, list(which(cn$sample=="PPPP"), which(cn$gene_fam=="sly"), which(cn$pID==paralog_thresh))), "copy_number"]
#Amplicone reports HAPLOID copy number, so autosomal gene family copy numbers in F1s should be the mean of the two parents
CP_cn_atakusan1<-cn[Reduce(intersect, list(which(cn$sample=="CCCC"), which(cn$gene_fam=="atakusan"), which(cn$pID==paralog_thresh))), "copy_number"]
CP_cn_atakusan2<-cn[Reduce(intersect, list(which(cn$sample=="PPPP"), which(cn$gene_fam=="atakusan"), which(cn$pID==paralog_thresh))), "copy_number"]
CP_cn_atakusan<-mean(c(CP_cn_atakusan1,CP_cn_atakusan2))

CPLY_cn_slx<-cn[Reduce(intersect, list(which(cn$sample=="CCCC"), which(cn$gene_fam=="slx"), which(cn$pID==paralog_thresh))), "copy_number"]
CPLY_cn_sly<-cn[Reduce(intersect, list(which(cn$sample=="PWKLY"), which(cn$gene_fam=="sly"), which(cn$pID==paralog_thresh))), "copy_number"]
CPLY_cn_atakusan1<-cn[Reduce(intersect, list(which(cn$sample=="CCCC"), which(cn$gene_fam=="atakusan"), which(cn$pID==paralog_thresh))), "copy_number"]
CPLY_cn_atakusan2<-cn[Reduce(intersect, list(which(cn$sample=="PWKLY"), which(cn$gene_fam=="atakusan"), which(cn$pID==paralog_thresh))), "copy_number"]
CPLY_cn_atakusan<-mean(c(CPLY_cn_atakusan1,CPLY_cn_atakusan2))

WL_cn_slx<-cn[Reduce(intersect, list(which(cn$sample=="WWWW"), which(cn$gene_fam=="slx"), which(cn$pID==paralog_thresh))), "copy_number"]
WL_cn_sly<-cn[Reduce(intersect, list(which(cn$sample=="LLLL"), which(cn$gene_fam=="sly"), which(cn$pID==paralog_thresh))), "copy_number"]
WL_cn_atakusan1<-cn[Reduce(intersect, list(which(cn$sample=="WWWW"), which(cn$gene_fam=="atakusan"), which(cn$pID==paralog_thresh))), "copy_number"]
WL_cn_atakusan2<-cn[Reduce(intersect, list(which(cn$sample=="LLLL"), which(cn$gene_fam=="atakusan"), which(cn$pID==paralog_thresh))), "copy_number"]
WL_cn_atakusan<-mean(c(WL_cn_atakusan1,WL_cn_atakusan2))

WLPY_cn_slx<-cn[Reduce(intersect, list(which(cn$sample=="WWWW"), which(cn$gene_fam=="slx"), which(cn$pID==paralog_thresh))), "copy_number"]
WLPY_cn_sly<-cn[Reduce(intersect, list(which(cn$sample=="LEWPY"), which(cn$gene_fam=="sly"), which(cn$pID==paralog_thresh))), "copy_number"]
WLPY_cn_atakusan1<-cn[Reduce(intersect, list(which(cn$sample=="WWWW"), which(cn$gene_fam=="atakusan"), which(cn$pID==paralog_thresh))), "copy_number"]
WLPY_cn_atakusan2<-cn[Reduce(intersect, list(which(cn$sample=="LEWPY"), which(cn$gene_fam=="atakusan"), which(cn$pID==paralog_thresh))), "copy_number"]
WLPY_cn_atakusan<-mean(c(WLPY_cn_atakusan1,WLPY_cn_atakusan2))

#EXP2
LP_cn_slx<-cn[Reduce(intersect, list(which(cn$sample=="LLLL"), which(cn$gene_fam=="slx"), which(cn$pID==paralog_thresh))), "copy_number"]
LP_cn_sly<-cn[Reduce(intersect, list(which(cn$sample=="PPPP"), which(cn$gene_fam=="sly"), which(cn$pID==paralog_thresh))), "copy_number"]
LP_cn_atakusan1<-cn[Reduce(intersect, list(which(cn$sample=="LLLL"), which(cn$gene_fam=="atakusan"), which(cn$pID==paralog_thresh))), "copy_number"]
LP_cn_atakusan2<-cn[Reduce(intersect, list(which(cn$sample=="PPPP"), which(cn$gene_fam=="atakusan"), which(cn$pID==paralog_thresh))), "copy_number"]
LP_cn_atakusan<-mean(c(LP_cn_atakusan1,LP_cn_atakusan2))

LPLY_cn_slx<-cn[Reduce(intersect, list(which(cn$sample=="LLLL"), which(cn$gene_fam=="slx"), which(cn$pID==paralog_thresh))), "copy_number"]
LPLY_cn_sly<-cn[Reduce(intersect, list(which(cn$sample=="PWKLY"), which(cn$gene_fam=="sly"), which(cn$pID==paralog_thresh))), "copy_number"]
LPLY_cn_atakusan1<-cn[Reduce(intersect, list(which(cn$sample=="LLLL"), which(cn$gene_fam=="atakusan"), which(cn$pID==paralog_thresh))), "copy_number"]
LPLY_cn_atakusan2<-cn[Reduce(intersect, list(which(cn$sample=="PWKLY"), which(cn$gene_fam=="atakusan"), which(cn$pID==paralog_thresh))), "copy_number"]
LPLY_cn_atakusan<-mean(c(LPLY_cn_atakusan1,LPLY_cn_atakusan2))

PL_cn_slx<-cn[Reduce(intersect, list(which(cn$sample=="PPPP"), which(cn$gene_fam=="slx"), which(cn$pID==paralog_thresh))), "copy_number"]
PL_cn_sly<-cn[Reduce(intersect, list(which(cn$sample=="LLLL"), which(cn$gene_fam=="sly"), which(cn$pID==paralog_thresh))), "copy_number"]
PL_cn_atakusan1<-cn[Reduce(intersect, list(which(cn$sample=="PPPP"), which(cn$gene_fam=="atakusan"), which(cn$pID==paralog_thresh))), "copy_number"]
PL_cn_atakusan2<-cn[Reduce(intersect, list(which(cn$sample=="LLLL"), which(cn$gene_fam=="atakusan"), which(cn$pID==paralog_thresh))), "copy_number"]
PL_cn_atakusan<-mean(c(PL_cn_atakusan1,PL_cn_atakusan2))

PLPY_cn_slx<-cn[Reduce(intersect, list(which(cn$sample=="PPPP"), which(cn$gene_fam=="slx"), which(cn$pID==paralog_thresh))), "copy_number"]
PLPY_cn_sly<-cn[Reduce(intersect, list(which(cn$sample=="LEWPY"), which(cn$gene_fam=="sly"), which(cn$pID==paralog_thresh))), "copy_number"]
PLPY_cn_atakusan1<-as.numeric(as.character(cn[Reduce(intersect, list(which(cn$sample=="PPPP"), which(cn$gene_fam=="atakusan"), which(cn$pID==paralog_thresh))), "copy_number"]))
PLPY_cn_atakusan2<-as.numeric(as.character(cn[Reduce(intersect, list(which(cn$sample=="LEWPY"), which(cn$gene_fam=="atakusan"), which(cn$pID==paralog_thresh))), "copy_number"]))
PLPY_cn_atakusan<-mean(c(PLPY_cn_atakusan1,PLPY_cn_atakusan2))

#Make copy number columns
copy_number_exp1<-c(rep(CP_cn_slx,4), rep(CPLY_cn_slx,4), rep(WL_cn_slx,4), rep(WLPY_cn_slx,4), rep(CP_cn_sly,4), rep(CPLY_cn_sly,4), rep(WL_cn_sly,4), rep(WLPY_cn_sly,4), rep(CP_cn_atakusan,4), rep(CPLY_cn_atakusan,4), rep(WL_cn_atakusan,4), rep(WLPY_cn_atakusan,4))
copy_number_exp2<-c(rep(LP_cn_slx,3), rep(LPLY_cn_slx,4), rep(PL_cn_slx,4), rep(PLPY_cn_slx,4), rep(LP_cn_sly,3), rep(LPLY_cn_sly,4), rep(PL_cn_sly,4), rep(PLPY_cn_sly,4), rep(LP_cn_atakusan,3), rep(LPLY_cn_atakusan,4), rep(PL_cn_atakusan,4), rep(PLPY_cn_atakusan,4))
copy_number_larson<-c(rep(CP_cn_slx,3), rep(LP_cn_slx,3), rep(PL_cn_slx,3),  rep(WL_cn_slx,3), rep(CP_cn_sly,3), rep(LP_cn_sly,3), rep(PL_cn_sly,3), rep(WL_cn_sly,3), rep(CP_cn_atakusan,3), rep(LP_cn_atakusan,3), rep(PL_cn_atakusan,3), rep(WL_cn_atakusan,3))

#Make cross type columns
cross_type_exp1<-sapply(rownames(exp1_gene_fam_tpms), function(x) substr(x,1,4))
cross_type_exp2<-sapply(rownames(exp2_gene_fam_tpms), function(x) substr(x,1,4))
cross_type_larson<-sapply(rownames(larson_gene_fam_tpms), function(x) substr(x,1,4))

#Make group columns (cross_dataset)
group_exp1<-sapply(cross_type_exp1, function(x) paste0(x,"_exp1"))
group_exp2<-sapply(cross_type_exp2, function(x) paste0(x,"_exp2"))
group_larson<-sapply(cross_type_larson, function(x) paste0(x,"_larson"))

#Add columns of things to test in linear models (XY mismatch vs not, sterile vs fertile, F1 vs non-F1 hybrid background)
mismatch_exp1<-sapply(cross_type_exp1, function(x) if(x %in% c("CPLY","WLPY","PPLL","LLPP")) "yes" else "no")
mismatch_exp2<-sapply(cross_type_exp2, function(x) if(x %in% c("CPLY","WLPY","PPLL","LLPP")) "yes" else "no")
mismatch_larson<-sapply(cross_type_larson, function(x) if(x %in% c("CPLY","WLPY","PPLL","LLPP")) "yes" else "no")
sterile_exp1<-sapply(cross_type_exp1, function(x) if(x %in% c("PLPY","PPLL")) "yes" else "no")
sterile_exp2<-sapply(cross_type_exp2, function(x) if(x %in% c("PLPY","PPLL")) "yes" else "no")
sterile_larson<-sapply(cross_type_larson, function(x) if(x %in% c("PLPY","PPLL")) "yes" else "no")
F1hybBG_exp1<-sapply(cross_type_exp1, function(x) if(x %in% c("LPLY","PLPY","PPLL","LLPP")) "yes" else "no")
F1hybBG_exp2<-sapply(cross_type_exp2, function(x) if(x %in% c("LPLY","PLPY","PPLL","LLPP")) "yes" else "no")
F1hybBG_larson<-sapply(cross_type_larson, function(x) if(x %in% c("LPLY","PLPY","PPLL","LLPP")) "yes" else "no")

exp1_df<-as.data.frame(cbind(exp1_gene_fam_tpms, copy_number=copy_number_exp1, cross_type=cross_type_exp1, group=group_exp1, mismatch=mismatch_exp1, sterile=sterile_exp1, F1hybBG=F1hybBG_exp1))
exp2_df<-as.data.frame(cbind(exp2_gene_fam_tpms, copy_number=copy_number_exp2, cross_type=cross_type_exp2, group=group_exp2, mismatch=mismatch_exp2, sterile=sterile_exp2, F1hybBG=F1hybBG_exp2))
larson_df<-as.data.frame(cbind(larson_gene_fam_tpms, copy_number=copy_number_larson, cross_type=cross_type_larson, group=group_larson, mismatch=mismatch_larson, sterile=sterile_larson, F1hybBG=F1hybBG_larson))
print(exp1_df)
print(exp2_df)
print(larson_df)


mydf<-as.data.frame(rbind(exp1_df,exp2_df,larson_df))
print(mydf)

#Plot
p_slx<-ggplot(mydf[which(mydf$gene_family=="slx"),], aes(x=as.numeric(as.character(copy_number)), y=as.numeric(as.character(expression)), group=group, shape=dataset, color=cross_type)) + geom_boxplot(width=3) + geom_quasirandom() #+ geom_text(vjust=0, nudge_y=0.1)
p_slx<-p_slx + labs(title="Expression vs Copy Number - Slx", x="Copy Number", y="Median Gene Family TPM")
p_slx<-p_slx + theme(axis.text=element_text(size=18), axis.title=element_text(size=21), plot.title=element_text(size=28))
p_slx<-p_slx + theme_minimal()

p_sly<-ggplot(mydf[which(mydf$gene_family=="sly"),], aes(x=as.numeric(as.character(copy_number)), y=as.numeric(as.character(expression)), group=group, shape=dataset, color=cross_type)) + geom_boxplot(width=10) + geom_quasirandom() #+ geom_text(vjust=0, nudge_y=0.1)
p_sly<-p_sly + labs(title="Expression vs Copy Number - Sly", x="Copy Number", y="Median Gene Family TPM")
p_sly<-p_sly + theme(axis.text=element_text(size=18), axis.title=element_text(size=21), plot.title=element_text(size=28))
p_sly<-p_sly + theme_minimal()

p_a<-ggplot(mydf[which(mydf$gene_family=="atakusan"),], aes(x=as.numeric(as.character(copy_number)), y=as.numeric(as.character(expression)), group=group, shape=dataset, color=cross_type)) + geom_boxplot(width=25) + geom_quasirandom() #+ geom_text(vjust=0, nudge_y=0.1)
p_a<-p_a + labs(title="Expression vs Copy Number - Alpha-takusan", x="Copy Number", y="Median Gene Family TPM")
p_a<-p_a + theme(axis.text=element_text(size=18), axis.title=element_text(size=21), plot.title=element_text(size=28))
p_a<-p_a + theme_minimal()

pdf("geneFam_expression_by_copyNumber.RS.pdf")
print(p_slx)
print(p_sly)
print(p_a)
dev.off()

#Linear models
print("Linear models: effect of XY mismatch\n")
print("lmer test for Slx:")
null_model<-lmer(as.numeric(as.character(expression)) ~ as.numeric(as.character(copy_number)) + (1|dataset), data=mydf[which(mydf$gene_family=="slx"),])
with_XYmismatch<-lmer(as.numeric(as.character(expression)) ~ as.numeric(as.character(copy_number)) + mismatch + (1|dataset), data=mydf[which(mydf$gene_family=="slx"),])
print(anova(with_XYmismatch, null_model))
print(summary(with_XYmismatch))
print("lmer test for Sly:")
null_model<-lmer(as.numeric(as.character(expression)) ~ as.numeric(as.character(copy_number)) + (1|dataset), data=mydf[which(mydf$gene_family=="sly"),])
with_XYmismatch<-lmer(as.numeric(as.character(expression)) ~ as.numeric(as.character(copy_number)) + mismatch + (1|dataset), data=mydf[which(mydf$gene_family=="sly"),])
print(anova(with_XYmismatch, null_model))
print(summary(with_XYmismatch))
print("lmer test for Alpha-takusan:")
null_model<-lmer(as.numeric(as.character(expression)) ~ as.numeric(as.character(copy_number)) + (1|dataset), data=mydf[which(mydf$gene_family=="atakusan"),])
with_XYmismatch<-lmer(as.numeric(as.character(expression)) ~ as.numeric(as.character(copy_number)) + mismatch + (1|dataset), data=mydf[which(mydf$gene_family=="atakusan"),])
print(anova(with_XYmismatch, null_model))
print(summary(with_XYmismatch))


print("\nLinear models: effect of sterile background (PPLL and PPLLPY)\n")
print("lmer test for Slx:")
null_model<-lmer(as.numeric(as.character(expression)) ~ as.numeric(as.character(copy_number)) + (1|dataset), data=mydf[which(mydf$gene_family=="slx"),])
with_sterile<-lmer(as.numeric(as.character(expression)) ~ as.numeric(as.character(copy_number)) + sterile + (1|dataset), data=mydf[which(mydf$gene_family=="slx"),])
print(anova(with_sterile, null_model))
print(summary(with_sterile))
print("lmer test for Sly:")
null_model<-lmer(as.numeric(as.character(expression)) ~ as.numeric(as.character(copy_number)) + (1|dataset), data=mydf[which(mydf$gene_family=="sly"),])
with_sterile<-lmer(as.numeric(as.character(expression)) ~ as.numeric(as.character(copy_number)) + sterile + (1|dataset), data=mydf[which(mydf$gene_family=="sly"),])
print(anova(with_sterile, null_model))
print(summary(with_sterile))
print("lmer test for Alpha-takusan:")
null_model<-lmer(as.numeric(as.character(expression)) ~ as.numeric(as.character(copy_number)) + (1|dataset), data=mydf[which(mydf$gene_family=="atakusan"),])
with_sterile<-lmer(as.numeric(as.character(expression)) ~ as.numeric(as.character(copy_number)) + sterile + (1|dataset), data=mydf[which(mydf$gene_family=="atakusan"),])
print(anova(with_sterile, null_model))
print(summary(with_sterile))

print("\nLinear models: effect of F1 hybrid background\n")
print("lmer test for Slx:")
null_model<-lmer(as.numeric(as.character(expression)) ~ as.numeric(as.character(copy_number)) + (1|dataset), data=mydf[which(mydf$gene_family=="slx"),])
with_F1hybBG<-lmer(as.numeric(as.character(expression)) ~ as.numeric(as.character(copy_number)) + F1hybBG + (1|dataset), data=mydf[which(mydf$gene_family=="slx"),])
print(anova(with_F1hybBG, null_model))
print(summary(with_F1hybBG))
print("lmer test for Sly:")
null_model<-lmer(as.numeric(as.character(expression)) ~ as.numeric(as.character(copy_number)) + (1|dataset), data=mydf[which(mydf$gene_family=="sly"),])
with_F1hybBG<-lmer(as.numeric(as.character(expression)) ~ as.numeric(as.character(copy_number)) + F1hybBG + (1|dataset), data=mydf[which(mydf$gene_family=="sly"),])
print(anova(with_F1hybBG, null_model))
print(summary(with_F1hybBG))
print("lmer test for Alpha-takusan:")
null_model<-lmer(as.numeric(as.character(expression)) ~ as.numeric(as.character(copy_number)) + (1|dataset), data=mydf[which(mydf$gene_family=="atakusan"),])
with_F1hybBG<-lmer(as.numeric(as.character(expression)) ~ as.numeric(as.character(copy_number)) + F1hybBG + (1|dataset), data=mydf[which(mydf$gene_family=="atakusan"),])
print(anova(with_F1hybBG, null_model))
print(summary(with_F1hybBG))

print("Done with 09_plot_exp_by_copyNum.r")

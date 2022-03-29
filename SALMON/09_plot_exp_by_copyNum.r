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
gene_fams<-c("slx","slxl1","sly","ssty1","ssty2","atakusan")
#gene_fams<-c("speer")
cell_type<-"RS"
paralog_thresh<-97

#Read in tpm data
exp1_tpms<-read.table(paste("gene_family_tpms.Yintro_exp1",cell_type,paralog_thresh,"txt", sep="."), header=TRUE)
exp2_tpms<-read.table(paste("gene_family_tpms.Yintro_exp2",cell_type,paralog_thresh,"txt", sep="."), header=TRUE)
#larson_tpms<-read.table(paste("../SALMON_LarsonEtAl2017/gene_family_tpms.LarsonEtal",cell_type,paralog_thresh,"txt", sep="."), header=TRUE)

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
#larson_gene_fam_tpms<-c()
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
#        my_tpms<-larson_tpms[which(rownames(larson_tpms) %in% gene_ids),]
#        if(is.vector(my_tpms)){ #Dealing with weird quirk if exactly one gene matches myDGE and gene_ids
#                my_tpms<-as.data.frame(rbind(my_tpms, rep(0,length(my_tpms))))
#        }
#        print(dim(my_tpms))
#        my_sums<-colSums(my_tpms)
#	temp_df<-as.data.frame(cbind(rep(g, length(my_sums)), rep("larson",length(my_sums)), my_sums))
#	larson_gene_fam_tpms<-as.data.frame(rbind(larson_gene_fam_tpms,temp_df))
}
colnames(exp1_gene_fam_tpms)<-c("gene_family","dataset","expression")
print(dim(exp1_gene_fam_tpms))
print(exp1_gene_fam_tpms)
colnames(exp2_gene_fam_tpms)<-c("gene_family","dataset","expression")
print(dim(exp2_gene_fam_tpms))
print(exp2_gene_fam_tpms)
#colnames(larson_gene_fam_tpms)<-c("gene_family","dataset","expression")
#print(dim(larson_gene_fam_tpms))
#print(larson_gene_fam_tpms)

#Get copy number for each cross type
#Each one is a little different so I think it will be easier to just have a lot of code instead of writing a loop unfortunately...
#EXP1
CP_cn_slx<-cn[Reduce(intersect, list(which(cn$sample=="CCCC"), which(cn$gene_fam=="slx"), which(cn$pID==paralog_thresh))), "copy_number"]
CP_cn_slxl1<-cn[Reduce(intersect, list(which(cn$sample=="CCCC"), which(cn$gene_fam=="slxl1"), which(cn$pID==paralog_thresh))), "copy_number"]
CP_cn_sly<-cn[Reduce(intersect, list(which(cn$sample=="PPPP"), which(cn$gene_fam=="sly"), which(cn$pID==paralog_thresh))), "copy_number"]
CP_cn_ssty1<-cn[Reduce(intersect, list(which(cn$sample=="PPPP"), which(cn$gene_fam=="ssty1"), which(cn$pID==paralog_thresh))), "copy_number"]
CP_cn_ssty2<-cn[Reduce(intersect, list(which(cn$sample=="PPPP"), which(cn$gene_fam=="ssty2"), which(cn$pID==paralog_thresh))), "copy_number"]
#Amplicone reports HAPLOID copy number, so autosomal gene family copy numbers in F1s should be the mean of the two parents
CP_cn_atakusan1<-cn[Reduce(intersect, list(which(cn$sample=="CCCC"), which(cn$gene_fam=="atakusan"), which(cn$pID==paralog_thresh))), "copy_number"]
CP_cn_atakusan2<-cn[Reduce(intersect, list(which(cn$sample=="PPPP"), which(cn$gene_fam=="atakusan"), which(cn$pID==paralog_thresh))), "copy_number"]
CP_cn_atakusan<-mean(c(CP_cn_atakusan1,CP_cn_atakusan2))

CPLY_cn_slx<-cn[Reduce(intersect, list(which(cn$sample=="CCCC"), which(cn$gene_fam=="slx"), which(cn$pID==paralog_thresh))), "copy_number"]
CPLY_cn_slxl1<-cn[Reduce(intersect, list(which(cn$sample=="CCCC"), which(cn$gene_fam=="slxl1"), which(cn$pID==paralog_thresh))), "copy_number"]
CPLY_cn_sly<-cn[Reduce(intersect, list(which(cn$sample=="PWKLY"), which(cn$gene_fam=="sly"), which(cn$pID==paralog_thresh))), "copy_number"]
CPLY_cn_ssty1<-cn[Reduce(intersect, list(which(cn$sample=="PWKLY"), which(cn$gene_fam=="ssty1"), which(cn$pID==paralog_thresh))), "copy_number"]
CPLY_cn_ssty2<-cn[Reduce(intersect, list(which(cn$sample=="PWKLY"), which(cn$gene_fam=="ssty2"), which(cn$pID==paralog_thresh))), "copy_number"]
CPLY_cn_atakusan1<-cn[Reduce(intersect, list(which(cn$sample=="CCCC"), which(cn$gene_fam=="atakusan"), which(cn$pID==paralog_thresh))), "copy_number"]
CPLY_cn_atakusan2<-cn[Reduce(intersect, list(which(cn$sample=="PWKLY"), which(cn$gene_fam=="atakusan"), which(cn$pID==paralog_thresh))), "copy_number"]
CPLY_cn_atakusan<-mean(c(CPLY_cn_atakusan1,CPLY_cn_atakusan2))

WL_cn_slx<-cn[Reduce(intersect, list(which(cn$sample=="WWWW"), which(cn$gene_fam=="slx"), which(cn$pID==paralog_thresh))), "copy_number"]
WL_cn_slxl1<-cn[Reduce(intersect, list(which(cn$sample=="WWWW"), which(cn$gene_fam=="slxl1"), which(cn$pID==paralog_thresh))), "copy_number"]
WL_cn_sly<-cn[Reduce(intersect, list(which(cn$sample=="LLLL"), which(cn$gene_fam=="sly"), which(cn$pID==paralog_thresh))), "copy_number"]
WL_cn_ssty1<-cn[Reduce(intersect, list(which(cn$sample=="LLLL"), which(cn$gene_fam=="ssty1"), which(cn$pID==paralog_thresh))), "copy_number"]
WL_cn_ssty2<-cn[Reduce(intersect, list(which(cn$sample=="LLLL"), which(cn$gene_fam=="ssty2"), which(cn$pID==paralog_thresh))), "copy_number"]
WL_cn_atakusan1<-cn[Reduce(intersect, list(which(cn$sample=="WWWW"), which(cn$gene_fam=="atakusan"), which(cn$pID==paralog_thresh))), "copy_number"]
WL_cn_atakusan2<-cn[Reduce(intersect, list(which(cn$sample=="LLLL"), which(cn$gene_fam=="atakusan"), which(cn$pID==paralog_thresh))), "copy_number"]
WL_cn_atakusan<-mean(c(WL_cn_atakusan1,WL_cn_atakusan2))

WLPY_cn_slx<-cn[Reduce(intersect, list(which(cn$sample=="WWWW"), which(cn$gene_fam=="slx"), which(cn$pID==paralog_thresh))), "copy_number"]
WLPY_cn_slxl1<-cn[Reduce(intersect, list(which(cn$sample=="WWWW"), which(cn$gene_fam=="slxl1"), which(cn$pID==paralog_thresh))), "copy_number"]
WLPY_cn_sly<-cn[Reduce(intersect, list(which(cn$sample=="LEWPY"), which(cn$gene_fam=="sly"), which(cn$pID==paralog_thresh))), "copy_number"]
WLPY_cn_ssty1<-cn[Reduce(intersect, list(which(cn$sample=="LEWPY"), which(cn$gene_fam=="ssty1"), which(cn$pID==paralog_thresh))), "copy_number"]
WLPY_cn_ssty2<-cn[Reduce(intersect, list(which(cn$sample=="LEWPY"), which(cn$gene_fam=="ssty2"), which(cn$pID==paralog_thresh))), "copy_number"]
WLPY_cn_atakusan1<-cn[Reduce(intersect, list(which(cn$sample=="WWWW"), which(cn$gene_fam=="atakusan"), which(cn$pID==paralog_thresh))), "copy_number"]
WLPY_cn_atakusan2<-cn[Reduce(intersect, list(which(cn$sample=="LEWPY"), which(cn$gene_fam=="atakusan"), which(cn$pID==paralog_thresh))), "copy_number"]
WLPY_cn_atakusan<-mean(c(WLPY_cn_atakusan1,WLPY_cn_atakusan2))

#EXP2
LP_cn_slx<-cn[Reduce(intersect, list(which(cn$sample=="LLLL"), which(cn$gene_fam=="slx"), which(cn$pID==paralog_thresh))), "copy_number"]
LP_cn_slxl1<-cn[Reduce(intersect, list(which(cn$sample=="LLLL"), which(cn$gene_fam=="slxl1"), which(cn$pID==paralog_thresh))), "copy_number"]
LP_cn_sly<-cn[Reduce(intersect, list(which(cn$sample=="PPPP"), which(cn$gene_fam=="sly"), which(cn$pID==paralog_thresh))), "copy_number"]
LP_cn_ssty1<-cn[Reduce(intersect, list(which(cn$sample=="PPPP"), which(cn$gene_fam=="ssty1"), which(cn$pID==paralog_thresh))), "copy_number"]
LP_cn_ssty2<-cn[Reduce(intersect, list(which(cn$sample=="PPPP"), which(cn$gene_fam=="ssty2"), which(cn$pID==paralog_thresh))), "copy_number"]
LP_cn_atakusan1<-cn[Reduce(intersect, list(which(cn$sample=="LLLL"), which(cn$gene_fam=="atakusan"), which(cn$pID==paralog_thresh))), "copy_number"]
LP_cn_atakusan2<-cn[Reduce(intersect, list(which(cn$sample=="PPPP"), which(cn$gene_fam=="atakusan"), which(cn$pID==paralog_thresh))), "copy_number"]
LP_cn_atakusan<-mean(c(LP_cn_atakusan1,LP_cn_atakusan2))

LPLY_cn_slx<-cn[Reduce(intersect, list(which(cn$sample=="LLLL"), which(cn$gene_fam=="slx"), which(cn$pID==paralog_thresh))), "copy_number"]
LPLY_cn_slxl1<-cn[Reduce(intersect, list(which(cn$sample=="LLLL"), which(cn$gene_fam=="slxl1"), which(cn$pID==paralog_thresh))), "copy_number"]
LPLY_cn_sly<-cn[Reduce(intersect, list(which(cn$sample=="PWKLY"), which(cn$gene_fam=="sly"), which(cn$pID==paralog_thresh))), "copy_number"]
LPLY_cn_ssty1<-cn[Reduce(intersect, list(which(cn$sample=="PWKLY"), which(cn$gene_fam=="ssty1"), which(cn$pID==paralog_thresh))), "copy_number"]
LPLY_cn_ssty2<-cn[Reduce(intersect, list(which(cn$sample=="PWKLY"), which(cn$gene_fam=="ssty2"), which(cn$pID==paralog_thresh))), "copy_number"]
LPLY_cn_atakusan1<-cn[Reduce(intersect, list(which(cn$sample=="LLLL"), which(cn$gene_fam=="atakusan"), which(cn$pID==paralog_thresh))), "copy_number"]
LPLY_cn_atakusan2<-cn[Reduce(intersect, list(which(cn$sample=="PWKLY"), which(cn$gene_fam=="atakusan"), which(cn$pID==paralog_thresh))), "copy_number"]
LPLY_cn_atakusan<-mean(c(LPLY_cn_atakusan1,LPLY_cn_atakusan2))

PL_cn_slx<-cn[Reduce(intersect, list(which(cn$sample=="PPPP"), which(cn$gene_fam=="slx"), which(cn$pID==paralog_thresh))), "copy_number"]
PL_cn_slxl1<-cn[Reduce(intersect, list(which(cn$sample=="PPPP"), which(cn$gene_fam=="slxl1"), which(cn$pID==paralog_thresh))), "copy_number"]
PL_cn_sly<-cn[Reduce(intersect, list(which(cn$sample=="LLLL"), which(cn$gene_fam=="sly"), which(cn$pID==paralog_thresh))), "copy_number"]
PL_cn_ssty1<-cn[Reduce(intersect, list(which(cn$sample=="LLLL"), which(cn$gene_fam=="ssty1"), which(cn$pID==paralog_thresh))), "copy_number"]
PL_cn_ssty2<-cn[Reduce(intersect, list(which(cn$sample=="LLLL"), which(cn$gene_fam=="ssty2"), which(cn$pID==paralog_thresh))), "copy_number"]
PL_cn_atakusan1<-cn[Reduce(intersect, list(which(cn$sample=="PPPP"), which(cn$gene_fam=="atakusan"), which(cn$pID==paralog_thresh))), "copy_number"]
PL_cn_atakusan2<-cn[Reduce(intersect, list(which(cn$sample=="LLLL"), which(cn$gene_fam=="atakusan"), which(cn$pID==paralog_thresh))), "copy_number"]
PL_cn_atakusan<-mean(c(PL_cn_atakusan1,PL_cn_atakusan2))

PLPY_cn_slx<-cn[Reduce(intersect, list(which(cn$sample=="PPPP"), which(cn$gene_fam=="slx"), which(cn$pID==paralog_thresh))), "copy_number"]
PLPY_cn_slxl1<-cn[Reduce(intersect, list(which(cn$sample=="PPPP"), which(cn$gene_fam=="slxl1"), which(cn$pID==paralog_thresh))), "copy_number"]
PLPY_cn_sly<-cn[Reduce(intersect, list(which(cn$sample=="LEWPY"), which(cn$gene_fam=="sly"), which(cn$pID==paralog_thresh))), "copy_number"]
PLPY_cn_ssty1<-cn[Reduce(intersect, list(which(cn$sample=="LEWPY"), which(cn$gene_fam=="ssty1"), which(cn$pID==paralog_thresh))), "copy_number"]
PLPY_cn_ssty2<-cn[Reduce(intersect, list(which(cn$sample=="LEWPY"), which(cn$gene_fam=="ssty2"), which(cn$pID==paralog_thresh))), "copy_number"]
PLPY_cn_atakusan1<-as.numeric(as.character(cn[Reduce(intersect, list(which(cn$sample=="PPPP"), which(cn$gene_fam=="atakusan"), which(cn$pID==paralog_thresh))), "copy_number"]))
PLPY_cn_atakusan2<-as.numeric(as.character(cn[Reduce(intersect, list(which(cn$sample=="LEWPY"), which(cn$gene_fam=="atakusan"), which(cn$pID==paralog_thresh))), "copy_number"]))
PLPY_cn_atakusan<-mean(c(PLPY_cn_atakusan1,PLPY_cn_atakusan2))

#Make copy number columns
copy_number_exp1<-c(rep(CP_cn_slx,4), rep(CPLY_cn_slx,4), rep(WL_cn_slx,4), rep(WLPY_cn_slx,4), rep(CP_cn_slxl1,4), rep(CPLY_cn_slxl1,4), rep(WL_cn_slxl1,4), rep(WLPY_cn_slxl1,4), rep(CP_cn_sly,4), rep(CPLY_cn_sly,4), rep(WL_cn_sly,4), rep(WLPY_cn_sly,4), rep(CP_cn_ssty1,4), rep(CPLY_cn_ssty1,4), rep(WL_cn_ssty1,4), rep(WLPY_cn_ssty1,4), rep(CP_cn_ssty2,4), rep(CPLY_cn_ssty2,4), rep(WL_cn_ssty2,4), rep(WLPY_cn_ssty2,4), rep(CP_cn_atakusan,4), rep(CPLY_cn_atakusan,4), rep(WL_cn_atakusan,4), rep(WLPY_cn_atakusan,4))
if(cell_type=="LZ"){
	copy_number_exp2<-c(rep(LP_cn_slx,4), rep(LPLY_cn_slx,4), rep(PL_cn_slx,4), rep(PLPY_cn_slx,4), rep(LP_cn_slxl1,4), rep(LPLY_cn_slxl1,4), rep(PL_cn_slxl1,4), rep(PLPY_cn_slxl1,4), rep(LP_cn_sly,4), rep(LPLY_cn_sly,4), rep(PL_cn_sly,4), rep(PLPY_cn_sly,4), rep(LP_cn_ssty1,4), rep(LPLY_cn_ssty1,4), rep(PL_cn_ssty1,4), rep(PLPY_cn_ssty1,4), rep(LP_cn_ssty2,4), rep(LPLY_cn_ssty2,4), rep(PL_cn_ssty2,4), rep(PLPY_cn_ssty2,4), rep(LP_cn_atakusan,4), rep(LPLY_cn_atakusan,4), rep(PL_cn_atakusan,4), rep(PLPY_cn_atakusan,4))
} else if(cell_type=="RS"){
	copy_number_exp2<-c(rep(LP_cn_slx,3), rep(LPLY_cn_slx,4), rep(PL_cn_slx,4), rep(PLPY_cn_slx,4), rep(LP_cn_slxl1,3), rep(LPLY_cn_slxl1,4), rep(PL_cn_slxl1,4), rep(PLPY_cn_slxl1,4), rep(LP_cn_sly,3), rep(LPLY_cn_sly,4), rep(PL_cn_sly,4), rep(PLPY_cn_sly,4), rep(LP_cn_ssty1,3), rep(LPLY_cn_ssty1,4), rep(PL_cn_ssty1,4), rep(PLPY_cn_ssty1,4), rep(LP_cn_ssty2,3), rep(LPLY_cn_ssty2,4), rep(PL_cn_ssty2,4), rep(PLPY_cn_ssty2,4), rep(LP_cn_atakusan,3), rep(LPLY_cn_atakusan,4), rep(PL_cn_atakusan,4), rep(PLPY_cn_atakusan,4))
}
#copy_number_larson<-c(rep(CP_cn_slx,3), rep(LP_cn_slx,3), rep(PL_cn_slx,3),  rep(WL_cn_slx,3), rep(CP_cn_sly,3), rep(LP_cn_sly,3), rep(PL_cn_sly,3), rep(WL_cn_sly,3), rep(CP_cn_atakusan,3), rep(LP_cn_atakusan,3), rep(PL_cn_atakusan,3), rep(WL_cn_atakusan,3))

#Make cross type columns
cross_type_exp1<-sapply(rownames(exp1_gene_fam_tpms), function(x) substr(x,1,4))
cross_type_exp2<-sapply(rownames(exp2_gene_fam_tpms), function(x) substr(x,1,4))
#cross_type_larson<-sapply(rownames(larson_gene_fam_tpms), function(x) substr(x,1,4))

#Make group columns (cross_dataset)
group_exp1<-sapply(cross_type_exp1, function(x) paste0(x,"_exp1"))
group_exp2<-sapply(cross_type_exp2, function(x) paste0(x,"_exp2"))
#group_larson<-sapply(cross_type_larson, function(x) paste0(x,"_larson"))

#Add columns of things to test in linear models (XY mismatch vs not, categories of XY mismatch direction, sterile vs fertile, F1 vs non-F1 hybrid background)
mismatch_exp1<-sapply(cross_type_exp1, function(x) if(x %in% c("CPLY","WLPY","PPLL","LLPP")) "yes" else "no")
mismatch_exp2<-sapply(cross_type_exp2, function(x) if(x %in% c("CPLY","WLPY","PPLL","LLPP")) "yes" else "no")
#mismatch_larson<-sapply(cross_type_larson, function(x) if(x %in% c("CPLY","WLPY","PPLL","LLPP")) "yes" else "no")
mismatchDir_exp1<-sapply(cross_type_exp1, function(x) if(x %in% c("CPLY","PPLL")) "musXdomY" else if(x %in% c("WLPY","LLPP")) "domXmusY" else "none")
mismatchDir_exp2<-sapply(cross_type_exp2, function(x) if(x %in% c("CPLY","PPLL")) "musXdomY" else if(x %in% c("WLPY","LLPP")) "domXmusY" else "none")
#mismatchDir_larson<-sapply(cross_type_larson, function(x) if(x %in% c("CPLY","PPLL")) "musXdomY" else if(x %in% c("WLPY","LLPP")) "domXmusY" else "none")
sterile_exp1<-sapply(cross_type_exp1, function(x) if(x %in% c("PLPY","PPLL")) "yes" else "no")
sterile_exp2<-sapply(cross_type_exp2, function(x) if(x %in% c("PLPY","PPLL")) "yes" else "no")
#sterile_larson<-sapply(cross_type_larson, function(x) if(x %in% c("PLPY","PPLL")) "yes" else "no")
F1hybBG_exp1<-sapply(cross_type_exp1, function(x) if(x %in% c("LPLY","PLPY","PPLL","LLPP")) "yes" else "no")
F1hybBG_exp2<-sapply(cross_type_exp2, function(x) if(x %in% c("LPLY","PLPY","PPLL","LLPP")) "yes" else "no")
#F1hybBG_larson<-sapply(cross_type_larson, function(x) if(x %in% c("LPLY","PLPY","PPLL","LLPP")) "yes" else "no")

exp1_df<-as.data.frame(cbind(exp1_gene_fam_tpms, copy_number=copy_number_exp1, cross_type=cross_type_exp1, group=group_exp1, mismatch=mismatch_exp1, mismatchDir=mismatchDir_exp1, sterile=sterile_exp1, F1hybBG=F1hybBG_exp1))
exp2_df<-as.data.frame(cbind(exp2_gene_fam_tpms, copy_number=copy_number_exp2, cross_type=cross_type_exp2, group=group_exp2, mismatch=mismatch_exp2, mismatchDir=mismatchDir_exp2, sterile=sterile_exp2, F1hybBG=F1hybBG_exp2))
#larson_df<-as.data.frame(cbind(larson_gene_fam_tpms, copy_number=copy_number_larson, cross_type=cross_type_larson, group=group_larson, mismatch=mismatch_larson, mismatchDir=mismatchDir_larson, sterile=sterile_larson, F1hybBG=F1hybBG_larson))
print(exp1_df)
print(exp2_df)
#print(larson_df)


mydf<-as.data.frame(rbind(exp1_df,exp2_df)) #,larson_df))
print(mydf)

#Plot
#Colors: CCPP 9B1D20, CCPPLY C86756, LLPP 7C898B (OLD: A5535A), LLPPLY DFAEB4 (OLD:813E5D), PPLL BBBDF6 (OLD: 391463), PPLLPY 053B06 (OLD: 3A0842), WWLL 023C8D, WWLLPY 59A5B1
pdf(paste("geneFam_expression_by_copyNumber",cell_type,"Yintro_expOnly.pdf", sep="."))
for(g in gene_fams){
	my_df<-as.data.frame(mydf[which(mydf$gene_family==g),])
	print(head(my_df))
	cn_range<-diff(range(my_df$copy_number))
	p<-ggplot(my_df, aes(x=copy_number, y=as.numeric(as.character(expression)), fill=cross_type, color=cross_type)) + geom_point(size=4, position = position_dodge(width=cn_range/25)) 
	p<-p + stat_summary(fun="mean", geom="point", size=10, shape="+", position = position_dodge(width=cn_range/25))
	p<-p + stat_summary(aes(x=copy_number, y=as.numeric(as.character(expression))), geom="errorbar", fun.data="mean_sdl", fun.args=list(mult=1), size=2, width=cn_range/10, position = position_dodge(width=cn_range/25))
	p<-p + labs(title=paste("Expression vs Copy Number - ",g), x="Copy Number", y="Median Gene Family TPM")
	p<-p + theme(axis.text=element_text(size=18), axis.title=element_text(size=21), plot.title=element_text(size=28))
	p<-p + theme_minimal() + scale_color_manual(values=c("#9B1D20","#C86756","#7C898B","#DFAEB4","#BBBDF6","#053B06","#59A5B1","#023C8D"))
	print(p)
}
dev.off()

#Linear models
print("Linear models: effect of XY mismatch\n")
for(g in gene_fams){
	print(paste0("lmer test for ",g,":"))
	null_model<-lmer(as.numeric(as.character(expression)) ~ as.numeric(as.character(copy_number)) + (1|dataset), data=mydf[which(mydf$gene_family==g),])
	with_XYmismatch<-lmer(as.numeric(as.character(expression)) ~ as.numeric(as.character(copy_number)) + mismatch + (1|dataset), data=mydf[which(mydf$gene_family==g),])
	print(anova(with_XYmismatch, null_model))
	print(summary(with_XYmismatch))
}

print("Linear models: effect of XY mismatch DIRECTION")
for(g in gene_fams){
	print(paste0("lmer test for ",g,":"))
	null_model<-lmer(as.numeric(as.character(expression)) ~ as.numeric(as.character(copy_number)) + (1|dataset), data=mydf[which(mydf$gene_family==g),])
	with_XYmismatchDir<-lmer(as.numeric(as.character(expression)) ~ as.numeric(as.character(copy_number)) + mismatchDir + (1|dataset), data=mydf[which(mydf$gene_family==g),])
	print(anova(with_XYmismatchDir, null_model))
	print(summary(with_XYmismatchDir))
}

print("\nLinear models: effect of sterile background (PPLL and PPLLPY)\n")
for(g in gene_fams){
	print(paste0("lmer test for ",g,":"))
	null_model<-lmer(as.numeric(as.character(expression)) ~ as.numeric(as.character(copy_number)) + (1|dataset), data=mydf[which(mydf$gene_family==g),])
	with_sterile<-lmer(as.numeric(as.character(expression)) ~ as.numeric(as.character(copy_number)) + sterile + (1|dataset), data=mydf[which(mydf$gene_family==g),])
	print(anova(with_sterile, null_model))
	print(summary(with_sterile))
}

print("\nLinear models: effect of F1 hybrid background\n")
for(g in gene_fams){
	print(paste0("lmer test for ",g,":"))
	null_model<-lmer(as.numeric(as.character(expression)) ~ as.numeric(as.character(copy_number)) + (1|dataset), data=mydf[which(mydf$gene_family==g),])
	with_F1hybBG<-lmer(as.numeric(as.character(expression)) ~ as.numeric(as.character(copy_number)) + F1hybBG + (1|dataset), data=mydf[which(mydf$gene_family==g),])
	print(anova(with_F1hybBG, null_model))
	print(summary(with_F1hybBG))
}

print("Null models:")
for(g in gene_fams){
	print(paste0("lmer test for ",g,":"))
	null_model<-lmer(as.numeric(as.character(expression)) ~ as.numeric(as.character(copy_number)) + (1|dataset), data=mydf[which(mydf$gene_family==g),])
	print(summary(null_model))
}

print("Done with 09_plot_exp_by_copyNum.r")

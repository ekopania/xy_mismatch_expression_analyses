#PURPOSE: Use the R package chromPlot to plot DE genes over chromosome images

library(chromPlot) #Note: loads biomaRt

dataset<-"Yintro_exp2_sameY"
#Options: Yintro_exp1_difY (CCPP vs CCPPLY and WWLL vs WWLLPY); Yintro_exp1_sameY (CCPP vs WWLLPY and WWLL vs CCPPLY); Yintro_exp2_difY (LLPP vs LLPPLY and PPLL vs PPLLPY); Yintro_exp2_sameY (LLPP vs PPLLPY and PPLL vs LLPPLY)
cell_type<-"RS"

#Load DE gene data
if(dataset == "Yintro_exp1_difY"){
	samp1_DEgenes<-read.table(paste0("topDEgenes.lrt.DE.CCPP_",cell_type,"vsCCPPLY_",cell_type,".txt"), header=TRUE)
	samp2_DEgenes<-read.table(paste0("topDEgenes.lrt.DE.WWLL_",cell_type,"vsWWLLPY_",cell_type,".txt"), header=TRUE)
} else if(dataset == "Yintro_exp1_sameY"){
	samp1_DEgenes<-read.table(paste0("topDEgenes.lrt.DE.CCPP_",cell_type,"vsWWLLPY_",cell_type,".txt"), header=TRUE)
	samp2_DEgenes<-read.table(paste0("topDEgenes.lrt.DE.CCPPLY_",cell_type,"vsWWLL_",cell_type,".txt"), header=TRUE)
} else if(dataset == "Yintro_exp2_difY"){
	samp1_DEgenes<-read.table(paste0("topDEgenes.lrt.DE.LLPP_",cell_type,"vsLLPPLY_",cell_type,".txt"), header=TRUE)
	samp2_DEgenes<-read.table(paste0("topDEgenes.lrt.DE.PPLL_",cell_type,"vsPPLLPY_",cell_type,".txt"), header=TRUE)
} else if(dataset == "Yintro_exp2_sameY"){
	samp1_DEgenes<-read.table(paste0("topDEgenes.lrt.DE.LLPPLY_",cell_type,"vsPPLL_",cell_type,".txt"), header=TRUE)
	samp2_DEgenes<-read.table(paste0("topDEgenes.lrt.DE.LLPP_",cell_type,"vsPPLLPY_",cell_type,".txt"), header=TRUE)
} else{
	stop("Invalid dataset: must be 'Yintro_exp1_difY', 'Yintro_exp1_sameY', 'Yintro_exp2_difY', or 'Yintro_exp2_sameY'")
}

#Get gene coordinates
ens_mus<-useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl")
all_genes<-getBM(attributes=c('ensembl_gene_id','chromosome_name','start_position','end_position'), mart=ens_mus)
colnames(all_genes)<-c("GeneID", "Chrom", "Start", "End")
print(head(all_genes))

#Append gene coordinates and logFC values for DE genes
samp1_stats<-c()
for(i in 1:nrow(samp1_DEgenes)){
	if(rownames(samp1_DEgenes)[i] %in% all_genes$GeneID){
		this_line<-cbind(all_genes[which(all_genes$GeneID == rownames(samp1_DEgenes)[i]),], samp1_DEgenes$logFC[i], samp1_DEgenes$direc[i])
		samp1_stats<-rbind(samp1_stats, this_line)
	}
}
samp1_df<-as.data.frame(lapply(samp1_stats, unlist))
colnames(samp1_df)<-c("GeneID", "Chrom", "Start", "End","LogFC", "Direc")
samp1_df$LogFC<-as.numeric(as.character(samp1_df$LogFC))

samp2_stats<-c()
for(i in 1:nrow(samp2_DEgenes)){
        if(rownames(samp2_DEgenes)[i] %in% all_genes$GeneID){
                this_line<-cbind(all_genes[which(all_genes$GeneID == rownames(samp2_DEgenes)[i]),], samp2_DEgenes$logFC[i], samp2_DEgenes$direc[i])
                samp2_stats<-rbind(samp2_stats, this_line)
        }
}
samp2_df<-as.data.frame(lapply(samp2_stats, unlist))
colnames(samp2_df)<-c("GeneID", "Chrom", "Start", "End","LogFC", "Direc")
samp2_df$LogFC<-as.numeric(as.character(samp2_df$LogFC))

print(dim(samp1_df))
print(dim(samp2_df))
print(head(samp1_df))
print(head(samp2_df))

#Generate chrom plots
#pdf("chromPlot.test.pdf")
#chromPlot(gaps=all_genes, chr=c(1:19,"X","Y"))
if(dataset == "Yintro_exp1_difY"){
	pdf(paste("chromPlot.CCPPvsCCPPLY.WWLLvsWWLLPY",cell_type,"pdf", sep="."))
	chromPlot(gaps=all_genes, stat=samp1_df, stat2=samp2_df, statCol="LogFC", statCol2="LogFC", statName="musVSmusDomY_LogFC", statName2="domVSdomMusY_LogFC", statTyp="p", chr=c(1:19,"X","Y"), statSumm="none")
} else if(dataset == "Yintro_exp1_sameY"){
        pdf(paste("chromPlot.CCPPvsWWLLPY.WWLLvsCCPPLY",cell_type,"pdf", sep="."))
	chromPlot(gaps=all_genes, stat=samp1_df, stat2=samp2_df, statCol="LogFC", statCol2="LogFC", statName="musVSdomMusY_LogFC", statName2="domVSmusDomY_LogFC", statTyp="p", chr=c(1:19,"X","Y"), statSumm="none")
} else if(dataset == "Yintro_exp2_difY"){
	pdf(paste("chromPlot.LLPPvsLLPPLY.PPLLvsPPLLPY",cell_type,"pdf", sep="."))
	chromPlot(gaps=all_genes, stat=samp1_df, stat2=samp2_df, statCol="LogFC", statCol2="LogFC", statName="domXmusVSdomXmusDomY_LogFC", statName2="musXdomVSmusXdomMusY_LogFC", statTyp="p", chr=c(1:19,"X","Y"), statSumm="none")
} else if(dataset == "Yintro_exp2_sameY"){
	pdf(paste("chromPlot.PPLLvsLLPPLY.LLPPvsPPLLPY",cell_type,"pdf", sep="."))
	chromPlot(gaps=all_genes, stat=samp1_df, stat2=samp2_df, statCol="LogFC", statCol2="LogFC", statName="musXdomVSdomXmusDomY_LogFC", statName2="domXmusVSmusXdomMusY_LogFC", statTyp="p", chr=c(1:19,"X","Y"), statSumm="none")
}
warnings() #Will generate a warning for each chromosome because there is no centromere info, but we don't need that so it's fine
dev.off()

#Repeat while removing things that are DE between mus and dom (so for CCPP vs CCPPLY, remove DE genes on the Y between CCPP and WWLL to tease apart things that are DE because of allelic differences vs disrupted expression)
background_DEgenes_musVdom<-read.table(paste0("topDEgenes.lrt.DE.CCPP_",cell_type,"vsWWLL_",cell_type,".txt"), header=TRUE)
background_DEgenes_hybrids<-read.table(paste0("topDEgenes.lrt.DE.LLPP_",cell_type,"vsPPLL_",cell_type,".txt"), header=TRUE)
if( (dataset == "Yintro_exp1_difY") || (dataset == "Yintro_exp2_difY") ){
	#Get background DE Y genes for comparisons with different Y
	background_DEgenes_musVdom_targetChr<-background_DEgenes_musVdom[which(background_DEgenes_musVdom$chrs == "Y"), ]
	background_DEgenes_hybrids_targetChr<-background_DEgenes_hybrids[which(background_DEgenes_hybrids$chrs == "Y"), ]
	#background_DEgenes_musVdom_targetChr<-background_DEgenes_musVdom
	#background_DEgenes_hybrids_targetChr<-background_DEgenes_hybrids
} else{
	#Get background DE X and auto genes for comparisons with same Y but different X/autos
	background_DEgenes_musVdom_targetChr<-background_DEgenes_musVdom[which(background_DEgenes_musVdom$chrs != "Y"), ]
	background_DEgenes_hybrids_targetChr<-background_DEgenes_hybrids[which(background_DEgenes_hybrids$chrs != "Y"), ]
	#background_DEgenes_musVdom_targetChr<-background_DEgenes_musVdom
	#background_DEgenes_hybrids_targetChr<-background_DEgenes_hybrids
}
#Identify genes DE in the same direction in background vs comparison of interest
#NEED TO MAKE SURE DIRECTIONALITY IS RIGHT
#Exp 1: background DE is CCPP vs WWLL (so + is higher in mus)
	#Comparisons are: CCPPvsCCPPLY, WWLLvsWWLLPY, CCPPvsWWLLPY, CCPPLYvsWWLL
#Exp 2: background DE is LLPP vs PPLL (so + is higher in domXmus)
	#Comparisons are: LLPPvsLLPPLY, PPLLvsPPLLPY, LLPPLYvsPPLL, LLPPvsPPLLPY 
if( (dataset == "Yintro_exp1_difY") || (dataset == "Yintro_exp1_sameY") || (dataset == "Yintro_exp2_difY") ){
	remove1<-c()
	remove2<-c()
	for(i in 1:nrow(background_DEgenes_musVdom_targetChr)){
		this_gene<-rownames(background_DEgenes_musVdom_targetChr)[i]
		if(dataset == "Yintro_exp1_difY"){
			#samp1 is CCPP vs CCPPLY and samp2 is WWLL vs WWLLPY
			if(this_gene %in% samp1_df$GeneID){
				samp1_line<-samp1_df[which(samp1_df$GeneID == this_gene),]
				#background set and experimental comparison are both mus v dom, ==
				if(background_DEgenes_musVdom_targetChr$direc == samp1_line$Direc){
					remove1<-c(remove1, this_gene)
				}
			}
			if(this_gene %in% samp2_df$GeneID){
				samp2_line<-samp2_df[which(samp2_df$GeneID == this_gene),]
				#background set is mus v dom but experimental set is dom v mus, !=
		                if(background_DEgenes_musVdom_targetChr$direc != samp2_line$Direc){
		                        remove2<-c(remove2, this_gene)
		                }
		        }
		}
		if(dataset == "Yintro_exp1_sameY"){
			#samp1 is CCPP vs WWLLPY and samp 2 is  CCPPLYvsWWLL
			if(this_gene %in% samp1_df$GeneID){
	                        samp1_line<-samp1_df[which(samp1_df$GeneID == this_gene),]
	                        #background set and experimental comparison are both mus v dom, ==
	                        if(background_DEgenes_musVdom_targetChr$direc == samp1_line$Direc){
	                                remove1<-c(remove1, this_gene)
	                        }
	                }
	                if(this_gene %in% samp2_df$GeneID){
	                        samp2_line<-samp2_df[which(samp2_df$GeneID == this_gene),]
	                        #background set and experimental comparison are both mus v dom, ==
	                        if(background_DEgenes_musVdom_targetChr$direc == samp2_line$Direc){
	                                remove2<-c(remove2, this_gene)
	                        }
	                }
		}
		if(dataset == "Yintro_exp2_difY"){
			#samp1 is LLPP vs LLPPLY and samp2 is PPLL vs PPLLPY
			#Backgrounds are the same and we just want to account for strain Y differences here, so using the mus vs dom DE instead of LLPP vs PPLL DE
			if(this_gene %in% samp1_df$GeneID){
                                samp1_line<-samp1_df[which(samp1_df$GeneID == this_gene),]
                                #background set and experimental comparison are both mus v dom, ==
                                if(background_DEgenes_musVdom_targetChr$direc == samp1_line$Direc){
                                        remove1<-c(remove1, this_gene)
                                }
                        }
			if(this_gene %in% samp2_df$GeneID){
                                samp2_line<-samp2_df[which(samp2_df$GeneID == this_gene),]
                                #background set is mus v dom but experimental set is dom v mus, !=
                                if(background_DEgenes_musVdom_targetChr$direc != samp2_line$Direc){
                                        remove2<-c(remove2, this_gene)
                                }
                        }
		}
	}
} else{
	#dataset is Yintro_exp2_sameY
	remove1<-c()
        remove2<-c()
        for(i in 1:nrow(background_DEgenes_hybrids_targetChr)){
                this_gene<-rownames(background_DEgenes_hybrids_targetChr)[i]
		#samp1 is LLPPLYvsPPLL and samp2 is LLPPvsPPLLPY
		if(this_gene %in% samp1_df$GeneID){
			samp1_line<-samp1_df[which(samp1_df$GeneID == this_gene),]
			#background set and experimental comparison are both LLPP v PPLL, ==
			if(background_DEgenes_hybrids_targetChr$direc == samp1_line$Direc){
				remove1<-c(remove1, this_gene)
			}
		}
		if(this_gene %in% samp2_df$GeneID){
			samp2_line<-samp2_df[which(samp2_df$GeneID == this_gene),]
			#background set and experimental comparison are both LLPP v PPLL, ==
			if(background_DEgenes_hybrids_targetChr$direc == samp2_line$Direc){
				remove2<-c(remove2, this_gene)
			}
		}
	}
}
	
#Remove background DE genes from input dataframes
samp1_df_filtered<-samp1_df[which(!(samp1_df$GeneID %in% remove1)),]
samp2_df_filtered<-samp2_df[which(!(samp2_df$GeneID %in% remove2)),]
print(dim(samp1_df_filtered))
print(dim(samp2_df_filtered))
print(head(samp1_df_filtered))
print(head(samp2_df_filtered))
#Generate chrom plots
if(dataset == "Yintro_exp1_difY"){
        pdf(paste("chromPlot.CCPPvsCCPPLY.WWLLvsWWLLPY",cell_type,"backgroundDEremoved.directional.pdf", sep="."))
        chromPlot(gaps=all_genes, stat=samp1_df_filtered, stat2=samp2_df_filtered, statCol="LogFC", statCol2="LogFC", statName="musVSmusDomY_LogFC", statName2="domVSdomMusY_LogFC", statTyp="p", chr=c(1:19,"X","Y"), statSumm="none")
} else if(dataset == "Yintro_exp1_sameY"){
        pdf(paste("chromPlot.CCPPvsWWLLPY.WWLLvsCCPPLY",cell_type,"backgroundDEremoved.directional.pdf", sep="."))
        chromPlot(gaps=all_genes, stat=samp1_df_filtered, stat2=samp2_df_filtered, statCol="LogFC", statCol2="LogFC", statName="musVSdomMusY_LogFC", statName2="domVSmusDomY_LogFC", statTyp="p", chr=c(1:19,"X","Y"), statSumm="none")
} else if(dataset == "Yintro_exp2_difY"){
        pdf(paste("chromPlot.LLPPvsLLPPLY.PPLLvsPPLLPY",cell_type,"backgroundDEremoved.directional.pdf", sep="."))
        chromPlot(gaps=all_genes, stat=samp1_df_filtered, stat2=samp2_df_filtered, statCol="LogFC", statCol2="LogFC", statName="domXmusVSdomXmusDomY_LogFC", statName2="musXdomVSmusXdomMusY_LogFC", statTyp="p", chr=c(1:19,"X","Y"), statSumm="none")
} else if(dataset == "Yintro_exp2_sameY"){
        pdf(paste("chromPlot.PPLLvsLLPPLY.LLPPvsPPLLPY",cell_type,"backgroundDEremoved.directional.pdf", sep="."))
        chromPlot(gaps=all_genes, stat=samp1_df_filtered, stat2=samp2_df_filtered, statCol="LogFC", statCol2="LogFC", statName="musXdomVSdomXmusDomY_LogFC", statName2="domXmusVSmusXdomMusY_LogFC", statTyp="p", chr=c(1:19,"X","Y"), statSumm="none")
}
warnings() #Will generate a warning for each chromosome because there is no centromere info, but we don't need that so it's fine
dev.off()
#Print #s of DE genes:
print("Numbers of DE genes, after controlling for background DE between subspecies:")
aDE1<-samp1_df_filtered[which(samp1_df_filtered$Chrom %in% c(1:19)), ]
aDE2<-samp2_df_filtered[which(samp2_df_filtered$Chrom %in% c(1:19)), ]
xDE1<-samp1_df_filtered[which(samp1_df_filtered$Chrom=="X"), ]
xDE2<-samp2_df_filtered[which(samp2_df_filtered$Chrom=="X"), ]
yDE1<-samp1_df_filtered[which(samp1_df_filtered$Chrom=="Y"), ]
yDE2<-samp2_df_filtered[which(samp2_df_filtered$Chrom=="Y"), ]
print(paste("Comparison 1 # autosomal DE genes:", nrow(aDE1)))
print(paste("	Higher in first cross type on autos:", length(which(aDE1$LogFC > 0))))
print(paste("   Higher in second cross type on autos:", length(which(aDE1$LogFC < 0))))
print(paste("Comparison 1 # X-linked DE genes:", nrow(xDE1)))
print(paste("   Higher in first cross type on X:", length(which(xDE1$LogFC > 0))))
print(paste("   Higher in second cross type on X:", length(which(xDE1$LogFC < 0))))
print(paste("Comparison 1 # Y-linked DE genes:", nrow(yDE1)))
print(paste("   Higher in first cross type on Y:", length(which(yDE1$LogFC > 0))))
print(paste("   Higher in second cross type on Y:", length(which(yDE1$LogFC < 0))))
print(paste("Comparison 2 # autosomal DE genes:", nrow(aDE2)))
print(paste("   Higher in first cross type on autos:", length(which(aDE2$LogFC > 0))))
print(paste("   Higher in second cross type on autos:", length(which(aDE2$LogFC < 0))))
print(paste("Comparison 2 # X-linked DE genes:", nrow(xDE2)))
print(paste("   Higher in first cross type on X:", length(which(xDE2$LogFC > 0))))
print(paste("   Higher in second cross type on X:", length(which(xDE2$LogFC < 0))))
print(paste("Comparison 2 # Y-linked DE genes:", nrow(yDE2)))
print(paste("   Higher in first cross type on Y:", length(which(yDE2$LogFC > 0))))
print(paste("   Higher in second cross type on Y:", length(which(yDE2$LogFC < 0))))

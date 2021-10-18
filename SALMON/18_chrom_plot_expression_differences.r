#PURPOSE: Generate chrom plot bar graphs showing expression differences for DE genes on each chromosome between XY mismatch and non-mismatch or between hybrids and not hybrids (Similar to Larson et al. 2017 Fig. 1)

#library(chromPlot)
library(biomaRt)
library(karyoploteR)
library(matrixStats)

dataset<-"Yintro_exp2"
cell_type<-"RS"

#Read in data
fpkm_data<-read.table(paste("fpkm_filtered_table",dataset,cell_type,"97.txt", sep="."), header=TRUE)

if(dataset=="Yintro_exp1"){
	#Set up with sex chr mismatch always first
	cross_types<-list(c("CPLY","CCPP"), c("WLPY","CCPP"), c("WLPY","WWLL"), c("CPLY","WWLL"))
	names(cross_types)<-c("musX","musY","domX","domY")
	if(cell_type=="LZ"){
                musX<-read.table("topDEgenes.lrt.DE.CCPP_LZvsCCPPLY_LZ.withGeneNames.txt", header=TRUE, fill=TRUE)
                musY<-read.table("topDEgenes.lrt.DE.CCPP_LZvsWWLLPY_LZ.withGeneNames.txt", header=TRUE, fill=TRUE)
                domX<-read.table("topDEgenes.lrt.DE.WWLL_LZvsWWLLPY_LZ.withGeneNames.txt", header=TRUE, fill=TRUE)
                domY<-read.table("topDEgenes.lrt.DE.CCPPLY_LZvsWWLL_LZ.withGeneNames.txt", header=TRUE, fill=TRUE)
        } else if(cell_type=="RS"){
		musX<-read.table("topDEgenes.lrt.DE.CCPP_RSvsCCPPLY_RS.withGeneNames.txt", header=TRUE, fill=TRUE)
		musY<-read.table("topDEgenes.lrt.DE.CCPP_RSvsWWLLPY_RS.withGeneNames.txt", header=TRUE, fill=TRUE)
		domX<-read.table("topDEgenes.lrt.DE.WWLL_RSvsWWLLPY_RS.withGeneNames.txt", header=TRUE, fill=TRUE)
		domY<-read.table("topDEgenes.lrt.DE.CCPPLY_RSvsWWLL_RS.withGeneNames.txt", header=TRUE, fill=TRUE)
        } else{
		stop("Invalid cell type: must be 'LZ' or 'RS'")
	}
} else if(dataset=="Yintro_exp2"){
	cross_types<-list(c("PPLL","PLPY"), c("LLPP","PLPY"), c("LLPP","LPLY"), c("LPLY","PPLL"))
	names(cross_types)<-c("musX","musY","domX","domY")
	if(cell_type=="LZ"){
                musX<-read.table("topDEgenes.lrt.DE.PPLL_LZvsPPLLPY_LZ.withGeneNames.txt", header=TRUE, fill=TRUE)
                musY<-read.table("topDEgenes.lrt.DE.LLPP_LZvsPPLLPY_LZ.withGeneNames.txt", header=TRUE, fill=TRUE)
                domX<-read.table("topDEgenes.lrt.DE.LLPP_LZvsLLPPLY_LZ.withGeneNames.txt", header=TRUE, fill=TRUE)
                domY<-read.table("topDEgenes.lrt.DE.LLPPLY_LZvsPPLL_LZ.withGeneNames.txt", header=TRUE, fill=TRUE)
        } else if(cell_type=="RS"){
		musX<-read.table("topDEgenes.lrt.DE.PPLL_RSvsPPLLPY_RS.withGeneNames.txt", header=TRUE, fill=TRUE)
                musY<-read.table("topDEgenes.lrt.DE.LLPP_RSvsPPLLPY_RS.withGeneNames.txt", header=TRUE, fill=TRUE)
                domX<-read.table("topDEgenes.lrt.DE.LLPP_RSvsLLPPLY_RS.withGeneNames.txt", header=TRUE, fill=TRUE)
                domY<-read.table("topDEgenes.lrt.DE.LLPPLY_RSvsPPLL_RS.withGeneNames.txt", header=TRUE, fill=TRUE)
        } else{
                stop("Invalid cell type: must be 'LZ' or 'RS'")
        }
} else{
	stop("Invalid dataset: must be 'Yintro_exp1' or 'Yintro_exp2'")
}

print("Here are the cross type comparisons:")
print(cross_types)
#print(cross_types[["musX"]])
#print(cross_types[["musX"]][1])

#Get gene coordinates
ens_mus<-useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl")
all_genes<-getBM(attributes=c('ensembl_gene_id','chromosome_name','start_position','end_position'), mart=ens_mus)
colnames(all_genes)<-c("geneID", "seqnames", "start", "end")
print(head(all_genes))
#Got chromosome lengths from Mus_musculus.GRCm38.dna.primary_assembly.dict
chr_lens<-c(195471971, 182113224, 160039680, 156508116, 151834684, 149736546, 145441459, 129401213, 124595110, 130694993, 122082543, 120129022, 120421639, 124902244, 104043685, 98207768, 94987271, 90702639, 61431566, 171031299, 91744698) 
custom.mus.genome <- toGRanges(data.frame(chr=c(1:19,"X","Y"), start=rep(1, 21), end=chr_lens))

#Loop through each comparison
pdf(paste("chromPlot.FPKMdiff",dataset,cell_type,"pdf", sep="."), height=8.5, width=11, onefile=TRUE)
for(c in names(cross_types)){
	print(paste("Working on comparison for:", c))
	#Get difference in median FPKM between XY mismatch and non mismatch cross types, for DE genes
	#Separate by direction (higher with XY mismatch or lower with XY mismatch)
	fpkm1<-as.matrix(fpkm_data[, grepl(cross_types[[c]][1], colnames(fpkm_data))])
	fpkm2<-as.matrix(fpkm_data[, grepl(cross_types[[c]][2], colnames(fpkm_data))])
	fpkm_med1<-rowMedians(fpkm1)
	fpkm_med2<-rowMedians(fpkm2)
	names(fpkm_med1)<-rownames(fpkm1)
	names(fpkm_med2)<-rownames(fpkm2)
	DE<-get(c)
	fpkm_med1_DE<-fpkm_med1[which(names(fpkm_med1) %in% rownames(DE))]
	fpkm_med2_DE<-fpkm_med2[which(names(fpkm_med2) %in% rownames(DE))]
	stopifnot(all.equal(names(fpkm_med1_DE), names(fpkm_med2_DE)))
	med_fpkm_diff<-fpkm_med1_DE - fpkm_med2_DE
	higher_mismatch<-med_fpkm_diff[which(med_fpkm_diff > 0)]
	lower_mismatch<-med_fpkm_diff[which(med_fpkm_diff < 0)]
	print(head(higher_mismatch))
	print(head(lower_mismatch))
	#Get chromosome coordinates for DE genes, append fpkm differences to make df for chromPlot
	my_gaps_higher<-all_genes[which(all_genes$geneID %in% names(higher_mismatch)),]
	my_gaps_lower<-all_genes[which(all_genes$geneID %in% names(lower_mismatch)),]
	print(head(my_gaps_higher))
	print(head(my_gaps_lower))
	my_gaps_higher_sorted<-my_gaps_higher[order(match(names(higher_mismatch), rownames(my_gaps_higher))),]
	my_gaps_lower_sorted<-my_gaps_lower[order(match(names(lower_mismatch), rownames(my_gaps_lower))),]
	my_gaps_higher_wDiff<-as.data.frame(cbind(my_gaps_higher_sorted, fpkm_dif=higher_mismatch))
	#Makes these positive numbers so that chromPlot puts them facing away from the chromosome, to the left
	my_gaps_lower_wDiff<-as.data.frame(cbind(my_gaps_lower_sorted, fpkm_dif=-lower_mismatch))
	print(head(my_gaps_higher_wDiff))
	print(head(my_gaps_lower_wDiff))
	my_gaps_higher_wDiff_noNA<-my_gaps_higher_wDiff[which(!(is.na(my_gaps_higher_wDiff$start))),]
	my_gaps_lower_wDiff_noNA<-my_gaps_lower_wDiff[which(!(is.na(my_gaps_lower_wDiff$start))),]
	#Generate chromPlot
	#chromPlot(gaps=all_genes, stat=my_gaps_lower_wDiff, stat2=my_gaps_higher_wDiff, statCol="fpkm_dif", statCol2="fpkm_dif", statName="lower_in_XYmismatch", statName2="higher_in_XYmismatch", statTyp="b", chr=c(1:19,"X","Y"), statSumm="none", figCols=1, title=paste(cross_types[[c]][1], "vs", cross_types[[c]][2]), legChrom=NA)
	#chromPlot(gaps=all_genes, annot1=my_gaps_lower_wDiff, annot2=my_gaps_higher_wDiff, statCol="fpkm_dif", statCol2="fpkm_dif", statName="lower_in_XYmismatch", statName2="higher_in_XYmismatch", statTyp="l", chr=c(1:19,"X","Y"), statSumm="none", figCols=1, title=paste(cross_types[[c]][1], "vs", cross_types[[c]][2]), legChrom="Y")
	#Generate karyoPlot
	#print(head(my_gaps_higher_wDiff_noNA))
	GRange_higher<-toGRanges(data.frame(chr=my_gaps_higher_wDiff_noNA$seqnames, start=my_gaps_higher_wDiff_noNA$start, end=my_gaps_higher_wDiff_noNA$end, y=my_gaps_higher_wDiff_noNA$fpkm_dif))
	print(head(GRange_higher))
	GRange_lower<-toGRanges(data.frame(chr=my_gaps_lower_wDiff_noNA$seqnames, start=my_gaps_lower_wDiff_noNA$start, end=my_gaps_lower_wDiff_noNA$end, y=my_gaps_lower_wDiff_noNA$fpkm_dif))
	print(head(GRange_lower))
	kp<-plotKaryotype(genome=custom.mus.genome, chromosomes=c(1:19, "X", "Y"), plot.type=3)
	kpAddMainTitle(kp, main=paste(cross_types[[c]][1],"vs",cross_types[[c]][2]))
	#kpDataBackground(kp, data.panel=1, r0=0, r1=10^ceiling(log10(max(GRange_higher$y))))
	kpAxis(kp, data.panel=1, ymax=10^ceiling(log10(max(GRange_higher$y))))
	kpBars(kp, data=GRange_higher, y1=GRange_higher$y, col="black", ymax=max(GRange_higher$y), data.panel=1)
	#kpDataBackground(kp, data.panel=2, r0=0, r1=10^ceiling(log10(max(GRange_lower$y))))
	kpAxis(kp, data.panel=2, ymax=10^ceiling(log10(max(GRange_lower$y))))
	kpBars(kp, data=GRange_lower, y1=GRange_lower$y, col="red", ymax=max(GRange_lower$y), data.panel=2)
}
dev.off()

#CHANGES: stats 1 and 2 will be negative and positive for same dataset; plot difference in FPKM (tpm?) not logFC
#chromPlot(gaps=all_genes, stat=samp1_df, stat2=samp2_df, statCol="LogFC", statCol2="LogFC", statName="musVSmusDomY_LogFC", statName2="domVSdomMusY_LogFC", statTyp="p", chr=c(1:19,"X","Y"), statSumm="none")

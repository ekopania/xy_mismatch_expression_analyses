#PURPOSE: Generate chrom plot bar graphs showing expression differences for DE genes on each chromosome between XY mismatch and non-mismatch or between hybrids and not hybrids (Similar to Larson et al. 2017 Fig. 1)

#library(chromPlot)
library(biomaRt)
library(karyoploteR)
library(matrixStats)
library(vegan)

dataset<-"LarsonEtal"
cell_type<-"LZ"

#Read in data
fpkm_data<-read.table(paste("fpkm_filtered_table",dataset,cell_type,"97.txt", sep="."), header=TRUE)

if(dataset=="LarsonEtal"){
	#Set up with hybrid always first
	cross_types<-list(c("PPLL","CCPP"), c("LLPP","CCPP"), c("LLPP","WWLL"), c("PPLL","WWLL"))
	names(cross_types)<-c("musX","musY","domX","domY")
	if(cell_type=="LZ"){
                musX<-read.table("topDEgenes.lrt.DE.CCPP_LZvsPPLL_LZ.txt", header=TRUE, fill=TRUE)
                musY<-read.table("topDEgenes.lrt.DE.CCPP_LZvsLLPP_LZ.txt", header=TRUE, fill=TRUE)
                domX<-read.table("topDEgenes.lrt.DE.LLPP_LZvsWWLL_LZ.txt", header=TRUE, fill=TRUE)
                domY<-read.table("topDEgenes.lrt.DE.PPLL_LZvsWWLL_LZ.txt", header=TRUE, fill=TRUE)
        } else if(cell_type=="RS"){
		musX<-read.table("topDEgenes.lrt.DE.CCPP_RSvsPPLL_RS.txt", header=TRUE, fill=TRUE)
                musY<-read.table("topDEgenes.lrt.DE.CCPP_RSvsLLPP_RS.txt", header=TRUE, fill=TRUE)
                domX<-read.table("topDEgenes.lrt.DE.LLPP_RSvsWWLL_RS.txt", header=TRUE, fill=TRUE)
                domY<-read.table("topDEgenes.lrt.DE.PPLL_RSvsWWLL_RS.txt", header=TRUE, fill=TRUE)
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
	#print(head(higher_mismatch))
	#print(head(lower_mismatch))
	#Get chromosome coordinates for DE genes, append fpkm differences to make df for chromPlot
	my_gaps_higher<-all_genes[which(all_genes$geneID %in% names(higher_mismatch)),]
	my_gaps_lower<-all_genes[which(all_genes$geneID %in% names(lower_mismatch)),]
	#print(head(my_gaps_higher))
	#print(head(my_gaps_lower))
	my_gaps_higher_sorted<-my_gaps_higher[match(names(higher_mismatch), my_gaps_higher$geneID),]
	my_gaps_lower_sorted<-my_gaps_lower[match(names(lower_mismatch), my_gaps_lower$geneID),]
	my_gaps_higher_wDiff<-as.data.frame(cbind(my_gaps_higher_sorted, fpkm_dif=higher_mismatch))
	#Makes these positive numbers so that chromPlot puts them facing away from the chromosome, to the left
	my_gaps_lower_wDiff<-as.data.frame(cbind(my_gaps_lower_sorted, fpkm_dif=-lower_mismatch))
	#print(head(my_gaps_higher_wDiff))
	#print(head(my_gaps_lower_wDiff))
	my_gaps_higher_wDiff_noNA<-my_gaps_higher_wDiff[which(!(is.na(my_gaps_higher_wDiff$start))),]
	my_gaps_lower_wDiff_noNA<-my_gaps_lower_wDiff[which(!(is.na(my_gaps_lower_wDiff$start))),]
	higher_final<-my_gaps_higher_wDiff_noNA<-my_gaps_higher_wDiff_noNA[which(my_gaps_higher_wDiff_noNA$seqnames %in% c(1:19,"X","Y")),]
        lower_final<-my_gaps_lower_wDiff_noNA<-my_gaps_lower_wDiff_noNA[which(my_gaps_lower_wDiff_noNA$seqnames %in% c(1:19,"X","Y")),]
        higher_norm<-decostand(higher_final$fpkm_dif, "range")[,1]
        lower_norm<-decostand(lower_final$fpkm_dif, "range")[,1]
	#Generate karyoPlot
	GRange_higher<-toGRanges(data.frame(chr=higher_final$seqnames, start=higher_final$start, end=higher_final$end, y=higher_norm))
	print(head(GRange_higher))
	GRange_lower<-toGRanges(data.frame(chr=lower_final$seqnames, start=lower_final$start, end=lower_final$end, y=lower_norm))
	print(head(GRange_lower))
	#max_any<-max(c(GRange_higher$y, GRange_lower$y))
	#nearest_hundred<-round(max_any+50,-2)
	kp<-plotKaryotype(genome=custom.mus.genome, chromosomes=c(1:19, "X", "Y"), plot.type=3)
	kpAddMainTitle(kp, main=paste(cross_types[[c]][1],"vs",cross_types[[c]][2]))
	kpAxis(kp, data.panel=1) #, ymax=nearest_hundred
	kpBars(kp, data=GRange_higher, y1=GRange_higher$y, col="black", data.panel=1) #, ymax=max_any
	kpAxis(kp, data.panel=2) #, ymax=nearest_hundred
	kpBars(kp, data=GRange_lower, y1=GRange_lower$y, col="black", data.panel=2) #, ymax=max_any
}
dev.off()

#Compare to F1 hybrids to both mus and dom
pdf(paste("chromPlot.FPKMdiff.vsMusAndDom",dataset,cell_type,"pdf", sep="."), height=8.5, width=11, onefile=TRUE)
mismatch<-c("PPLL","LLPP")
for(m in mismatch){
	fpkm_mismatch<-as.matrix(fpkm_data[, grepl(m, colnames(fpkm_data))])
	fpkm_pure<-as.matrix(fpkm_data[, grepl("CCPP|WWLL", colnames(fpkm_data))])
	fpkm_med_mismatch<-rowMedians(fpkm_mismatch)
	fpkm_med_pure<-rowMedians(fpkm_pure)
	names(fpkm_med_mismatch)<-rownames(fpkm_mismatch)
	names(fpkm_med_pure)<-rownames(fpkm_pure)
	if(m=="PPLL"){
		#Get DE genes in PPLL compared to BOTH CCPP and WWLL
		DE_mismatch_higher<-intersect(rownames(musX[which(musX$direc=="-"),]), rownames(domY[which(domY$direc=="+"),]))
		DE_mismatch_lower<-intersect(rownames(musX[which(musX$direc=="+"),]), rownames(domY[which(domY$direc=="-"),]))
		DE_mismatch<-c(DE_mismatch_higher, DE_mismatch_lower)
	} else{
		#Get DE genes in LLPP compared to BOTH CCPP and WWLL
		DE_mismatch_higher<-intersect(rownames(musY[which(musY$direc=="-"),]), rownames(domX[which(domX$direc=="+"),]))
		DE_mismatch_lower<-intersect(rownames(musY[which(musY$direc=="+"),]), rownames(domX[which(domX$direc=="-"),]))
		DE_mismatch<-c(DE_mismatch_higher, DE_mismatch_lower)
	}
	fpkm_med_mismatch_DE<-fpkm_med_mismatch[which(names(fpkm_med_mismatch) %in% DE_mismatch)]
	fpkm_med_pure_VSmismatch<-fpkm_med_pure[which(names(fpkm_med_pure) %in% DE_mismatch)]
	stopifnot(all.equal(names(fpkm_med_mismatch_DE), names(fpkm_med_pure_VSmismatch)))
	med_fpkm_diff_mismatch<-fpkm_med_mismatch_DE - fpkm_med_pure_VSmismatch
	higher_mismatch<-med_fpkm_diff_mismatch[which(med_fpkm_diff_mismatch > 0)]
	lower_mismatch<-med_fpkm_diff_mismatch[which(med_fpkm_diff_mismatch < 0)]
	print(head(higher_mismatch))
	print(head(lower_mismatch))
	#Get chromosome coordinates for DE genes, append fpkm differences to make df for chromPlot
	my_gaps_higher<-all_genes[which(all_genes$geneID %in% names(higher_mismatch)),]
	my_gaps_lower<-all_genes[which(all_genes$geneID %in% names(lower_mismatch)),]
	print(head(my_gaps_higher))
	print(head(my_gaps_lower))
	my_gaps_higher_sorted<-my_gaps_higher[match(names(higher_mismatch), my_gaps_higher$geneID),]
	my_gaps_lower_sorted<-my_gaps_lower[match(names(lower_mismatch), my_gaps_lower$geneID),]
	my_gaps_higher_wDiff<-as.data.frame(cbind(my_gaps_higher_sorted, fpkm_dif=higher_mismatch))
	#Makes these positive numbers so that chromPlot puts them facing away from the chromosome, to the left
	my_gaps_lower_wDiff<-as.data.frame(cbind(my_gaps_lower_sorted, fpkm_dif=-lower_mismatch))
	print(head(my_gaps_higher_wDiff))
	print(head(my_gaps_lower_wDiff))
	my_gaps_higher_wDiff_noNA<-my_gaps_higher_wDiff[which(!(is.na(my_gaps_higher_wDiff$start))),]
	my_gaps_lower_wDiff_noNA<-my_gaps_lower_wDiff[which(!(is.na(my_gaps_lower_wDiff$start))),]
	higher_final<-my_gaps_higher_wDiff_noNA<-my_gaps_higher_wDiff_noNA[which(my_gaps_higher_wDiff_noNA$seqnames %in% c(1:19,"X","Y")),]
	lower_final<-my_gaps_lower_wDiff_noNA<-my_gaps_lower_wDiff_noNA[which(my_gaps_lower_wDiff_noNA$seqnames %in% c(1:19,"X","Y")),]
	higher_norm<-decostand(higher_final$fpkm_dif, "range")[,1]
	lower_norm<-decostand(lower_final$fpkm_dif, "range")[,1]
	GRange_higher<-toGRanges(data.frame(chr=higher_final$seqnames, start=higher_final$start, end=higher_final$end, y=higher_norm))
	print(head(GRange_higher))
	GRange_lower<-toGRanges(data.frame(chr=lower_final$seqnames, start=lower_final$start, end=lower_final$end, y=lower_norm))
	print(head(GRange_lower))
	#max_any<-max(c(GRange_higher$y, GRange_lower$y))
	#nearest_hundred<-round(max_any+50,-2)
	kp<-plotKaryotype(genome=custom.mus.genome, chromosomes=c(1:19, "X", "Y"), plot.type=3)
	kpAddMainTitle(kp, main=paste(m,"vs mus AND dom"))
	kpAxis(kp, data.panel=1) #, ymax=nearest_hundred)
	if(length(GRange_higher) > 0){
		kpBars(kp, data=GRange_higher, y1=GRange_higher$y, col="black", data.panel=1) #, ymax=max_any
	}
	kpAxis(kp, data.panel=2) #, ymax=nearest_hundred)
	if(length(GRange_lower)){
		kpBars(kp, data=GRange_lower, y1=GRange_lower$y, col="red", data.panel=2) #, ymax=max_any
	}
}
dev.off()

print("Done with 18_chrom_plot_expression_differences.r")
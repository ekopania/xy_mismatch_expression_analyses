#PURPOSE: Many windows are only slightly negative and are probably error, not introgression; determine threshold cutoff for negative value to consider a window introgressed

win_size<-100
print(paste("Working on window size:", win_size, "kb"))
pdf(paste0("variant_count_differences_hist.",win_size,"kbWindows.pdf"), onefile=TRUE)

#Read in negative differences
mydata<-read.table(paste0("negative_windows.LEWPY.",win_size,"kbWindows.txt"), header=TRUE)
print(head(mydata))
print(dim(mydata))

#Read in original count data
toPWK_counts<-read.table(paste0("counts.LEWPY.toPWK.",win_size,"kbWindows.txt"), col.names=c("chr","start","stop","counts"))
print(dim(toPWK_counts))
print(head(toPWK_counts))
print(nrow(mydata)/nrow(toPWK_counts))
toLEWES_counts<-read.table(paste0("counts.LEWPY.toLEWES.",win_size,"kbWindows.txt"), col.names=c("chr","start","stop","counts"))
print(dim(toLEWES_counts))
print(head(toLEWES_counts))
print(nrow(mydata)/nrow(toLEWES_counts))

#Histogram of negative differences
hist(mydata$diff, main="Hist of negative values")

#Deal with some outliers
length(which(mydata$diff < -5000))
mydata_noOutlier<-mydata[which(mydata$diff > -5000),]
print(dim(mydata))
print(dim(mydata_noOutlier))
hist(mydata_noOutlier$diff, breaks=500, main="Hist of negative values, exclude < -5000")
mydata_lt500<-mydata[which(mydata$diff > -500),]
print(dim(mydata_lt500))
hist(mydata_lt500$diff, breaks=500, main="Hist of negative values, exclude < -500")
print(length(which(mydata$diff==-1)))
print(length(which(mydata$diff>-5)))
print(length(which(mydata$diff>=-4)))

#Proportion, not absolute difference (thanks, Kelsie!)
print("Proportions of negative SNPs in windows")
win_labels<-paste(mydata$chr, format(mydata$start, scientific=FALSE), sep="_")
rownames(mydata)<-win_labels
print(head(mydata))
print(tail(mydata))

LEWES_win_labels<-paste(toLEWES_counts$chr, format(toLEWES_counts$start, scientific=FALSE), sep="_")
print(length(LEWES_win_labels))
print(dim(toLEWES_counts))
rownames(toLEWES_counts)<-LEWES_win_labels
print(head(toLEWES_counts))
print(tail(toLEWES_counts))

#weird<-rownames(mydata)[which(!(rownames(mydata) %in% rownames(toLEWES_counts)))]
#print(weird)

toLEWES_counts_neg<-toLEWES_counts[which(rownames(toLEWES_counts) %in% rownames(mydata)),]
stopifnot(all.equal(rownames(toLEWES_counts_neg), rownames(mydata)))
print(head(toLEWES_counts_neg))

PWK_win_labels<-paste(toPWK_counts$chr, format(toPWK_counts$start, scientific=FALSE), sep="_")
print(length(PWK_win_labels))
print(dim(toPWK_counts))
rownames(toPWK_counts)<-PWK_win_labels

toPWK_counts_neg<-toPWK_counts[which(rownames(toPWK_counts) %in% rownames(mydata)),]
stopifnot(all.equal(rownames(toPWK_counts_neg), rownames(mydata)))
print(head(toPWK_counts_neg))

prop_toLEWES<-mydata$diff/toLEWES_counts_neg$counts
prop_toPWK<-mydata$diff/toPWK_counts_neg$counts
print(head(prop_toLEWES))
hist(prop_toLEWES, breaks=500, main="Proportion vs LEWES counts")
print(head(prop_toPWK))
hist(prop_toPWK, breaks=300, main="Proportion vs PWK counts")

stopifnot(all.equal(rownames(toPWK_counts_neg), rownames(toLEWES_counts_neg)))
avg_count<-mean(c(toPWK_counts_neg$counts, toLEWES_counts_neg$counts))
prop_toMean<-mydata$diff/avg_count
print(head(prop_toMean))
print(min(prop_toMean))
hist(prop_toMean, breaks=500, main="Proportion vs mean PWK and LEWES counts")

introgressed<-mydata_withPropToMean[which(mydata_withPropToMean$proportion < -0.1),]
dim(introgressed)
write.table(introgressed, "putative_introgressed_regions.txt", quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")

dev.off()

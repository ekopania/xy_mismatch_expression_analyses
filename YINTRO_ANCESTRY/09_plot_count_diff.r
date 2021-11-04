#PURPOSE: Get the difference in number of SNPs for each window and plot along chromosomes

#Read in counts files
sample<-"PWKLY" #"PWKLY" #"LEWPY"
ref1<-"LEWES"
ref2<-"PWK"
win_size<-10

counts1<-read.table(paste0("counts.",sample,".to",ref1,".",win_size,"kbWindows.txt"), header=FALSE, col.names=c("chr", "start", "stop", "count"))
counts2<-read.table(paste0("counts.",sample,".to",ref2,".",win_size,"kbWindows.txt"), header=FALSE, col.names=c("chr", "start", "stop", "count"))

#Loop through each window to get difference
diffs<-c()
for(i in 1:nrow(counts1)){
	this_chr<-counts1$chr[i]
	this_start<-counts1$start[i]
	this_stop<-counts1$stop[i]
	this_count1<-counts1$count[i]
	#print(paste(this_chr, this_start, this_stop, this_count1))
	if(length(intersect(which(counts2$chr==this_chr), which(counts2$start==this_start)) > 0)){
		this_count2<-counts2$count[intersect(which(counts2$chr==this_chr), which(counts2$start==this_start))]
		#print(this_count2)
		this_diff<-this_count1-this_count2
		diffs<-rbind(diffs, c(as.character(this_chr), this_start, this_stop, this_diff))
	}
	#print(dim(diffs))
}
diff_df<-as.data.frame(diffs)
#print(head(diff_df))
colnames(diff_df)<-c("chr","start","stop","diff")
diff_df$start<-as.numeric(as.character(diff_df$start))
diff_df$stop<-as.numeric(as.character(diff_df$stop))
diff_df$diff<-as.numeric(as.character(diff_df$diff))
print(dim(counts1))
print(dim(counts2))
print(dim(diff_df))
print(head(diff_df))

print("Plotting by chromosome...")
pdf(paste0("variant_count_differences_plot.",sample,".",win_size,"kbWindows.pdf"), onefile=TRUE)
chrs<-unique(diff_df$chr)
for(c in chrs){
        print(paste("Working on chromosome", c))
	temp_vars<-diff_df[which(diff_df$chr==c),]
        barplot(temp_vars$diff, names.arg=formatC(temp_vars$start, format="e", digits=2), main=paste("Difference in number of variants in",win_size,"windows along chromosome", c), xlab="Window start position", ylab="Number of variants in window")
}
dev.off()

#Extract negative (i.e., "introgressed" regions)
introgression<-diff_df[which(diff_df$diff < 0),]
print("Proportion negative windows:")
print(nrow(introgression)/nrow(diff_df))
write.table(introgression, paste0("negative_windows.",sample,".",win_size,"kbWindows.txt"), col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

print("Done with 09_plot_count_diff.r") 

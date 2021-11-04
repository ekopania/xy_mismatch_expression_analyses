#PURPOSE: Plot number of variants along X chr

#Read in data
print("Reading in variants along windows")
args<-commandArgs(TRUE)
if(length(args) != 1){
        stop("Missing command line argument. Enter the name of a file counting variants along chromosome windows (ex: counts.1kbWindows.txt)")
}

myVars<-read.table(args[1], header=FALSE, col.names=c("chr", "start", "stop", "count"))
win_size<-gsub("Windows\\.txt","",gsub("counts\\.[A-Z]*\\.to[A-Z]*\\.","",args[1]))
sample<-gsub("counts\\.","",gsub("\\.to.*","",args[1]))
ref<-gsub(".*to","",gsub("\\.[0-9].*","",args[1]))
print(paste("Plotting sample", sample, "mapped to", ref, "with window size",win_size))

#Split by chromosome
print("Plotting by chromosome...")
pdf(paste0("variant_count_plot.",sample,".to",ref,".",win_size,"Windows.pdf"), onefile=TRUE)
chrs<-unique(myVars$chr)
for(c in chrs){
	#assign(paste0("myVars_",c), myVars[which(myVars$chr==c),])
	temp_vars<-myVars[which(myVars$chr==c),]
	#plot(temp_vars$start, temp_vars$count, main=paste("Number of variants in",win_size,"windows along chromosome", c), type="l", xlab="Window start position", ylab="Number of variants in window")
	barplot(temp_vars$count, names.arg=formatC(temp_vars$start, format="e", digits=2), main=paste("Number of variants in",win_size,"windows along chromosome", c), xlab="Window start position", ylab="Number of variants in window")
}
dev.off()

#print("Plotting")
#pdf(paste("variant_count_plot",win_size,"pdf", sep="."))
#plot(myVars$start, myVars$count, main=paste("Number of variants in",win_size,"windows along X chromosome"), type="l", xlab="Window start position", ylab="Number of variants in window")
#dev.off()

print("Done with 03_plot_vars.r")

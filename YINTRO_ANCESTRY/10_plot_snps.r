#PURPOSE: Plot sample SNPs in corresponding windows when mapped to two different genomes (i.e. check if high variant windows high because of introgression or because those are just highly variable regions)
win_size<-100 #Only works for 100kb

#Read in data
print(paste("Reading in SNP counts for", win_size, "kb windows"))
wsb<-read.table(paste0("counts.toWSB.",win_size,"kbWindows.txt"))
print(dim(wsb))
print(head(wsb))
lew<-read.table(paste0("counts.toLEWES.",win_size,"kbWindows.txt"))
print(dim(lew))
print(head(lew))

#Make sure same windows and chromosomes exist for both genomes
my_chrs_wsb<-unique(wsb$V1)
my_chrs_lew<-unique(lew$V1)
print("Same chromosomes present in both genomes?") 
all.equal(my_chrs_wsb, my_chrs_lew)
for(i in my_chrs_wsb){
	print(paste("Checking lengths for chr:", i))
	print(length(which(lew$V1==i)))
	print(length(which(wsb$V1==i)))
}

#Chromosome X is slighly longer in LEWES than in WSB; last window has no SNPs so just deleting it
lew_removeExtraX<-lew[-16578,]
#These should be the same now
print(dim(lew_removeExtraX))
print(dim(wsb))

#Plot
pdf(paste0("LLLLvariants.WSBvsLEWES.", win_size,"kbWindows.pdf"))
plot(x=lew_removeExtraX$V4, y=wsb$V4, main="SNP counts in 100kb windows", xlab="# SNPs - LEWES genome", ylab="# SNPs - WSB genome")
abline(a=0, b=1)
dev.off()

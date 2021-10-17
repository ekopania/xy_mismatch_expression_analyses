#PURPOSE: Look at genes DE in different comparisons and check if they are paralogs of ampliconic gene families

library(biomaRt)

#Use biomaRt to convert between gene name and gene ID
ens_mus<-useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl")
all_genes<-getBM(attributes=c('ensembl_gene_id','external_gene_name','chromosome_name','start_position','end_position','strand','gene_biotype'), mart=ens_mus)
colnames(all_genes)<-c("gene_id","gene_name","seqnames","start","end","strand","gene_biotype")

#Read in paralogs
#sep="\t" is necessary otherwise it tries to split on " " and things with multiple gene names per blast hit will get split up, creating the wrong number of columns
ssty1<-read.table("/mnt/beegfs/ek112884/cnvs/ssty1.pID97.gene_family_paralogs.bed", sep="\t")
#Split up ones where there are 2 gene names per blast hit (e.g., Gm21440, Gm214541)
ssty1_paralogs<-unlist(sapply(ssty1$V6, function(x) strsplit(x, ", ")))
#Get rid of NAs
ssty1_paralogs_noNA<-ssty1_paralogs[which(!(is.na(ssty1_paralogs)))]
#Repeat for every Y gene family
ssty2<-read.table("/mnt/beegfs/ek112884/cnvs/ssty2.pID97.gene_family_paralogs.bed", sep="\t")
ssty2_paralogs<-unlist(sapply(ssty2$V6, function(x) strsplit(x, ", ")))
ssty2_paralogs_noNA<-ssty2_paralogs[which(!(is.na(ssty2_paralogs)))]
srsy<-read.table("/mnt/beegfs/ek112884/cnvs/srsy.pID97.gene_family_paralogs.bed", sep="\t")
srsy_paralogs<-unlist(sapply(srsy$V6, function(x) strsplit(x, ", ")))
srsy_paralogs_noNA<-srsy_paralogs[which(!(is.na(srsy_paralogs)))]
sly<-read.table("/mnt/beegfs/ek112884/cnvs/sly.pID97.gene_family_paralogs.bed", sep="\t")
sly_paralogs<-unlist(sapply(sly$V6, function(x) strsplit(x, ", ")))
sly_paralogs_noNA<-sly_paralogs[which(!(is.na(sly_paralogs)))]
asty<-read.table("/mnt/beegfs/ek112884/cnvs/asty.pID97.gene_family_paralogs.bed", sep="\t")
asty_paralogs<-unlist(sapply(asty$V6, function(x) strsplit(x, ", ")))
asty_paralogs_noNA<-asty_paralogs[which(!(is.na(asty_paralogs)))]
sstx<-read.table("/mnt/beegfs/ek112884/cnvs/sstx.pID97.gene_family_paralogs.bed", sep="\t")
sstx_paralogs<-unlist(sapply(sstx$V6, function(x) strsplit(x, ", ")))
sstx_paralogs_noNA<-sstx_paralogs[which(!(is.na(sstx_paralogs)))]
srsx<-read.table("/mnt/beegfs/ek112884/cnvs/srsx.pID97.gene_family_paralogs.bed", sep="\t")
srsx_paralogs<-unlist(sapply(srsx$V6, function(x) strsplit(x, ", ")))
srsx_paralogs_noNA<-srsx_paralogs[which(!(is.na(srsx_paralogs)))]
slx<-read.table("/mnt/beegfs/ek112884/cnvs/slx.pID97.gene_family_paralogs.bed", sep="\t")
slx_paralogs<-unlist(sapply(slx$V6, function(x) strsplit(x, ", ")))
slx_paralogs_noNA<-slx_paralogs[which(!(is.na(slx_paralogs)))]
slxl1<-read.table("/mnt/beegfs/ek112884/cnvs/slxl1.pID97.gene_family_paralogs.bed", sep="\t")
slxl1_paralogs<-unlist(sapply(slxl1$V6, function(x) strsplit(x, ", ")))
slxl1_paralogs_noNA<-slxl1_paralogs[which(!(is.na(slxl1_paralogs)))]
astx<-read.table("/mnt/beegfs/ek112884/cnvs/astx.pID97.gene_family_paralogs.bed", sep="\t")
astx_paralogs<-unlist(sapply(astx$V6, function(x) strsplit(x, ", ")))
astx_paralogs_noNA<-astx_paralogs[which(!(is.na(astx_paralogs)))]
atakusan<-read.table("/mnt/beegfs/ek112884/cnvs/atakusan.pID97.gene_family_paralogs.bed", sep="\t")
atakusan_paralogs<-unlist(sapply(atakusan$V6, function(x) strsplit(x, ", ")))
atakusan_paralogs_noNA<-atakusan_paralogs[which(!(is.na(atakusan_paralogs)))]
speer<-read.table("/mnt/beegfs/ek112884/cnvs/speer.pID97.gene_family_paralogs.bed", sep="\t")
speer_paralogs<-unlist(sapply(speer$V6, function(x) strsplit(x, ", ")))
speer_paralogs_noNA<-speer_paralogs[which(!(is.na(speer_paralogs)))]

myContrasts<-c("CCPP_RSvsCCPPLY_RS", "CCPP_RSvsWWLLPY_RS", "WWLL_RSvsWWLLPY_RS", "CCPPLY_RSvsWWLL_RS", "PPLL_RSvsPPLLPY_RS", "LLPP_RSvsPPLLPY_RS", "LLPP_RSvsLLPPLY_RS", "LLPPLY_RSvsPPLL_RS")
myAmps<-c("ssty1","ssty2","srsy","sly","asty","sstx","srsx","slx","slxl1","astx","atakusan","speer")
for(c in myContrasts){
	print(paste("Working on", c))
	mytab<-read.table(paste("topDEgenes.lrt.DE", c, "txt", sep="."), header=TRUE)
	mynames<-all_genes$gene_name[which(all_genes$gene_id %in% rownames(mytab))]
	pos<-mytab[which(mytab$direc == "+"),]
	neg<-mytab[which(mytab$direc == "-"),]
	posnames<-all_genes$gene_name[which(all_genes$gene_id %in% rownames(pos))]
	negnames<-all_genes$gene_name[which(all_genes$gene_id %in% rownames(neg))]
	for(a in myAmps){
		print(paste("Gene family:", a))
		myparalogs<-get(paste(a,"paralogs_noNA", sep="_"))
		print(length(intersect(mynames, myparalogs)))
		print(paste("Positive:", length(intersect(posnames, myparalogs))))
		print(paste("Negative:", length(intersect(negnames, myparalogs))))
	}

	#write new DE genes file with gene names added
	allnames<-c()
	for(r in rownames(mytab)){
		if(r %in% all_genes$gene_id){
			allnames<-c(allnames, all_genes$gene_name[which(all_genes$gene_id==r)])
		} else{
			allnames<-c(allnames, NA)
		}
	}
	newtab<-as.data.frame(cbind(mytab, gene_name=allnames))
	write.table(newtab, paste("topDEgenes.lrt.DE", c, "withGeneNames.txt", sep="."), row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")
}

print("Done with 16_get_amplicon_DE.r")

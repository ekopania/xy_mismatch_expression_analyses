#PURPOSE: Plot logFC values for DE genes on different chromosomes

library(biomaRt)

proco<-TRUE

#Get some gene info
ens_mus<-useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl")
all_genes<-getBM(attributes=c('ensembl_gene_id','external_gene_name','chromosome_name','start_position','end_position','strand','gene_biotype'), mart=ens_mus)
colnames(all_genes)<-c("gene_id","gene_name","seqnames","start","end","strand","gene_biotype")
proco_genes<-all_genes$gene_id[which(all_genes$gene_biotype=="protein_coding")]
#Loop through each comparison and chromosome
if(proco){
        pdf("DE_genes_logFC.barplots.protein_coding_only.pdf")
} else{
        pdf("DE_genes_logFC.barplots.pdf")
}
comparisons<-c("CCPP_RSvsWWLL_RS", "CCPP_RSvsCCPPLY_RS", "WWLL_RSvsWWLLPY_RS", "CCPP_RSvsWWLLPY_RS", "CCPPLY_RSvsWWLL_RS", "LLPP_RSvsPPLL_RS", "LLPP_RSvsLLPPLY_RS", "PPLL_RSvsPPLLPY_RS", "LLPPLY_RSvsPPLL_RS", "LLPP_RSvsPPLLPY_RS")
for(i in comparisons){
	print(i)
	mydata<-read.table(paste("topDEgenes.lrt.DE",i,"txt", sep="."), header=TRUE)
	mydata_proco<-mydata[which(rownames(mydata) %in% proco_genes), ]
	par(mfrow=c(2,2))
	mychrs<-c("5","14","X","Y")
	for(j in mychrs){
		print(j)
		mydata_chr<-mydata[which(mydata$chrs==j),]
		if(proco){
                        mydata_chr<-mydata_proco[which(mydata_proco$chrs==j),]
                } else{
                        mydata_chr<-mydata[which(mydata$chrs==j),]
                }
		if(nrow(mydata_chr) > 0){
			barplot(mydata_chr$logFC, main=paste(i, "DE genes,", j, "chr, RS"), xlab="gene", ylab="logFC")
		}
	}
}
dev.off()

#Erica did a "conservative" DE gene list - genes DE between all pairwise contrasts of PPLL vs fertile
#So to be in this conservative DE list, a gene must be DE in PPLL vs LLPP, PPLL vs CCPP, and PPLL vs WWLL
#Need to figure out what the equivalent should be for the Y introgression dataset...

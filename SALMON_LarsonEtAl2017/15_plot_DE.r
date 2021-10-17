#PURPOSE: Plot logFC values for DE genes on different chromosomes

library(biomaRt)

proco<-FALSE

#Get some gene info
ens_mus<-useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl")
all_genes<-getBM(attributes=c('ensembl_gene_id','external_gene_name','chromosome_name','start_position','end_position','strand','gene_biotype'), mart=ens_mus)
colnames(all_genes)<-c("gene_id","gene_name","seqnames","start","end","strand","gene_biotype")
proco_genes<-all_genes$gene_id[which(all_genes$gene_biotype=="protein_coding")]
#Loop through each comparison and chromosome
if(proco){
	pdf("DE_genes_logFC.barplots.LarsonEtal.protein_coding_only.pdf")
} else{
	pdf("DE_genes_logFC.barplots.LarsonEtal.pdf")
}
comparisons<-c("CCPP_RSvsWWLL_RS", "CCPP_RSvsLLPP_RS", "CCPP_RSvsPPLL_RS", "LLPP_RSvsPPLL_RS", "LLPP_RSvsWWLL_RS", "PPLL_RSvsWWLL_RS")
for(i in comparisons){
	print(i)
	mydata<-read.table(paste("topDEgenes.lrt.DE",i,"txt", sep="."), header=TRUE)
	mydata_proco<-mydata[which(rownames(mydata) %in% proco_genes), ]
	mychrs<-c("5","14","X","Y")
	par(mfrow=c(2,2))
	for(j in mychrs){
		print(j)
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
sterileVfertile<-read.table("topDEgenes.lrt.DE.LLPP_RSvsPPLL_RS.txt", header=TRUE)
sterileVfertile.up<-sterileVfertile[which(sterileVfertile$direc=="-"), ] #"-" is up in sterile hybrid
sterileVfertile.down<-sterileVfertile[which(sterileVfertile$direc=="+"), ]
sterileVmus<-read.table("topDEgenes.lrt.DE.CCPP_RSvsPPLL_RS.txt", header=TRUE)
sterileVmus.up<-sterileVmus[which(sterileVmus$direc=="-"), ] #"-" is up in sterile hybrid
sterileVmus.down<-sterileVmus[which(sterileVmus$direc=="+"), ]
sterileVdom<-read.table("topDEgenes.lrt.DE.PPLL_RSvsWWLL_RS.txt", header=TRUE)
sterileVdom.up<-sterileVdom[which(sterileVdom$direc=="+"), ] #"+" is up in sterile hybrid
sterileVdom.down<-sterileVdom[which(sterileVdom$direc=="-"), ]
DEup_genes<-Reduce(intersect, list(rownames(sterileVfertile.up), rownames(sterileVmus.up), rownames(sterileVdom.up)))
DEdown_genes<-Reduce(intersect, list(rownames(sterileVfertile.down), rownames(sterileVmus.down), rownames(sterileVdom.down)))
print(length(DEup_genes))
print(length(DEdown_genes))
print("Significant proportion of conservative DE genes upregulated in sterile hybrid?")
prop.test(length(DEup_genes), length(DEup_genes)+length(DEdown_genes))
print("Conservative DE gene set - upregulated in sterile hybrid (from sterile vs fertile hybrid):")
print(sterileVfertile[which(rownames(sterileVfertile) %in% DEup_genes), ])
print("Conservative DE gene set - downregulated in sterile hybrid (from sterile vs fertile hybrid):")
print(sterileVfertile[which(rownames(sterileVfertile) %in% DEdown_genes), ])

#DEup_table<-

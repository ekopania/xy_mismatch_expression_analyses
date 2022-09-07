#PURPOSE: Add module membership for each gene to DE genes table

#myDE<-read.csv("../SALMON/DEgenes_auto.forSupp.removeIntrogressed.txt", header=TRUE)
myDE<-read.csv("../SALMON/DEgenes_sexChr.forSupp.txt", header=TRUE)

exp1_WGCNA<-read.table("DE_modules.exp1_RS.txt", header=TRUE)
exp2_WGCNA<-read.table("DE_modules.exp2_RS.txt", header=TRUE)

myMods<-c()
for(i in 1:nrow(myDE)){
	myGene<-myDE$ensembl_gene_ID[i]
	if(myDE$comparison[i] %in% c("musdomY vs mus", "musdomY vs dom", "dommusY vs mus", "dommusY vs dom")){
		if(myGene %in% exp1_WGCNA$geneID){	
			module<-unique(exp1_WGCNA$module[which(exp1_WGCNA$geneID==myGene)])
		} else{
			module<-NA #Some genes filtered out of WGCNA
		}
	} else{
		if(myGene %in% exp2_WGCNA$geneID){
                        module<-unique(exp2_WGCNA$module[which(exp2_WGCNA$geneID==myGene)])
                } else{
                        module<-NA
                }
	}
	if(length(module) == 1){
		myMods<-c(myMods, module)
	} else{
		print(paste("ERROR: module length not 1!"))
	}
}

newDE<-as.data.frame(cbind(myDE, WGCNA_module=myMods))
#write.csv(newDE, file="DEgenes_auto.forSupp.removeIntrogressed.withMods.csv", row.names=FALSE)
write.csv(newDE, file="DEgenes_sexChr.forSupp.withMods.csv", row.names=FALSE)

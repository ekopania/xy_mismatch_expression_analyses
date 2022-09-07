#PURPOSE: Read in DE genes file and print out number of over and under expressed genes for each chromosome type

comparisons<-c("CCPP_RSvsCCPPLY_RS", "CCPPLY_RSvsWWLL_RS", "CCPP_RSvsWWLLPY_RS", "WWLL_RSvsWWLLPY_RS", "PPLL_RSvsPPLLPY_RS", "LLPPLY_RSvsPPLL_RS", "LLPP_RSvsPPLLPY_RS", "LLPP_RSvsLLPPLY_RS")

#mus vs dom and musxdom vs domxmus
#comparisons<-c("CCPP_RSvsWWLL_RS", "LLPP_RSvsPPLL_RS")

DEgenes_auto<-c()
DEgenes_sexChr<-c()

for(c in comparisons){
	print(paste("Working on comparison:", c))
	#EXCLUDING putatively introgressed regions
	#mydata<-read.table(paste("topDEgenes.lrt.DE",c,"removeIntrogressed.bed", sep="."), header=TRUE)
	#INCLUDING putatively introgressed regions
	mydata<-read.table(paste("topDEgenes.lrt.DE",c,"txt", sep="."), header=TRUE)
	if(c %in% c("CCPP_RSvsCCPPLY_RS","CCPP_RSvsWWLLPY_RS","WWLL_RSvsWWLLPY_RS","LLPPLY_RSvsPPLL")){
		over<-mydata[which(mydata$direc=="-"),]
		under<-mydata[which(mydata$direc=="+"),]
	} else{
		over<-mydata[which(mydata$direc=="+"),]
                under<-mydata[which(mydata$direc=="-"),]
	}
	print("Dimensions of all DE genes, overexpressed DE genes, and underexpressed DE genes:")
	print("'overexpressed' and 'underexpressed' are always in XY mismatch cross type relative to non-mismatch")
	print(dim(mydata))
	print(dim(over))
	print(dim(under))
	print("Numbers of DE genes:")
	print("Order is: auto overexpressed, auto underexpressed, X overexpressed, X underexpressed, Y overexpressed, Y underexpressed")
	autos_over<-over[which(!(over$chrs %in% c("X","Y"))),]
	autos_under<-under[which(!(under$chrs %in% c("X","Y"))),]
	x_over<-over[which(over$chrs=="X"), ]
	x_under<-under[which(under$chrs=="X"), ]
	y_over<-over[which(over$chrs=="Y"), ]
	y_under<-under[which(under$chrs=="Y"), ]

	print(nrow(autos_over))
	print(nrow(autos_under))
	print(nrow(x_over))
	print(nrow(x_under))
	print(nrow(y_over))
	print(nrow(y_under))
	
	myOvers_auto<-as.data.frame(cbind(rownames(autos_over), rep("overexpressed", nrow(autos_over)), rep(c, nrow(autos_over))))
	myUnders_auto<-as.data.frame(cbind(rownames(autos_under), rep("underexpressed", nrow(autos_under)), rep(c, nrow(autos_under))))
	DEgenes_auto<-as.data.frame(rbind(DEgenes_auto, myOvers_auto))
	DEgenes_auto<-as.data.frame(rbind(DEgenes_auto, myUnders_auto))

	myOvers_sexChr<-as.data.frame(cbind(c(rownames(x_over), rownames(y_over)), rep("overexpressed", nrow(x_over)+nrow(y_over)), rep(c, nrow(x_over)+nrow(y_over))))
	myUnders_sexChr<-as.data.frame(cbind(c(rownames(x_under), rownames(y_under)), rep("underexpressed", nrow(x_under)+nrow(y_under)), rep(c, nrow(x_under)+nrow(y_under))))
	DEgenes_sexChr<-as.data.frame(rbind(DEgenes_sexChr, myOvers_sexChr))
        DEgenes_sexChr<-as.data.frame(rbind(DEgenes_sexChr, myUnders_sexChr))
}

colnames(DEgenes_auto)<-c("ensembl_gene_ID", "direction", "comparison")
DEgenes_auto$comparison<-unlist(sapply(DEgenes_auto$comparison, function(x) sub("CCPP_RSvsCCPPLY_RS", "musdomY vs mus", x)))
DEgenes_auto$comparison<-unlist(sapply(DEgenes_auto$comparison, function(x) sub("CCPP_RSvsWWLLPY_RS", "dommusY vs mus", x)))
DEgenes_auto$comparison<-unlist(sapply(DEgenes_auto$comparison, function(x) sub("CCPPLY_RSvsWWLL_RS", "musdomY vs dom", x)))
DEgenes_auto$comparison<-unlist(sapply(DEgenes_auto$comparison, function(x) sub("WWLL_RSvsWWLLPY_RS", "dommusY vs dom", x)))
DEgenes_auto$comparison<-unlist(sapply(DEgenes_auto$comparison, function(x) sub("LLPP_RSvsLLPPLY_RS", "domxmus vs domxmusdomY", x)))
DEgenes_auto$comparison<-unlist(sapply(DEgenes_auto$comparison, function(x) sub("LLPP_RSvsPPLLPY_RS", "domxmus vs musxdommusY", x)))
DEgenes_auto$comparison<-unlist(sapply(DEgenes_auto$comparison, function(x) sub("LLPPLY_RSvsPPLL_RS", "musxdom vs domxmusdomY", x)))
DEgenes_auto$comparison<-unlist(sapply(DEgenes_auto$comparison, function(x) sub("PPLL_RSvsPPLLPY_RS", "musxdom vs musxdommusY", x)))

colnames(DEgenes_sexChr)<-c("ensembl_gene_ID", "direction", "comparison")
DEgenes_sexChr$comparison<-unlist(sapply(DEgenes_sexChr$comparison, function(x) sub("CCPP_RSvsCCPPLY_RS", "musdomY vs mus", x)))
DEgenes_sexChr$comparison<-unlist(sapply(DEgenes_sexChr$comparison, function(x) sub("CCPP_RSvsWWLLPY_RS", "dommusY vs mus", x)))
DEgenes_sexChr$comparison<-unlist(sapply(DEgenes_sexChr$comparison, function(x) sub("CCPPLY_RSvsWWLL_RS", "musdomY vs dom", x)))
DEgenes_sexChr$comparison<-unlist(sapply(DEgenes_sexChr$comparison, function(x) sub("WWLL_RSvsWWLLPY_RS", "dommusY vs dom", x)))
DEgenes_sexChr$comparison<-unlist(sapply(DEgenes_sexChr$comparison, function(x) sub("LLPP_RSvsLLPPLY_RS", "domxmus vs domxmusdomY", x)))
DEgenes_sexChr$comparison<-unlist(sapply(DEgenes_sexChr$comparison, function(x) sub("LLPP_RSvsPPLLPY_RS", "domxmus vs musxdommusY", x)))
DEgenes_sexChr$comparison<-unlist(sapply(DEgenes_sexChr$comparison, function(x) sub("LLPPLY_RSvsPPLL_RS", "musxdom vs domxmusdomY", x)))
DEgenes_sexChr$comparison<-unlist(sapply(DEgenes_sexChr$comparison, function(x) sub("PPLL_RSvsPPLLPY_RS", "musxdom vs musxdommusY", x)))

write.csv(DEgenes_auto, file="DEgenes_auto.forSupp.txt", row.names=FALSE)
write.csv(DEgenes_sexChr, file="DEgenes_sexChr.forSupp.txt", row.names=FALSE)

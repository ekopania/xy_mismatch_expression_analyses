#PURPOSE: Read in DE genes file and print out number of over and under expressed genes for each chromosome type


comparisons<-c("CCPP_RSvsPPLL_RS", "PPLL_RSvsWWLL_RS", "CCPP_RSvsLLPP_RS", "LLPP_RSvsWWLL_RS")

for(c in comparisons){
	print(paste("Working on comparison:", c))
	mydata<-read.table(paste("topDEgenes.lrt.DE",c,"txt", sep="."), header=TRUE)
	if(c %in% c("CCPP_RSvsPPLL_RS","CCPP_RSvsLLPP_RS")){
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
	print(length(which(!(over$chrs %in% c("X","Y")))))
	print(length(which(!(under$chrs %in% c("X","Y")))))
	print(length(which(over$chrs=="X")))
	print(length(which(under$chrs=="X")))
	print(length(which(over$chrs=="Y")))
	print(length(which(under$chrs=="Y")))
}

#PURPOSE: Log transform normalize expression data; generate table of "traits" to associate with expression modules

print("Reading in normalized expressiond data...")
#OLD - these aren't filtered by expression level
#exp1LZ<-read.table("/mnt/beegfs/ek112884/amplicon_expression_analysis/SALMON/fpkm_full_table.Yintro_exp1.LZ.95.txt", header=TRUE)
#exp1RS<-read.table("/mnt/beegfs/ek112884/amplicon_expression_analysis/SALMON/fpkm_full_table.Yintro_exp1.RS.95.txt", header=TRUE)
#exp2LZ<-read.table("/mnt/beegfs/ek112884/amplicon_expression_analysis/SALMON/fpkm_full_table.Yintro_exp2.LZ.95.txt", header=TRUE)
#exp2RS<-read.table("/mnt/beegfs/ek112884/amplicon_expression_analysis/SALMON/fpkm_full_table.Yintro_exp2.RS.95.txt", header=TRUE)
#These are filtered by expression level
exp1LZ<-read.table("/mnt/beegfs/ek112884/amplicon_expression_analysis/SALMON/fpkm_filtered_table.Yintro_exp1.LZ.97.txt", header=TRUE)
exp1RS<-read.table("/mnt/beegfs/ek112884/amplicon_expression_analysis/SALMON/fpkm_filtered_table.Yintro_exp1.RS.97.txt", header=TRUE)
exp2LZ<-read.table("/mnt/beegfs/ek112884/amplicon_expression_analysis/SALMON/fpkm_filtered_table.Yintro_exp2.LZ.97.txt", header=TRUE)
exp2RS<-read.table("/mnt/beegfs/ek112884/amplicon_expression_analysis/SALMON/fpkm_filtered_table.Yintro_exp2.RS.97.txt", header=TRUE)
#Merge, inputting NAs for genes that are expressed in some experiments/cell types but not all
all_genes<-Reduce(union, list(rownames(exp1LZ), rownames(exp1RS), rownames(exp2LZ), rownames(exp2RS)))
print(length(all_genes))
print(head(all_genes))
myExp<-c()
for(g in all_genes){
#for(g in all_genes[1:50]){ #testing...
	if(g %in% rownames(exp1LZ)){
		this_exp1LZ<-exp1LZ[g,]
	} else{
		this_exp1LZ<-rep(NA, ncol(exp1LZ))
	}
	if(g %in% rownames(exp1RS)){
                this_exp1RS<-exp1RS[g,]
        } else{
                this_exp1RS<-rep(NA, ncol(exp1RS))
        }
	if(g %in% rownames(exp2LZ)){
                this_exp2LZ<-exp2LZ[g,]
        } else{
                this_exp2LZ<-rep(NA, ncol(exp2LZ))
        }
        if(g %in% rownames(exp2RS)){
                this_exp2RS<-exp2RS[g,]
        } else{
                this_exp2RS<-rep(NA, ncol(exp2RS))
        }
	myExp<-rbind(myExp, c(this_exp1LZ, this_exp1RS, this_exp2LZ, this_exp2RS))
	print(dim(myExp))
}
print(myExp[1:5,1:5])
#myExp<-as.data.frame(myExp)
rownames(myExp)<-all_genes
#rownames(myExp)<-all_genes[1:50] #testing....
colnames(myExp)<-c(colnames(exp1LZ), colnames(exp1RS), colnames(exp2LZ), colnames(exp2RS))
#myExp<-read.table("fpkm_full_table.txt", header=TRUE)
#myExp<-as.data.frame(cbind(exp1LZ, exp1RS, exp2LZ, exp2RS))
print(dim(myExp))
write.table(myExp, file="fpkm_full_table.txt", quote=FALSE, sep="\t")

print("Performing log transformation: log2(x+1)...")
myLogData<-apply(myExp, c(1,2), function(x) log2(as.numeric(x)+1))
print(dim(myLogData))
write.table(myLogData, file="logTransformed_expData.txt", quote=FALSE, sep="\t")

#print("Generating trait table...")
#samples<-colnames(myExp)
#cross_type<-sapply(samples, function(x) substr(x, 1, 4))
#ct<-sapply(samples, function(x) substr(x, nchar(x)-1, nchar(x)))
#yMismatch<-c()
#XYdirection<-c()
#for(i in samples){
#	#Test for XY mismatch in either direction
#	if(grepl("CPLY|WLPY|LLPP|PPLL", i)){
#		yMismatch<-c(yMismatch,"yes")
#	} else{
#		yMismatch<-c(yMismatch,"no")
#	}
#	#Test for too much X (mus X and dom Y), too little X (dom X and mus Y), or matching XY
#	if(grepl("CPLY|PPLL", i)){
#		XYdirection<-c(XYdirection, "musXdomY")
#	} else if(grepl("WLPY|LLPP", i)){
#		XYdirection<-c(XYdirection, "domXmusY")
#	} else{
#		XYdirection<-c(XYdirection, "XYmatch")
#	}
#}
#
#mytraits<-as.data.frame(cbind(cross_type=cross_type, cell_type=ct, XYmismatch=yMismatch, XYdirection=XYdirection))
#rownames(mytraits)<-samples
#write.table(mytraits, file="traits.txt", quote=FALSE, sep="\t")

print("Done with 01_process_exp_data.r")

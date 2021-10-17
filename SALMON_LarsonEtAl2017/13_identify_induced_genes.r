#PURPOSE: Identify genes induced in a cell type
#INDUCED: (median expression in one cell type) > (2*median expression in other cell type)

#change these to test different parameters
min_rpkm<-1 #1
min_samples<-4 #4 #16
#min_logFC<-0
induced_cutoff<-2

#Using Larson et al. 2017 cell sorted F1s as a test at first; will eventually use the Y introgression data
dataset<-"LarsonEtal"

print("Loading R libraries...")
library(edgeR)
library(ggplot2)
library(gridExtra)
library(EnsDb.Mmusculus.v79)
library(biomaRt)
library(matrixStats)

print("Setting up DGE object...")
mydata<-read.table("salmon_output.counts.txt")
mylengths<-read.table("salmon_output.lengths.txt")
#Groups: 1=CCPP_LZ, 2=CCPP_RS, 3=WWLL_LZ 4=WWLL_RS, 5=LLPP_LZ, 6=LLPP_RS, 7=PPLL_LZ, 8=PPLL_RS
if(dataset=="LarsonEtal"){ #Always include both LZ and RS because we're trying to compare the two
	these_data<-mydata
	these_lengths<-mylengths
	mygroups<-c(1,2,1,2,1,2,6,5,6,6,5,5,8,7,8,7,8,7,4,3,4,3,4,3)
	geno_order<-c("CCPP_LZ","CCPP_RS","WWLL_LZ","WWLL_RS","LLPP_LZ","LLPP_RS","PPLL_LZ","PPLL_RS")
} else{
        stop("Invalid dataset")
}

#NEED TO ACCOUNT FOR GENE LENGTHS - they may be different among samples because we originally came from txpt level counts
#SEE TXIMPORT DOCUMENTATION: http://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html#Downstream_DGE_in_Bioconductor
normMat<-these_lengths/exp(rowMeans(log(these_lengths)))
normCts<-these_data/normMat
eff.lib <- calcNormFactors(normCts) * colSums(normCts)
normMat <- sweep(normMat, 2, eff.lib, "*")
normMat <- log(normMat)
num_normMat<-data.matrix(normMat) #Convert data frame to numeric matrix to keep EdgeR happy
print(dim(these_data))
print(length(mygroups))
myDGE<-DGEList(counts=these_data, group=mygroups)
myDGE<-scaleOffset(myDGE, num_normMat)
print("Dimensions of myDGE counts:")
print(dim(myDGE$counts))
print("myDGE samples - CHECK GROUP LABELS:")
print(myDGE$samples)

print("Filtering by expression level...")
#EdgeR rpkm function seems to require same gene length across all samples; taking median of tximport output gene lengths and using that
new_lengths<-rowMedians(as.matrix(these_lengths))
names(new_lengths)<-rownames(these_lengths)
print(head(new_lengths))
my_rpkm<-rpkm(myDGE, gene.length=new_lengths)

keep<-rowSums(my_rpkm > min_rpkm) >= min_samples
keep[is.na(keep)]<-FALSE #NAs show up for gene IDs that aren't in myDGE; replace then with false otherwise myDGE gets confused

myDGE<-myDGE[keep, , keep.lib.sizes=FALSE]
myDGE<-calcNormFactors(myDGE)
print(myDGE$samples)

print("Dimensions of filtered myDGE counts:")
print(dim(myDGE$counts))

#Get new fpkm table based on filtered myDGE counts
newOrdered_lengths<-new_lengths[match(rownames(myDGE$counts),names(new_lengths))] #filtered myDGE order
new_rpkm<-rpkm(myDGE,gene.length=newOrdered_lengths) #Based on filtered gene set
print(dim(new_rpkm))


#NOTE: When I looked for induced genes I first identified genes that were expressed in a cell type in all 4 (sub)species; I don't think that's necessary here because I'm only dealing with dom and mus so most genes should be expressed in both if they're expressed at all according to my original DGE filtering

#Get cell type specific fpkm tables
#FPKM<-fpkm
#PHYLO DATA - two cell types
#LZ<-c(1,3,5,7)
#RS<-c(2,4,6,8)
#PHYLO DATA - NO PAHARI, two cell types,
#LZ<-c(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,42,44,46,48,50,52,54,56,58)
#RS<-c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,41,43,45,47,49,51,53,55,57,59)
FPKM_LZ<-new_rpkm[,grepl("LZ",colnames(new_rpkm))]
FPKM_RS<-new_rpkm[,grepl("RS",colnames(new_rpkm))]

#determine genes induced (not just expressed) in each cell type
#induced: median FPKM in one cell type > induced_cutoff*median in other cell type
#induced_cutoff set at beginning of script
LZ_meds<-apply(FPKM_LZ,1,median)
if(length(LZ_meds) != nrow(new_rpkm)){
    print("ERROR: number of LZ medians not equal to number of genes in analysis")
}
LZ_meds<-as.table(LZ_meds)
RS_meds<-apply(FPKM_RS,1,median)
if(length(RS_meds) != nrow(new_rpkm)){
    print("ERROR: number of RS medians not equal to number of genes in analysis")
}
RS_meds<-as.table(RS_meds)
#Two cell types
LZ_induced<-c()
for(i in rownames(LZ_meds)){
    if(i %in% rownames(RS_meds)){
        if((LZ_meds[i]>induced_cutoff*RS_meds[i]) & (LZ_meds[i]>1)){ # & (min(FPKM_LZ[i,])>1)){
            LZ_induced<-c(LZ_induced,i)
        }
    }
    else if((!(i %in% rownames(RS_meds)))  & (LZ_meds[i]>1)){ # & (min(FPKM_LZ[i,])>1)){
        LZ_induced<-c(LZ_induced,i)
    }
}
RS_induced<-c()
for(i in rownames(RS_meds)){
    if((i %in% rownames(LZ_meds)) & (RS_meds[i]>1)){ # & (min(FPKM_RS[i,])>1)){
        if(RS_meds[i]>induced_cutoff*LZ_meds[i]){
            RS_induced<-c(RS_induced,i)
        }
    }
    else if((!(i %in% rownames(LZ_meds))) & (RS_meds[i]>1)){ # & (min(FPKM_RS[i,])>1)){
        RS_induced<-c(RS_induced,i)
    }
}
print("Here is length and head for LZ and RS:")
print(length(LZ_induced))
print(length(RS_induced))
print(head(LZ_induced))
print(head(RS_induced))
for(i in LZ_induced){
    if(i %in% RS_induced){
        print(paste("ERROR: ", i, "cannot be induced in both LZ and RS"))
    }
}
for(i in RS_induced){
    if(i %in% LZ_induced){
        print(paste("ERROR: ", i, "cannot be induced in both LZ and RS"))
    }
}
#WITH (min(FPKM_RS[i,])>1)) cutoff
#write(LZ_induced, file=paste("gene_list_LZinduced_edgeR",dataset,"txt",sep="."),append=FALSE)
#write(RS_induced, file=paste("gene_list_RSinduced_edgeR",dataset,"txt",sep="."),append=FALSE)
#WITHOUT (min(FPKM_RS[i,])>1))
write(LZ_induced, file=paste("gene_list_LZinduced_edgeR.lessStrict",dataset,"txt",sep="."),append=FALSE)
write(RS_induced, file=paste("gene_list_RSinduced_edgeR.lessStrict",dataset,"txt",sep="."),append=FALSE)

print("Done with 13_identify_induced_genes.r") 

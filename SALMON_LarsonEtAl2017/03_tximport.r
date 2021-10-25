#PURPOSE: Convert from transcript level counts to gene level counts using tximport

#Read in libraries
library(tximport)
library(rhdf5)
library(EnsDb.Mmusculus.v79)

print("Reading in salmon output *.sf files...")
myfiles<-list.files(path=Sys.glob("*salmon_output_libTypeIU_premadeIndex"), pattern="quant.sf", recursive=TRUE, full.names=TRUE)
print(myfiles)
sample_names<-c("CCPP211MLZ" ,"CCPP211MRS" ,"CCPP212MLZ" ,"CCPP212MRS" ,"CCPP213MLZ" ,"CCPP213MRS" ,"LLPP172MRS" ,"LLPP181MLZ" ,"LLPP192MRS","LLPP193MRS","LLPP227MLZ","LLPP228MLZ","PPLL152MRS","PPLL161MLZ","PPLL161MRS","PPLL171MLZ","PPLL171MRS","PPLL173MLZ","WWLL31MRS","WWLL41MLZ","WWLL61MRS","WWLL72MLZ","WWLL72MRS","WWLL73MLZ")
names(myfiles)<-sample_names

edb<-EnsDb.Mmusculus.v79
mytxpts<-transcripts(edb, return.type="DataFrame")
mytx2gene<-as.data.frame(cbind(TXNAME=mytxpts$tx_id,GENEID=mytxpts$gene_id))

print("Running tximport...")
txi.salmon<-tximport(myfiles, type="salmon", tx2gene=mytx2gene, ignoreTxVersion=TRUE)
print(dim(txi.salmon$counts))
print(head(txi.salmon$counts))

print("Writing to output file...")
write.table(txi.salmon$counts, file="salmon_output.counts.txt", quote=FALSE, sep="\t")
write.table(txi.salmon$length, file="salmon_output.lengths.txt", quote=FALSE, sep="\t")

print("Done with 03_tximport.r")

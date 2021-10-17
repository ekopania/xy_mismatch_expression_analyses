#PURPOSE: Convert from transcript level counts to gene level counts using tximport

#Read in libraries
library(tximport)
library(rhdf5)
library(EnsDb.Mmusculus.v79)

print("Reading in salmon output *.h5 files...")
myfiles<-list.files(pattern="quant.sf",recursive=TRUE)
print(myfiles)
sample_names<-c("CCPP451MLZ","CCPP451MRS","CCPP452MLZ","CCPP452MRS","CCPP453MLZ","CCPP453MRS","CCPP454MLZ","CCPP454MRS","CPLY186MLZ","CPLY186MRS","CPLY213MLZ","CPLY213MRS","CPLY214MLZ","CPLY214MRS","CPLY215MLZ","CPLY215MRS","LLPP292MLZ","LLPP292MRS","LLPP293MLZ","LLPP293MRS","LLPP294MLZ","LLPP294MRS","LLPP306MLZ","LLPP306MRS","LPLY15MLZ","LPLY15MRS","LPLY16MLZ","LPLY16MRS","LPLY17MLZ","LPLY17MRS","LPLY35MLZ","LPLY35MRS","PPLL306MLZ","PPLL306MRS","PPLL307MLZ","PPLL307MRS","PPLL326MLZ","PPLL326MRS","PPLL327MLZ","PPLL327MRS","PLPY23MLZ","PLPY23MRS","PLPY26MLZ","PLPY26MRS","PLPY32MLZ","PLPY32MRS","PLPY46MLZ","PLPY46MRS","WWLL505MLZ","WWLL505MRS","WWLL506MLZ","WWLL506MRS","WWLL531MLZ","WWLL531MRS","WWLL532MLZ","WWLL532MRS","WLPY105MLZ","WLPY105MRS","WLPY125MLZ","WLPY125MRS","WLPY154MLZ","WLPY154MRS","WLPY94MLZ","WLPY94MRS")
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

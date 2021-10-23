#PURPOSE: Plot logFC values as violin plots and compare with Wilcoxon rank sum test

library(ggplot2)
library(ggbeeswarm)

dataset<-"Yintro_exp2"
cell_type<-"LZ"

#Read in data
fpkm_data<-read.table(paste("fpkm_filtered_table",dataset,cell_type,"97.txt", sep="."), header=TRUE)

if(dataset=="Yintro_exp1"){
        #Set up with sex chr mismatch always first
        cross_types<-list(c("CPLY","CCPP"), c("WLPY","CCPP"), c("WLPY","WWLL"), c("CPLY","WWLL"))
        names(cross_types)<-c("musX","musY","domX","domY")
        if(cell_type=="LZ"){
                musX<-read.table("topDEgenes.lrt.DE.CCPP_LZvsCCPPLY_LZ.txt", header=TRUE, fill=TRUE)
                musY<-read.table("topDEgenes.lrt.DE.CCPP_LZvsWWLLPY_LZ.txt", header=TRUE, fill=TRUE)
                domX<-read.table("topDEgenes.lrt.DE.WWLL_LZvsWWLLPY_LZ.txt", header=TRUE, fill=TRUE)
                domY<-read.table("topDEgenes.lrt.DE.CCPPLY_LZvsWWLL_LZ.txt", header=TRUE, fill=TRUE)
        } else if(cell_type=="RS"){
                musX<-read.table("topDEgenes.lrt.DE.CCPP_RSvsCCPPLY_RS.withGeneNames.txt", header=TRUE, fill=TRUE)
                musY<-read.table("topDEgenes.lrt.DE.CCPP_RSvsWWLLPY_RS.withGeneNames.txt", header=TRUE, fill=TRUE)
                domX<-read.table("topDEgenes.lrt.DE.WWLL_RSvsWWLLPY_RS.withGeneNames.txt", header=TRUE, fill=TRUE)
                domY<-read.table("topDEgenes.lrt.DE.CCPPLY_RSvsWWLL_RS.withGeneNames.txt", header=TRUE, fill=TRUE)
        } else{
                stop("Invalid cell type: must be 'LZ' or 'RS'")
        }
} else if(dataset=="Yintro_exp2"){
        cross_types<-list(c("PPLL","PLPY"), c("LLPP","PLPY"), c("LLPP","LPLY"), c("PPLL","LPLY"))
        names(cross_types)<-c("musX","musY","domX","domY")
        if(cell_type=="LZ"){
                musX<-read.table("topDEgenes.lrt.DE.PPLL_LZvsPPLLPY_LZ.txt", header=TRUE, fill=TRUE)
                musY<-read.table("topDEgenes.lrt.DE.LLPP_LZvsPPLLPY_LZ.txt", header=TRUE, fill=TRUE)
                domX<-read.table("topDEgenes.lrt.DE.LLPP_LZvsLLPPLY_LZ.txt", header=TRUE, fill=TRUE)
                domY<-read.table("topDEgenes.lrt.DE.LLPPLY_LZvsPPLL_LZ.txt", header=TRUE, fill=TRUE)
        } else if(cell_type=="RS"){
                musX<-read.table("topDEgenes.lrt.DE.PPLL_RSvsPPLLPY_RS.withGeneNames.txt", header=TRUE, fill=TRUE)
                musY<-read.table("topDEgenes.lrt.DE.LLPP_RSvsPPLLPY_RS.withGeneNames.txt", header=TRUE, fill=TRUE)
                domX<-read.table("topDEgenes.lrt.DE.LLPP_RSvsLLPPLY_RS.withGeneNames.txt", header=TRUE, fill=TRUE)
                domY<-read.table("topDEgenes.lrt.DE.LLPPLY_RSvsPPLL_RS.withGeneNames.txt", header=TRUE, fill=TRUE)
        } else{
                stop("Invalid cell type: must be 'LZ' or 'RS'")
        }
} else{
        stop("Invalid dataset: must be 'Yintro_exp1' or 'Yintro_exp2'")
}

print("Here are the cross type comparisons:")
print(cross_types)

print("Appending chromosome type...")
chr_type_musX<-sapply(musX$chr, function(x) if(x %in% c(1:19)){"auto"} else{x})
#print(head(chr_type_musX))
#print(tail(chr_type_musX))
musX_wChr<-as.data.frame(cbind(musX, chr_type=chr_type_musX))
musX_df<-musX_wChr[which(musX_wChr$chr_type %in% c("auto","X","Y")),]
#print(head(musX_df))
chr_type_musY<-sapply(musY$chr, function(x) if(x %in% c(1:19)){"auto"} else{x})
musY_wChr<-as.data.frame(cbind(musY, chr_type=chr_type_musY))
musY_df<-musY_wChr[which(musY_wChr$chr_type %in% c("auto","X","Y")),]
chr_type_domX<-sapply(domX$chr, function(x) if(x %in% c(1:19)){"auto"} else{x})
domX_wChr<-as.data.frame(cbind(domX, chr_type=chr_type_domX))
domX_df<-domX_wChr[which(domX_wChr$chr_type %in% c("auto","X","Y")),]
chr_type_domY<-sapply(domY$chr, function(x) if(x %in% c(1:19)){"auto"} else{x})
domY_wChr<-as.data.frame(cbind(domY, chr_type=chr_type_domY))
domY_df<-domY_wChr[which(domY_wChr$chr_type %in% c("auto","X","Y")),]

print("Plotting...")
pdf(paste("logFC_violins",dataset,cell_type,"pdf", sep="."), onefile=TRUE)
for(c in names(cross_types)){
	print(paste("Working on comparison for:", c))
	this_df<-get(paste(c,"df", sep="_"))
	print("Wilcoxon test:")
	print(pairwise.wilcox.test(this_df$logFC, this_df$chr_type, method="fdr"))
	#Might need to do negative logFC such that positive numbers on Y-axis always correspond to higher expression in XY mismatch crosses
	if( (dataset=="Yintro_exp1") && (c %in% c("musX","musY","domX"))){
		p<-ggplot(this_df, aes(x=chr_type, y=-logFC, color=chr_type)) + geom_violin() + geom_boxplot(width=0.1) + geom_quasirandom()
	} else if( (dataset=="Yintro_exp2") && (c=="domY") ){
		p<-ggplot(this_df, aes(x=chr_type, y=-logFC, color=chr_type)) + geom_violin() + geom_boxplot(width=0.1) + geom_quasirandom()
	} else{
		p<-ggplot(this_df, aes(x=chr_type, y=logFC, color=chr_type)) + geom_violin() + geom_boxplot(width=0.1) + geom_quasirandom()
	}
	p<-p + labs(title=paste("LogFC by chromosome type:", cross_types[[c]][1], "vs", cross_types[[c]][2]), x="", y="logFC")
	p<-p + theme(axis.text=element_text(size=18), axis.title=element_text(size=21), plot.title=element_text(size=28))
	p<-p + theme_minimal() + ylim(-15,15) + scale_color_manual(values=c("lightgrey","darkgrey","black"))
	print(p)
}
dev.off()

print("Done with 20_logFC_violin_plots.r")
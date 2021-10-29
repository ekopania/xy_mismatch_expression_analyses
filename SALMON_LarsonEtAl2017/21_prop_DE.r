#PURPOSE: Generate barplots of the proportion of expressed genes that are DE between different cross type comparisons

library(biomaRt)

dataset<-"LarsonEtal"
cell_type<-"RS"
procoOnly<-FALSE

#Read in data
fpkm_data<-read.table(paste("fpkm_filtered_table",dataset,cell_type,"97.txt", sep="."), header=TRUE)

#Set up with hybrid always first
cross_types<-list(c("PPLL","CCPP"), c("LLPP","CCPP"), c("LLPP","WWLL"), c("PPLL","WWLL"))
names(cross_types)<-c("musX","musY","domX","domY")
if(cell_type=="LZ"){
        musX<-read.table("topDEgenes.lrt.DE.CCPP_LZvsPPLL_LZ.txt", header=TRUE, fill=TRUE)
        musY<-read.table("topDEgenes.lrt.DE.CCPP_LZvsLLPP_LZ.txt", header=TRUE, fill=TRUE)
        domX<-read.table("topDEgenes.lrt.DE.LLPP_LZvsWWLL_LZ.txt", header=TRUE, fill=TRUE)
        domY<-read.table("topDEgenes.lrt.DE.PPLL_LZvsWWLL_LZ.txt", header=TRUE, fill=TRUE)
} else if(cell_type=="RS"){
        musX<-read.table("topDEgenes.lrt.DE.CCPP_RSvsPPLL_RS.txt", header=TRUE, fill=TRUE)
        musY<-read.table("topDEgenes.lrt.DE.CCPP_RSvsLLPP_RS.txt", header=TRUE, fill=TRUE)
        domX<-read.table("topDEgenes.lrt.DE.LLPP_RSvsWWLL_RS.txt", header=TRUE, fill=TRUE)
        domY<-read.table("topDEgenes.lrt.DE.PPLL_RSvsWWLL_RS.txt", header=TRUE, fill=TRUE)
} else{
        stop("Invalid cell type: must be 'LZ' or 'RS'")
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

#Separate expressed genes by chr type
ens_mus<-useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl")
all_genes<-getBM(attributes=c('ensembl_gene_id','chromosome_name','start_position','end_position'), mart=ens_mus)
colnames(all_genes)<-c("geneID", "chr", "start", "end")
print(head(all_genes))
auto_genes<-all_genes$geneID[which(all_genes$chr %in% c(1:19))]
x_genes<-all_genes$geneID[which(all_genes$chr=="X")]
y_genes<-all_genes$geneID[which(all_genes$chr=="Y")]

fpkm_auto<-fpkm_data[which(rownames(fpkm_data) %in% auto_genes),]
fpkm_x<-fpkm_data[which(rownames(fpkm_data) %in% x_genes),]
fpkm_y<-fpkm_data[which(rownames(fpkm_data) %in% y_genes),]

if(procoOnly){
        pdf(paste("prop_DE_genes.barplots",dataset,cell_type,"procoOnly.pdf", sep="."))
} else{
        pdf(paste("prop_DE_genes.barplots",dataset,cell_type,"pdf", sep="."))
}
#Loop through each comparison
for(c in names(cross_types)){
        print(paste("Working on comparison for:", c))
        this_df<-get(paste(c,"df", sep="_"))
	#Get expression levels for samples in this comparison
	fpkm_auto_1<-fpkm_auto[, grepl(cross_types[[c]][1], colnames(fpkm_auto))]
	fpkm_auto_2<-fpkm_auto[, grepl(cross_types[[c]][2], colnames(fpkm_auto))]
	fpkm_x_1<-fpkm_x[, grepl(cross_types[[c]][1], colnames(fpkm_x))]
        fpkm_x_2<-fpkm_x[, grepl(cross_types[[c]][2], colnames(fpkm_x))]
	fpkm_y_1<-fpkm_y[, grepl(cross_types[[c]][1], colnames(fpkm_y))]
        fpkm_y_2<-fpkm_y[, grepl(cross_types[[c]][2], colnames(fpkm_y))]
	#Get total expressed genes for each chromosome type, separate by higher in hybrid or lower in hybrid
	auto_higher<-c()
	auto_lower<-c()
	for(i in 1:nrow(fpkm_auto_1)){
		if(median(as.numeric(as.character(fpkm_auto_1[i,]))) > median(as.numeric(as.character(fpkm_auto_2[i,])))){
			auto_higher<-c(auto_higher, rownames(fpkm_auto_1)[i])
		} else{
			auto_lower<-c(auto_lower, rownames(fpkm_auto_1)[i])
		}
	}
	x_higher<-c()
        x_lower<-c()
        for(i in 1:nrow(fpkm_x_1)){
                if(median(as.numeric(as.character(fpkm_x_1[i,]))) > median(as.numeric(as.character(fpkm_x_2[i,])))){
                        x_higher<-c(x_higher, rownames(fpkm_x_1)[i])
                } else{
                        x_lower<-c(x_lower, rownames(fpkm_x_1)[i])
                }
        }
	y_higher<-c()
        y_lower<-c()
        for(i in 1:nrow(fpkm_y_1)){
                if(median(as.numeric(as.character(fpkm_y_1[i,]))) > median(as.numeric(as.character(fpkm_y_2[i,])))){
                        y_higher<-c(y_higher, rownames(fpkm_y_1)[i])
                } else{
                        y_lower<-c(y_lower, rownames(fpkm_y_1)[i])
                }
        }
	#Separate DE genes by chromosome type and by higher or lower in hybrid
	if(c %in% c("musX","musY")){
        	DE_auto_higher0<-this_df[intersect(which(this_df$chrs %in% c(1:19)), which(this_df$direc=="-")),]
		DE_auto_lower0<-this_df[intersect(which(this_df$chrs %in% c(1:19)), which(this_df$direc=="+")),]
		DE_x_higher0<-this_df[intersect(which(this_df$chrs=="X"), which(this_df$direc=="-")),]
		DE_x_lower0<-this_df[intersect(which(this_df$chrs=="X"), which(this_df$direc=="+")),]
		DE_y_higher0<-this_df[intersect(which(this_df$chrs=="Y"), which(this_df$direc=="-")),]
                DE_y_lower0<-this_df[intersect(which(this_df$chrs=="Y"), which(this_df$direc=="+")),]
	} else{
		DE_auto_higher0<-this_df[intersect(which(this_df$chrs %in% c(1:19)), which(this_df$direc=="+")),]
                DE_auto_lower0<-this_df[intersect(which(this_df$chrs %in% c(1:19)), which(this_df$direc=="-")),]
                DE_x_higher0<-this_df[intersect(which(this_df$chrs=="X"), which(this_df$direc=="+")),]
                DE_x_lower0<-this_df[intersect(which(this_df$chrs=="X"), which(this_df$direc=="-")),]
                DE_y_higher0<-this_df[intersect(which(this_df$chrs=="Y"), which(this_df$direc=="+")),]
                DE_y_lower0<-this_df[intersect(which(this_df$chrs=="Y"), which(this_df$direc=="-")),]
	}
	#Check that DE genes for each direction are in the right higher and lower groups
	DE_auto_higher<-DE_auto_higher0[which(rownames(DE_auto_higher0) %in% auto_higher),]
	DE_auto_lower<-DE_auto_lower0[which(rownames(DE_auto_lower0) %in% auto_lower),]
	DE_x_higher<-DE_x_higher0[which(rownames(DE_x_higher0) %in% x_higher),]
        DE_x_lower<-DE_x_lower0[which(rownames(DE_x_lower0) %in% x_lower),]
	DE_y_higher<-DE_y_higher0[which(rownames(DE_y_higher0) %in% y_higher),]
        DE_y_lower<-DE_y_lower0[which(rownames(DE_y_lower0) %in% y_lower),]
	#Get proportion of expressed genes that are DE for each category
	prop_DE_higher<-c(nrow(DE_auto_higher)/length(auto_higher), nrow(DE_x_higher)/length(x_higher), nrow(DE_y_higher)/length(y_higher))
	prop_DE_lower<-c(nrow(DE_auto_lower)/length(auto_lower), nrow(DE_x_lower)/length(x_lower), nrow(DE_y_lower)/length(y_lower))
	names(prop_DE_higher)<-c("auto","X","Y")
        names(prop_DE_lower)<-c("auto","X","Y")
	#Plot
	barplot(prop_DE_higher, main=paste("Proportion DE genes:", cross_types[[c]][1], "vs", cross_types[[c]][2]), ylim=c(-1,1))
        barplot(-prop_DE_lower, add = TRUE)
	#Pearson's chi-squared tests
	print("X vs autos:")
        print(prop.test(c(nrow(DE_x_higher)+nrow(DE_x_lower), nrow(DE_auto_higher)+nrow(DE_auto_lower)), c(length(x_higher)+length(x_lower), length(auto_higher)+length(auto_lower))))
        print("Y vs autos:")
        print(prop.test(c(nrow(DE_y_higher)+nrow(DE_y_lower), nrow(DE_auto_higher)+nrow(DE_auto_lower)), c(length(y_higher)+length(y_lower), length(auto_higher)+length(auto_lower))))
	if((nrow(DE_auto_higher) > 0) || (nrow(DE_auto_lower) > 0) ){
                print("Up vs down, autos:")
                print(prop.test(c(nrow(DE_auto_higher), nrow(DE_auto_lower)), c(nrow(DE_auto_higher)+nrow(DE_auto_lower), nrow(DE_auto_higher)+nrow(DE_auto_lower))))
        }
	if((nrow(DE_x_higher) > 0) || (nrow(DE_x_lower) > 0) ){
                print("Up vs down, x:")
                print(prop.test(c(nrow(DE_x_higher), nrow(DE_x_lower)), c(nrow(DE_x_higher)+nrow(DE_x_lower), nrow(DE_x_higher)+nrow(DE_x_lower))))
        }
	if((nrow(DE_y_higher) > 0) || (nrow(DE_y_lower) > 0) ){
                print("Up vs down, y:")
                print(prop.test(c(nrow(DE_y_higher), nrow(DE_y_lower)), c(nrow(DE_y_higher)+nrow(DE_y_lower), nrow(DE_y_higher)+nrow(DE_y_lower))))
        }
}
dev.off()

print("Done with 21_prop_DE.r") 

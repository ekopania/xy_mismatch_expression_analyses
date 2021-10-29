#PURPOSE: Use EdgeR to identify genes most DE between Y intro and controls (exp 1) or between F1 hybrid sex chromosome mismatch and match (exp2)

#change these to test different parameters
#Normally I'd filter out low exp things but given that we want to catch everything that's part of these gene families it might be best to leave everything zero for now; one individual paralog might have a low expression level that looks like noise but it still might be contributing to the overall expression of the gene family
#Also, some of these gene families (srsx I believe) have FPKM <1 for every paralog across all samples
min_rpkm<-1 #1
min_samples<-4 #4 #16
#min_logFC<-0
#induced_cutoff<-2

#Choose which dataset and cell type to analyze
dataset<-"Yintro_exp1"
#dataset<-"Yintro_exp2"
#cell_type<-"LZandRS"
#cell_type<-"LZ"
cell_type<-"RS"
procoOnly<-FALSE

print(paste("Running EdgeR DE analysis for: ", dataset, cell_type))
print(paste("Protein coding only?", procoOnly))
print(paste("Cutoff of RPKM >", min_rpkm, "in", min_samples, "samples"))

print("Loading R libraries...")
library(edgeR)
#library(EnsDb.Mmusculus.v79)
library(ggplot2)
library(gridExtra)
library(biomaRt)
library(matrixStats)

print("Setting up DGE object...")
#OLD - No decoys for salmon mapping
#mydata<-read.table("salmon_output.counts.txt")
#mylengths<-read.table("salmon_output.lengths.txt")
#New version with Salmon's premade mm10 index that includes decoys
mydata<-read.table("salmon_output.counts.premadeIndex.txt")
mylengths<-read.table("salmon_output.lengths.premadeIndex.txt")
#For Y intro data, 1=CCPP_LZ, 2=CCPP_RS, 3=CCPPLY_LZ, 4=CCPPLY_RS, 5=LLPP_LZ, 6=LLPP_RS, 7=LLPPLY_LZ, 8=LLPPLY_RS, 9=PPLLPY_LZ, 10=PPLLPY_RS, 11=PPLL_LZ, 12=PPLL_RS, 13=WWLLPY_LZ, 14=WWLLPY_RS, 15=WWLL_LZ, 16=WWLL_RS
#For appropriate comparisons in Erica's dataset (2017 MBE paper)
contrast_to_larson<-c("CCPP_RSvsPPLL_RS", "PPLL_RSvsWWLL_RS", "LLPP_RSvsWWLL_RS", "CCPP_RSvsLLPP_RS")
if(dataset=="Yintro_exp1"){
        if(cell_type=="LZandRS"){
                these_data<-mydata[,grepl("CP|WL",colnames(mydata))]
                these_lengths<-mylengths[,grepl("CP|WL",colnames(mydata))]
                mygroups<-c(1,2,1,2,1,2,1,2,3,4,3,4,3,4,3,4,13,14,13,14,13,14,13,14,15,16,15,16,15,16,15,16)
                geno_order<-c("CCPP_LZ","CCPP_RS","CPLY_LZ","CPLY_RS","WWLL_LZ","WWLL_RS","WLPY_LZ","WLPY_RS")
        } else if(cell_type=="LZ"){
                these_data<-mydata[,grepl("CP.*LZ|WL.*LZ",colnames(mydata))]
                these_lengths<-mylengths[,grepl("CP.*LZ|WL.*LZ",colnames(mylengths))]
                mygroups<-c(1,1,1,1,3,3,3,3,13,13,13,13,15,15,15,15) #LZ only
                geno_order<-c("CCPP_LZ","CPLY_LZ","WWLL_LZ","WLPY_LZ")
		names(contrast_to_larson)<-c("CCPP_LZvsCCPPLY_LZ", "CCPPLY_LZvsWWLL_LZ", "WWLL_LZvsWWLLPY_LZ", "CCPP_LZvsWWLLPY_LZ")
        } else if(cell_type=="RS"){
                these_data<-mydata[,grepl("CP.*RS|WL.*RS",colnames(mydata))]
                these_lengths<-mylengths[,grepl("CP.*RS|WL.*RS",colnames(mylengths))]
                mygroups<-c(2,2,2,2,4,4,4,4,14,14,14,14,16,16,16,16) #RS only
                geno_order<-c("CCPP_RS","CPLY_RS","WWLL_RS","WLPY_RS")
		names(contrast_to_larson)<-c("CCPP_RSvsCCPPLY_RS", "CCPPLY_RSvsWWLL_RS", "WWLL_RSvsWWLLPY_RS", "CCPP_RSvsWWLLPY_RS")
        } else{
                stop("Invalid cell type")
        }
} else if(dataset=="Yintro_exp2"){
        if(cell_type=="LZandRS"){
                these_data<-mydata[,grepl("LLPP|PPLL|LPLY|PLPY",colnames(mydata))]
                these_lengths<-mylengths[,grepl("LLPP|PPLL|LPLY|PLPY",colnames(mylengths))]
                mygroups<-c(5,6,5,6,5,6,5,7,8,7,8,7,8,7,8,9,10,9,10,9,10,9,10,11,12,11,12,11,12,11,12)
                geno_order<-c("LLPP_LZ","LLPP_RS","LPLY_LZ","LPLY_RS","PPLL_LZ","PPLL_RS","PLPY_LZ","PLPY_RS")
        } else if(cell_type=="LZ"){
                these_data<-mydata[,grepl("LLPP.*LZ|PPLL.*LZ|LPLY.*LZ|PLPY.*LZ",colnames(mydata))]
                these_lengths<-mylengths[,grepl("LLPP.*LZ|PPLL.*LZ|LPLY.*LZ|PLPY.*LZ",colnames(mylengths))]
                mygroups<-c(5,5,5,5,7,7,7,7,9,9,9,9,11,11,11,11) #LZ only
                geno_order<-c("LLPP_LZ","LPLY_LZ","PPLL_LZ","PLPY_LZ")
		names(contrast_to_larson)<-c("PPLL_LZvsPPLLPY_LZ", "LLPPLY_LZvsPPLL_LZ", "LLPP_LZvsLLPPLY_LZ", "LLPP_LZvsPPLLPY_LZ")
        } else if(cell_type=="RS"){
                these_data<-mydata[,grepl("LLPP.*RS|PPLL.*RS|LPLY.*RS|PLPY.*RS",colnames(mydata))]
                these_lengths<-mylengths[,grepl("LLPP.*RS|PPLL.*RS|LPLY.*RS|PLPY.*RS",colnames(mylengths))]
                mygroups<-c(6,6,6,8,8,8,8,10,10,10,10,12,12,12,12) #RS only
                geno_order<-c("LLPP_RS","LPLY_RS","PPLL_RS","PLPY_RS")
		names(contrast_to_larson)<-c("PPLL_RSvsPPLLPY_RS", "LLPPLY_RSvsPPLL_RS", "LLPP_RSvsLLPPLY_RS", "LLPP_RSvsPPLLPY_RS")
        } else{
                stop("Invalid cell type")
        }
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
#print(dim(these_data))
#print(length(mygroups))
#print(head(these_data))
#print(is.numeric(these_data))
myDGE<-DGEList(counts=these_data, group=mygroups)
myDGE<-scaleOffset(myDGE, num_normMat)
print("Dimensions of myDGE counts:")
print(dim(myDGE$counts))
print(myDGE$counts[1:5,1:5])
print("myDGE samples - CHECK GROUP LABELS:")
print(myDGE$samples)
#This will be useful later to put meaningful labels on graphs, instead of the numeric group labels EdgeR uses
group_to_cross<-c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16)
names(group_to_cross)<-c("CCPP_LZ", "CCPP_RS", "CCPPLY_LZ", "CCPPLY_RS", "LLPP_LZ", "LLPP_RS", "LLPPLY_LZ", "LLPPLY_RS", "PPLLPY_LZ", "PPLLPY_RS", "PPLL_LZ", "PPLL_RS", "WWLL_LZ", "WWLL_RS", "WWLLPY_LZ", "WWLLPY_RS")

print("Filtering by expression level...")
#EdgeR rpkm function seems to require same gene length across all samples; taking median of tximport output gene lengths and using that
new_lengths<-rowMedians(as.matrix(these_lengths))
names(new_lengths)<-rownames(these_lengths)
print(head(new_lengths))
print(dim(myDGE$offset))
print(myDGE$samples$lib.size)
my_rpkm<-rpkm(myDGE, gene.length=new_lengths)
print(dim(my_rpkm))

keep<-rowSums(my_rpkm > min_rpkm) >= min_samples
keep[is.na(keep)]<-FALSE #NAs show up for gene IDs that aren't in myDGE; replace then with false otherwise myDGE gets confused

myDGE<-myDGE[keep, , keep.lib.sizes=FALSE]
myDGE<-calcNormFactors(myDGE)
print(myDGE$samples)
print("Dimensions of filtered myDGE counts:")
print(dim(myDGE$counts))

#Set up contrasts
print("Setting up contrasts...")
cross_type<-c() #Get cross names corresponding to groups
for(i in myDGE$samples$group){
	cross_type<-c(cross_type, names(group_to_cross)[which(group_to_cross == i)]) #Replace numbers w/ names
}
design_table<-as.data.frame(cbind(rownames(myDGE$samples), cross_type))
colnames(design_table)<-c("file","cross.type")
print(design_table)
design_matrix<-model.matrix(~0+design_table$cross.type, data=myDGE$samples)
design_table$cross.type<-as.factor(design_table$cross.type) #Needed in R version 4 or later
colnames(design_matrix)<-levels(design_table$cross.type)
print(design_matrix)
myDGE<-estimateDisp(myDGE,design_matrix)

myContrasts<-makeContrasts(c1=paste0(colnames(design_matrix)[1],"-",colnames(design_matrix)[2]), c2=paste0(colnames(design_matrix)[3],"-",colnames(design_matrix)[4]), c3=paste0(colnames(design_matrix)[1],"-",colnames(design_matrix)[3]), c4=paste0(colnames(design_matrix)[2],"-",colnames(design_matrix)[4]), c5=paste0(colnames(design_matrix)[1],"-",colnames(design_matrix)[4]), c6=paste0(colnames(design_matrix)[2],"-",colnames(design_matrix)[3]), levels=design_table$cross.type)
print(myContrasts)

#Identify most DE genes using qlf test
qlfit<-glmQLFit(myDGE,design_matrix)
qlf_numDE_table<-c()
#qlf_numDE_table_X<-c()
#qlf_numDE_table_auto<-c()
lfc_cutoff_qlf_numDE_table<-c()
qlf_tables<-NULL
for(i in 1:ncol(myContrasts)){
        result<-glmQLFTest(qlfit,contrast=myContrasts[,i])
        assign(paste("qlf",colnames(myContrasts)[i], sep="."),result)

        compare<-sub("-", "vs", colnames(myContrasts)[i])
        print(compare)

        tt<-topTags(result,n=nrow(result))
        assign(paste("qlf.tt",compare, sep="."), tt)
        assign(paste("qlf.DE",compare, sep="."), tt$table[tt$table$FDR<0.05,])
        assign(paste("qlf.numDE",compare, sep="."), length(which(tt$table$FDR<0.05)))
        qlf_numDE_table<-rbind(qlf_numDE_table,cbind(paste("qlf.numDE",compare, sep="."),length(which(tt$table$FDR<0.05))))
#        lfc_cutoff_qlf_numDE_table<-rbind(lfc_cutoff_qlf_numDE_table, cbind(paste("lfc.qlf.numDE",colnames(myContrasts)[i], sep="."), length(intersect(intersect(which(tt$table$FDR<0.05),which(rownames(tt$table) %in% pairwise_expressed[[compare]])), which(abs(tt$table$logFC)>min_logFC)))))
        qlf_tables[[colnames(myContrasts)[i]]]<-tt$table
}

#Identify most DE genes using the GLM and LRT
fit<-glmFit(myDGE,design_matrix)
lrt_numDE_table<-c()
lrt_numDE_table_X<-c()
lrt_numDE_table_auto<-c()
lfc_cutoff_lrt_numDE_table<-c()
lrt_tables<-NULL
for(i in 1:ncol(myContrasts)){
        result<-glmLRT(fit,contrast=myContrasts[,i])
        assign(paste("lrt",colnames(myContrasts)[i], sep="."),result)

        compare<-sub("-", "vs", colnames(myContrasts)[i])
        print(compare)

        tt<-topTags(result,n=nrow(result))
        assign(paste("lrt.tt",compare, sep="."), tt)
        assign(paste("lrt.DE",compare, sep="."), tt$table[tt$table$FDR<0.05,])
        assign(paste("lrt.numDE",compare, sep="."), length(which(tt$table$FDR<0.05)))
        lrt_numDE_table<-rbind(lrt_numDE_table,cbind(paste("lrt.numDE",compare, sep="."),length(which(tt$table$FDR<0.05))))
#        lfc_cutoff_lrt_numDE_table<-rbind(lfc_cutoff_lrt_numDE_table, cbind(paste("lfc.lrt.numDE",colnames(myContrasts)[i], sep="."), length(intersect(intersect(which(tt$table$FDR<0.05),which(rownames(tt$table) %in% pairwise_expressed[[compare]])), which(abs(tt$table$logFC)>min_logFC)))))
        lrt_tables[[colnames(myContrasts)[i]]]<-tt$table
}

print("DE genes:")
for(i in ls(pat="\\.DE\\.")){
	print(head(get(i)))
}

print("Genes with lowest FDR-corrected p-value:")
for(i in ls(pat="\\.tt\\.")){
	print(head(get(i)))
}

print("Number of DE genes:")
print(qlf_numDE_table)
print(lrt_numDE_table)

#Get chromosomes for all DE genes
#edb<-EnsDb.Mmusculus.v79
#all_genes<-genes(edb)
ens_mus<-useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl")
all_genes<-getBM(attributes=c('ensembl_gene_id','external_gene_name','chromosome_name','start_position','end_position','strand','gene_biotype'), mart=ens_mus)
colnames(all_genes)<-c("gene_id","gene_name","seqnames","start","end","strand","gene_biotype")
auto_genes<-all_genes[which(all_genes$seqname %in% c(1:19)),]
proco_genes<-all_genes$gene_id[which(all_genes$gene_biotype=="protein_coding")]

#Generate volcano plots
pdf(paste("volcano_plots",dataset,cell_type,"pdf", sep="."), onefile=TRUE)
for(c in names(contrast_to_larson)){
	i<-sub("vs","-",c)
	print(i)
	this_table<-lrt_tables[[i]]
	#print(dim(this_table))
	this_chrs<-c()
	chr_type<-c()
	larson_DE_data<-read.table(paste0("/mnt/beegfs/ek112884/amplicon_expression_analysis/SALMON_LarsonEtAl2017/topDEgenes.lrt.DE.",contrast_to_larson[c],".txt"), header=TRUE)
	#print(head(larson_DE_data))
	larson_DE<-c()
	for(j in rownames(this_table)){
		if(j %in% all_genes$gene_id){
			myChr<-all_genes$seqnames[which(all_genes$gene_id == j)]
			if(myChr %in% c(1:19)){
				chr_type<-c(chr_type, "auto")
			} else if(myChr == "X"){
				chr_type<-c(chr_type, "X")
			} else if(myChr == "Y"){
				chr_type<-c(chr_type, "Y")
			} else{
				chr_type<-c(chr_type, NA)
			}
		} else{
			myChr<-NA
			chr_type<-c(chr_type, NA)
		}
		this_chrs<-c(this_chrs, myChr)
		if(j %in% rownames(larson_DE_data)){
			if( (larson_DE_data$logFC[which(rownames(larson_DE_data)==j)] < 0) && (contrast_to_larson[c] %in% c("CCPP_RSvsPPLL_RS", "CCPP_RSvsLLPP_RS")) ){
				larson_DE<-c(larson_DE, "higher_in_hybrid")
			} else if( (larson_DE_data$logFC[which(rownames(larson_DE_data)==j)] > 0) && (contrast_to_larson[c] %in% c("CCPP_RSvsPPLL_RS", "CCPP_RSvsLLPP_RS")) ){
				larson_DE<-c(larson_DE, "lower_in_hybrid")
			} else if( (larson_DE_data$logFC[which(rownames(larson_DE_data)==j)] < 0) && (contrast_to_larson[c] %in% c("PPLL_RSvsWWLL_RS", "LLPP_RSvsWWLL_RS")) ){
				larson_DE<-c(larson_DE, "lower_in_hybrid")
			} else if( (larson_DE_data$logFC[which(rownames(larson_DE_data)==j)] > 0) && (contrast_to_larson[c] %in% c("PPLL_RSvsWWLL_RS", "LLPP_RSvsWWLL_RS")) ){
				larson_DE<-c(larson_DE, "higher_in_hybrid")
			} else{
				larson_DE<-c(larson_DE, "not_DE")
			}
		} else{
			larson_DE<-c(larson_DE, "not_DE")
		}
	}
	this_df<-as.data.frame(cbind(this_table, chr=this_chrs, chr_type, larson_DE))
	#print(head(this_df))
	#All chrs together
#	if(c %in% c("CCPP_RSvsCCPPLY_RS", "WWLL_RSvsWWLLPY_RS", "CCPP_RSvsWWLLPY_RS", "LLPPLY_RSvsPPLL_RS")){
#		#using -logFC on x-axis because edgeR comparison was done in the opposite direction from what we want
#		p<-ggplot(this_df, aes(x=-logFC, y=-log10(FDR), color=chr_type, shape=larson_DE)) + geom_point()
#		p<-p + labs(title=paste("Volcano plot:", i, "- REVERSE"), x="logFC", y="-log10(FDR-corrected P-value)")
#	} else{
#		p<-ggplot(this_df, aes(x=logFC, y=-log10(FDR), color=chr_type, shape=larson_DE)) + geom_point()
#		p<-p + labs(title=paste("Volcano plot:", i), x="logFC", y="-log10(FDR-corrected P-value)")
#	}
#	p<-p + theme(axis.text = element_text(size=18), axis.title = element_text(size=21), plot.title = element_text(size=36))
#	p<-p + theme_minimal() + scale_color_manual(values=c('grey','coral4', 'cornflowerblue'))
#	p<-p + geom_hline(yintercept = -log10(0.05))
#	print(p)
	#X and Y separate
	for(j in c("X","Y","auto")){
		if(c %in% c("CCPP_RSvsCCPPLY_RS", "WWLL_RSvsWWLLPY_RS", "CCPP_RSvsWWLLPY_RS", "LLPPLY_RSvsPPLL_RS")){
			#using -logFC on x-axis because edgeR comparison was done in the opposite direction from what we want
			p<-ggplot(this_df[which(this_df$chr_type==j),], aes(x=-logFC, y=-log10(FDR), shape=larson_DE)) + geom_point()
			p<-p + labs(title=paste("Volcano plot:", i, "- REVERSE; chromosome", j), x="logFC", y="-log10(FDR-corrected P-value)")
		} else{
			p<-ggplot(this_df[which(this_df$chr_type==j),], aes(x=logFC, y=-log10(FDR), shape=larson_DE)) + geom_point()
			p<-p + labs(title=paste("Volcano plot:", i, "; chromosome", j), x="logFC", y="-log10(FDR-corrected P-value)")
		}
		p<-p + theme(axis.text = element_text(size=18), axis.title = element_text(size=21), plot.title = element_text(size=36))
		p<-p + theme_minimal()
		p<-p + geom_hline(yintercept = -log10(0.05))
		print(p)
	}
}
dev.off()

#Generate barplots of DE genes and test for chromosomes significantly over- or under-enriched for DE genes
if(procoOnly){
        pdf(paste("prop_DE_genes.barplots.REVERSE",dataset,cell_type,"procoOnly.pdf", sep="."))
} else{
        pdf(paste("prop_DE_genes.barplots.REVERSE",dataset,cell_type,"pdf", sep="."))
}
all_chrs<-as.character(all_genes$seqnames)
prop_DE_up_all<-c()
prop_DE_down_all<-c()
for(i in ls(pat="\\.DE\\.")){
	print(i)
	if(procoOnly){
                myDE<-get(i)[which(rownames(get(i)) %in% proco_genes), ]
        } else{
                myDE<-get(i)
        }
	chrs<-c()
	direc<-c()
	print(paste("# higher in control:", length(which(myDE$logFC>0))))
	print(paste("# higher in Y introgression:", length(which(myDE$logFC<0))))
	pos_genes<-rownames(myDE)[which(myDE$logFC>0)]
	neg_genes<-rownames(myDE)[which(myDE$logFC<0)]
	for(j in rownames(myDE)){
		if(j %in% all_genes$gene_id){
			chrs<-c(chrs, all_chrs[which(all_genes$gene_id==j)])
		} else{
			chrs<-c(chrs, NA)
		}
		if(j %in% pos_genes){
			direc<-c(direc,"+")
		} else if(j %in% neg_genes){
			direc<-c(direc,"-")
		}
	}
	assign(paste0(i,".withChrs"), as.data.frame(cbind(myDE,chrs,direc)))
	print(head(get(paste0(i,".withChrs"))))
	if(procoOnly){
		write.table(get(paste0(i,".withChrs")), file=paste("topDEgenes",i,"procoOnly.txt",sep="."), quote=FALSE, sep="\t")
	} else{
		write.table(get(paste0(i,".withChrs")), file=paste("topDEgenes",i,"txt",sep="."), quote=FALSE, sep="\t")
	}
	print(paste("total DE:",length(chrs)))
	print(paste("# DE on X:",length(which(chrs=="X"))))
	print(paste("# DE on Y:",length(which(chrs=="Y"))))
	print(paste("# DE on 5:",length(which(chrs=="5"))))
	print(paste("# DE on 14:",length(which(chrs=="14"))))
	print(paste("# higher in control on X:",length(intersect(which(chrs=="X"),which(direc=="+")))))
	print(paste("# higher in Yinro on X:",length(intersect(which(chrs=="X"),which(direc=="-")))))
	print(paste("# higher in control on Y:",length(intersect(which(chrs=="Y"),which(direc=="+")))))
        print(paste("# higher in Yinro on Y:",length(intersect(which(chrs=="Y"),which(direc=="-")))))
	print(paste("# higher in control on 5:",length(intersect(which(chrs=="5"),which(direc=="+")))))
        print(paste("# higher in Yinro on 5:",length(intersect(which(chrs=="5"),which(direc=="-")))))
	print(paste("# higher in control on 14:",length(intersect(which(chrs=="14"),which(direc=="+")))))
        print(paste("# higher in Yinro on 14:",length(intersect(which(chrs=="14"),which(direc=="-")))))
	
	#Run hypergeometric tests and proportion tests
	#Hypergeometric for chromosomes enriched for DE genes
	#Proportion tests for difference in reciprocal cross directions for specific chromosomes
	#Number of DE genes across all autos
	auto_DE<-length(chrs) - length(which(chrs=="X")) - length(which(chrs=="Y"))
	#Number of expressed but not DE genes across all autos
	if(procoOnly){
		all_exp<-rownames(myDGE$counts)[which(rownames(myDGE$counts) %in% proco_genes)]
	} else{
		all_exp<-rownames(myDGE$counts) #all expressed genes
	}
	all_auto<-all_exp[which(all_exp %in% auto_genes$gene_id)] #autosomal only
	auto_notDE<-length(all_auto[which(!(all_auto %in% rownames(myDE)))]) #NOT DE
	print(auto_notDE)
	pvalues_lower<-c()
	pvalues_higher<-c()
	major_chrs<-c(1:19,"X","Y")
	prop_DE<-c()
	prop_DE_up<-c()
        prop_DE_down<-c()
        DE_up_autos<-0
        DE_down_autos<-0
	DE_count<-c()
        exp_count<-c()
	for(j in major_chrs){
		print(j)
        	#Number of DE genes on chromosome i
	        chr_DE<-length(which(chrs==j))
        	chr_DE_up<-length(intersect(which(chrs==j), which(direc=="+")))
                chr_DE_down<-length(intersect(which(chrs==j), which(direc=="-")))
		#Total number of expressed genes on chromosome i
	        chr_exp<-length(intersect(rownames(myDGE), all_genes$gene_id[which(all_genes$seqnames == j)]))
		print(paste(chr_DE, auto_DE, auto_notDE, chr_exp))
        	pvalues_lower<-c(pvalues_lower,phyper(chr_DE,auto_DE,auto_notDE,chr_exp))
	        pvalues_higher<-c(pvalues_higher,phyper(chr_DE,auto_DE,auto_notDE,chr_exp,lower.tail=FALSE))
		prop_DE<-c(prop_DE, chr_DE/chr_exp)
		prop_DE_up<-c(prop_DE_up, chr_DE_up/chr_exp)
                prop_DE_down<-c(prop_DE_down, chr_DE_down/chr_exp)
                if( (j != "X") && (j != "Y") ){
                        DE_up_autos<-DE_up_autos + chr_DE_up
                        DE_down_autos<-DE_down_autos + chr_DE_down
                }
		DE_count<-c(DE_count, chr_DE)
		exp_count<-c(exp_count, chr_exp)
	}
	names(prop_DE)<-major_chrs
	names(prop_DE_up)<-major_chrs
	names(prop_DE_down)<-major_chrs
	names(DE_count)<-major_chrs
	names(exp_count)<-major_chrs
        prop_DE_up_mergeAutos<-c(DE_up_autos/length(all_auto), prop_DE_up["X"], prop_DE_up["Y"])
        prop_DE_down_mergeAutos<-c(DE_down_autos/length(all_auto), prop_DE_down["X"], prop_DE_down["Y"])
        names(prop_DE_up_mergeAutos)<-c("auto","X","Y")
        names(prop_DE_down_mergeAutos)<-c("auto","X","Y")
        prop_DE_up_all<-c(prop_DE_up_all, prop_DE_up_mergeAutos)
	prop_DE_down_all<-c(prop_DE_down_all, prop_DE_down_mergeAutos)
	print("X vs autos:")
	print(prop.test(c(DE_count["X"], DE_up_autos+DE_down_autos), c(exp_count["X"], length(all_auto))))
	print("Y vs autos:")
	print(prop.test(c(DE_count["Y"], DE_up_autos+DE_down_autos), c(exp_count["Y"], length(all_auto))))
	if((DE_up_autos > 0) || (DE_down_autos > 0) ){
		print("Up vs down, autos:")
		print(prop.test(c(DE_up_autos, DE_down_autos), c(DE_up_autos+DE_down_autos, DE_up_autos+DE_down_autos)))
	}
	if((prop_DE_up["X"] > 0) || (prop_DE_down["X"] > 0) ){
		print("Up vs down, X:")
		print(prop.test(c(prop_DE_up["X"]*exp_count["X"], prop_DE_down["X"]*exp_count["X"]), c(DE_count["X"], DE_count["X"])))
	}
	if((prop_DE_up["Y"] > 0) || (prop_DE_down["Y"] > 0) ){
		print("Up vs down, Y:")
		print(prop.test(c(prop_DE_up["Y"]*exp_count["Y"], prop_DE_down["Y"]*exp_count["Y"]), c(DE_count["Y"], DE_count["Y"])))
	}
        barplot(-prop_DE_up_mergeAutos, main=paste(i, "- REVERSE"), ylim=c(-1,1))
        barplot(prop_DE_down_mergeAutos, add = TRUE)
	#Control for multiple test and print which cell type and chromosome significantly differs from expected
	padjust_all<-p.adjust(c(pvalues_lower, pvalues_higher), "fdr")
	padjust_lower<-padjust_all[1:length(pvalues_lower)]
	padjust_higher<-padjust_all[(length(pvalues_lower)+1):length(padjust_all)]
	#Print results
	print("Under-enriched:")
	for(i in 1:length(padjust_lower)){
		if(padjust_lower[i] < 0.05){
			print(paste0(major_chrs[i],": ",padjust_lower[i]))
		}
	}
	print("Over-enriched:")
	for(i in 1:length(padjust_higher)){
		if(padjust_higher[i] < 0.05){
			print(paste0(major_chrs[i],": ",padjust_higher[i]))
		}
	}
}
dev.off()

print("Done with 11_edgeR_DE.r")

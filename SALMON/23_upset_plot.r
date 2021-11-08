#PURPOSE: Generate upset plots for DE genes shared among cross types

library(UpSetR)

#Read in data
cell_type<-"RS"

print(paste("Reading in DE genes for cell type:", cell_type))

exp1_musX<-read.table(paste0("topDEgenes.lrt.DE.CCPP_",cell_type,"vsCCPPLY_",cell_type,".txt"), header=TRUE)
exp1_musY<-read.table(paste0("topDEgenes.lrt.DE.CCPP_",cell_type,"vsWWLLPY_",cell_type,".txt"), header=TRUE)
exp1_domX<-read.table(paste0("topDEgenes.lrt.DE.WWLL_",cell_type,"vsWWLLPY_",cell_type,".txt"), header=TRUE)
exp1_domY<-read.table(paste0("topDEgenes.lrt.DE.CCPPLY_",cell_type,"vsWWLL_",cell_type,".txt"), header=TRUE)

exp2_musX<-read.table(paste0("topDEgenes.lrt.DE.PPLL_",cell_type,"vsPPLLPY_",cell_type,".txt"), header=TRUE)
exp2_musY<-read.table(paste0("topDEgenes.lrt.DE.LLPP_",cell_type,"vsPPLLPY_",cell_type,".txt"), header=TRUE)
exp2_domX<-read.table(paste0("topDEgenes.lrt.DE.LLPP_",cell_type,"vsLLPPLY_",cell_type,".txt"), header=TRUE)
exp2_domY<-read.table(paste0("topDEgenes.lrt.DE.LLPPLY_",cell_type,"vsPPLL_",cell_type,".txt"), header=TRUE)

larson_musX<-read.table(paste0("../SALMON_LarsonEtAl2017/topDEgenes.lrt.DE.CCPP_",cell_type,"vsPPLL_",cell_type,".txt"), header=TRUE)
larson_musY<-read.table(paste0("../SALMON_LarsonEtAl2017/topDEgenes.lrt.DE.CCPP_",cell_type,"vsLLPP_",cell_type,".txt"), header=TRUE)
larson_domX<-read.table(paste0("../SALMON_LarsonEtAl2017/topDEgenes.lrt.DE.LLPP_",cell_type,"vsWWLL_",cell_type,".txt"), header=TRUE)
larson_domY<-read.table(paste0("../SALMON_LarsonEtAl2017/topDEgenes.lrt.DE.PPLL_",cell_type,"vsWWLL_",cell_type,".txt"), header=TRUE)

print("Separating by direction and chromosome type")
#Get that are consistently DE in the same direction; higher means higher exp in hybrid or XY mismatch
exp1_musX_higher<-rownames(exp1_musX)[which(exp1_musX$direc=="-")]
exp2_musX_higher<-rownames(exp2_musX)[which(exp2_musX$direc=="+")]
larson_musX_higher<-rownames(larson_musX)[which(larson_musX$direc=="-")]
exp1_musX_lower<-rownames(exp1_musX)[which(exp1_musX$direc=="+")]
exp2_musX_lower<-rownames(exp2_musX)[which(exp2_musX$direc=="-")]
larson_musX_lower<-rownames(larson_musX)[which(larson_musX$direc=="+")]
exp1_domX_higher<-rownames(exp1_domX)[which(exp1_domX$direc=="-")]
exp2_domX_higher<-rownames(exp2_domX)[which(exp2_domX$direc=="+")]
larson_domX_higher<-rownames(larson_domX)[which(larson_domX$direc=="+")]
exp1_domX_lower<-rownames(exp1_domX)[which(exp1_domX$direc=="+")]
exp2_domX_lower<-rownames(exp2_domX)[which(exp2_domX$direc=="-")]
larson_domX_lower<-rownames(larson_domX)[which(larson_domX$direc=="-")]
#Get list of genes that ure upregulated in one cross type but downregulated in other, depsite having the same sex chromosome combination
exclude_musX<-c(intersect(exp1_musX_higher, union(exp2_musX_lower, larson_musX_lower)), intersect(exp2_musX_higher, union(exp1_musX_lower, larson_musX_lower)), intersect(larson_musX_higher, union(exp1_musX_lower, exp2_musX_lower)), intersect(exp1_musX_lower, union(exp2_musX_higher, larson_musX_higher)), intersect(exp2_musX_lower, union(exp1_musX_higher, larson_musX_higher)), intersect(larson_musX_lower, union(exp1_musX_higher, exp2_musX_higher)))


exp1_musY_higher<-rownames(exp1_musY)[which(exp1_musY$direc=="-")]
exp2_musY_higher<-rownames(exp2_musY)[which(exp2_musY$direc=="+")]
larson_musY_higher<-rownames(larson_musY)[which(larson_musY$direc=="-")]
exp1_musY_lower<-rownames(exp1_musY)[which(exp1_musY$direc=="+")]
exp2_musY_lower<-rownames(exp2_musY)[which(exp2_musY$direc=="-")]
larson_musY_lower<-rownames(larson_musY)[which(larson_musY$direc=="+")]
exp1_domY_higher<-rownames(exp1_domY)[which(exp1_domY$direc=="+")]
exp2_domY_higher<-rownames(exp2_domY)[which(exp2_domY$direc=="-")]
larson_domY_higher<-rownames(larson_domY)[which(larson_domY$direc=="+")]
exp1_domY_lower<-rownames(exp1_domY)[which(exp1_domY$direc=="-")]
exp2_domY_lower<-rownames(exp2_domY)[which(exp2_domY$direc=="+")]
larson_domY_lower<-rownames(larson_domY)[which(larson_domY$direc=="-")]

#DE in same direction, regardless of direction
#This is going to be a bit complicated; need to remove only those genes that are DE in different directions but keep ones that are DE in one but not the other

#Separate by chromosome type
#Auto
exp1_musX_higher_auto<-intersect(exp1_musX_higher, rownames(exp1_musX)[which(exp1_musX$chrs %in% c(1:19))])
exp2_musX_higher_auto<-intersect(exp2_musX_higher, rownames(exp2_musX)[which(exp2_musX$chrs %in% c(1:19))])
larson_musX_higher_auto<-intersect(larson_musX_higher, rownames(larson_musX)[which(larson_musX$chrs %in% c(1:19))])
exp1_musX_lower_auto<-intersect(exp1_musX_lower, rownames(exp1_musX)[which(exp1_musX$chrs %in% c(1:19))])
exp2_musX_lower_auto<-intersect(exp2_musX_lower, rownames(exp2_musX)[which(exp2_musX$chrs %in% c(1:19))])
larson_musX_lower_auto<-intersect(larson_musX_lower, rownames(larson_musX)[which(larson_musX$chrs %in% c(1:19))])
exp1_domX_higher_auto<-intersect(exp1_domX_higher, rownames(exp1_domX)[which(exp1_domX$chrs %in% c(1:19))])
exp2_domX_higher_auto<-intersect(exp2_domX_higher, rownames(exp2_domX)[which(exp2_domX$chrs %in% c(1:19))])
larson_domX_higher_auto<-intersect(larson_domX_higher, rownames(larson_domX)[which(larson_domX$chrs %in% c(1:19))])
exp1_domX_lower_auto<-intersect(exp1_domX_lower, rownames(exp1_domX)[which(exp1_domX$chrs %in% c(1:19))])
exp2_domX_lower_auto<-intersect(exp2_domX_lower, rownames(exp2_domX)[which(exp2_domX$chrs %in% c(1:19))])
larson_domX_lower_auto<-intersect(larson_domX_lower, rownames(larson_domX)[which(larson_domX$chrs %in% c(1:19))])
exp1_musY_higher_auto<-intersect(exp1_musY_higher, rownames(exp1_musY)[which(exp1_musY$chrs %in% c(1:19))])
exp2_musY_higher_auto<-intersect(exp2_musY_higher, rownames(exp2_musY)[which(exp2_musY$chrs %in% c(1:19))])
larson_musY_higher_auto<-intersect(larson_musY_higher, rownames(larson_musY)[which(larson_musY$chrs %in% c(1:19))])
exp1_musY_lower_auto<-intersect(exp1_musY_lower, rownames(exp1_musY)[which(exp1_musY$chrs %in% c(1:19))])
exp2_musY_lower_auto<-intersect(exp2_musY_lower, rownames(exp2_musY)[which(exp2_musY$chrs %in% c(1:19))])
larson_musY_lower_auto<-intersect(larson_musY_lower, rownames(larson_musY)[which(larson_musY$chrs %in% c(1:19))])
exp1_domY_higher_auto<-intersect(exp1_domY_higher, rownames(exp1_domY)[which(exp1_domY$chrs %in% c(1:19))])
exp2_domY_higher_auto<-intersect(exp2_domY_higher, rownames(exp2_domY)[which(exp2_domY$chrs %in% c(1:19))])
larson_domY_higher_auto<-intersect(larson_domY_higher, rownames(larson_domY)[which(larson_domY$chrs %in% c(1:19))])
exp1_domY_lower_auto<-intersect(exp1_domY_lower, rownames(exp1_domY)[which(exp1_domY$chrs %in% c(1:19))])
exp2_domY_lower_auto<-intersect(exp2_domY_lower, rownames(exp2_domY)[which(exp2_domY$chrs %in% c(1:19))])
larson_domY_lower_auto<-intersect(larson_domY_lower, rownames(larson_domY)[which(larson_domY$chrs %in% c(1:19))])

#X
exp1_musX_higher_X<-intersect(exp1_musX_higher, rownames(exp1_musX)[which(exp1_musX$chrs=="X")])
exp2_musX_higher_X<-intersect(exp2_musX_higher, rownames(exp2_musX)[which(exp2_musX$chrs=="X")])
larson_musX_higher_X<-intersect(larson_musX_higher, rownames(larson_musX)[which(larson_musX$chrs=="X")])
exp1_musX_lower_X<-intersect(exp1_musX_lower, rownames(exp1_musX)[which(exp1_musX$chrs=="X")])
exp2_musX_lower_X<-intersect(exp2_musX_lower, rownames(exp2_musX)[which(exp2_musX$chrs=="X")])
larson_musX_lower_X<-intersect(larson_musX_lower, rownames(larson_musX)[which(larson_musX$chrs=="X")])
exp1_domX_higher_X<-intersect(exp1_domX_higher, rownames(exp1_domX)[which(exp1_domX$chrs=="X")])
exp2_domX_higher_X<-intersect(exp2_domX_higher, rownames(exp2_domX)[which(exp2_domX$chrs=="X")])
larson_domX_higher_X<-intersect(larson_domX_higher, rownames(larson_domX)[which(larson_domX$chrs=="X")])
exp1_domX_lower_X<-intersect(exp1_domX_lower, rownames(exp1_domX)[which(exp1_domX$chrs=="X")])
exp2_domX_lower_X<-intersect(exp2_domX_lower, rownames(exp2_domX)[which(exp2_domX$chrs=="X")])
larson_domX_lower_X<-intersect(larson_domX_lower, rownames(larson_domX)[which(larson_domX$chrs=="X")])

#Y
exp1_musY_higher_Y<-intersect(exp1_musY_higher, rownames(exp1_musY)[which(exp1_musY$chrs=="Y")])
exp2_musY_higher_Y<-intersect(exp2_musY_higher, rownames(exp2_musY)[which(exp2_musY$chrs=="Y")])
larson_musY_higher_Y<-intersect(larson_musY_higher, rownames(larson_musY)[which(larson_musY$chrs=="Y")])
exp1_musY_lower_Y<-intersect(exp1_musY_lower, rownames(exp1_musY)[which(exp1_musY$chrs=="Y")])
exp2_musY_lower_Y<-intersect(exp2_musY_lower, rownames(exp2_musY)[which(exp2_musY$chrs=="Y")])
larson_musY_lower_Y<-intersect(larson_musY_lower, rownames(larson_musY)[which(larson_musY$chrs=="Y")])
exp1_domY_higher_Y<-intersect(exp1_domY_higher, rownames(exp1_domY)[which(exp1_domY$chrs=="Y")])
exp2_domY_higher_Y<-intersect(exp2_domY_higher, rownames(exp2_domY)[which(exp2_domY$chrs=="Y")])
larson_domY_higher_Y<-intersect(larson_domY_higher, rownames(larson_domY)[which(larson_domY$chrs=="Y")])
exp1_domY_lower_Y<-intersect(exp1_domY_lower, rownames(exp1_domY)[which(exp1_domY$chrs=="Y")])
exp2_domY_lower_Y<-intersect(exp2_domY_lower, rownames(exp2_domY)[which(exp2_domY$chrs=="Y")])
larson_domY_lower_Y<-intersect(larson_domY_lower, rownames(larson_domY)[which(larson_domY$chrs=="Y")])

#Plot
print("Plotting...")
pdf(paste("upset_plot",cell_type,"pdf", sep="."), onefile=TRUE)
listInput_x_higher_auto<-list(mus_domYVSmus=exp1_musX_higher_auto, musXdomVSmusXdom_musY=exp2_musX_higher_auto, musXdomVSmus=larson_musX_higher_auto, dom_musYVSdom=exp1_domX_higher_auto, domXmusVSdomXmus_domY=exp2_domX_higher_auto, domXmusVSdom=larson_domX_higher_auto)
upset(fromList(listInput_x_higher_auto), sets=names(listInput_x_higher_auto), keep.order=TRUE, mainbar.y.label="Number of DE genes - autos - higher in mismatch or hybrid") #order.by="freq"
listInput_x_lower_auto<-list(mus_domYVSmus=exp1_musX_lower_auto, musXdomVSmusXdom_musY=exp2_musX_lower_auto, musXdomVSmus=larson_musX_lower_auto, dom_musYVSdom=exp1_domX_lower_auto, domXmusVSdomXmus_domY=exp2_domX_lower_auto, domXmusVSdom=larson_domX_lower_auto)
upset(fromList(listInput_x_lower_auto), sets=names(listInput_x_lower_auto), keep.order=TRUE, mainbar.y.label="Number of DE genes - autos - lower in mismatch or hybrid")
listInput_x_higher_X<-list(mus_domYVSmus=exp1_musX_higher_X, musXdomVSmusXdom_musY=exp2_musX_higher_X, musXdomVSmus=larson_musX_higher_X, dom_musYVSdom=exp1_domX_higher_X, domXmusVSdomXmus_domY=exp2_domX_higher_X, domXmusVSdom=larson_domX_higher_X)
upset(fromList(listInput_x_higher_X), sets=names(listInput_x_higher_X), keep.order=TRUE, mainbar.y.label="Number of DE genes - X chr - higher in mismatch or hybrid")
listInput_x_lower_X<-list(mus_domYVSmus=exp1_musX_lower_X, musXdomVSmusXdom_musY=exp2_musX_lower_X, musXdomVSmus=larson_musX_lower_X, dom_musYVSdom=exp1_domX_lower_X, domXmusVSdomXmus_domY=exp2_domX_lower_X, domXmusVSdom=larson_domX_lower_X)
#If there are no DE genes in this category, list will have a character(0), which upset can't process, so remove them here
listInput_x_lower_X_no0<-listInput_x_lower_X[lapply(listInput_x_lower_X, length)>0]
upset(fromList(listInput_x_lower_X_no0), sets=names(listInput_x_lower_X_no0), keep.order=TRUE, mainbar.y.label="Number of DE genes - X chr - lower in mismatch or hybrid")
listInput_y_higher_Y<-list(dom_musYVSmus=exp1_musY_higher_Y, domXmusVSmusXdom_musY=exp2_musY_higher_Y, domXmusVSmus=larson_musY_higher_Y, mus_domYVSdom=exp1_domY_higher_Y, musXdomVSdomXmus_domY=exp2_domY_higher_Y, musXdomVSdom=larson_domY_higher_Y)
upset(fromList(listInput_y_higher_Y), sets=names(listInput_y_higher_Y), keep.order=TRUE, mainbar.y.label="Number of DE genes - Y chr - higher in mismatch or hybrid")
listInput_y_lower_Y<-list(dom_musYVSmus=exp1_musY_lower_Y, domXmusVSmusXdom_musY=exp2_musY_lower_Y, domXmusVSmus=larson_musY_lower_Y, mus_domYVSdom=exp1_domY_lower_Y, musXdomVSdomXmus_domY=exp2_domY_lower_Y, musXdomVSdom=larson_domY_lower_Y)
upset(fromList(listInput_y_lower_Y), sets=names(listInput_y_lower_Y), keep.order=TRUE, mainbar.y.label="Number of DE genes - Y chr - lower in mismatch or hybrid")
dev.off()

#Repeat WITHOUT directionality (all DE genes combined in one plot, regardless of if they're over- or underexpressed)
#Get list of genes that ure upregulated in one cross type but downregulated in other, depsite having the same sex chromosome combination
exclude_musX<-c(intersect(exp1_musX_higher, union(exp2_musX_lower, larson_musX_lower)), intersect(exp2_musX_higher, union(exp1_musX_lower, larson_musX_lower)), intersect(larson_musX_higher, union(exp1_musX_lower, exp2_musX_lower)), intersect(exp1_musX_lower, union(exp2_musX_higher, larson_musX_higher)), intersect(exp2_musX_lower, union(exp1_musX_higher, larson_musX_higher)), intersect(larson_musX_lower, union(exp1_musX_higher, exp2_musX_higher)))
exclude_domX<-c(intersect(exp1_domX_higher, union(exp2_domX_lower, larson_domX_lower)), intersect(exp2_domX_higher, union(exp1_domX_lower, larson_domX_lower)), intersect(larson_domX_higher, union(exp1_domX_lower, exp2_domX_lower)), intersect(exp1_domX_lower, union(exp2_domX_higher, larson_domX_higher)), intersect(exp2_domX_lower, union(exp1_domX_higher, larson_domX_higher)), intersect(larson_domX_lower, union(exp1_domX_higher, exp2_domX_higher)))
exclude_musY<-c(intersect(exp1_musY_higher, union(exp2_musY_lower, larson_musY_lower)), intersect(exp2_musY_higher, union(exp1_musY_lower, larson_musY_lower)), intersect(larson_musY_higher, union(exp1_musY_lower, exp2_musY_lower)), intersect(exp1_musY_lower, union(exp2_musY_higher, larson_musY_higher)), intersect(exp2_musY_lower, union(exp1_musY_higher, larson_musY_higher)), intersect(larson_musY_lower, union(exp1_musY_higher, exp2_musY_higher)))
exclude_domY<-c(intersect(exp1_domY_higher, union(exp2_domY_lower, larson_domY_lower)), intersect(exp2_domY_higher, union(exp1_domY_lower, larson_domY_lower)), intersect(larson_domY_higher, union(exp1_domY_lower, exp2_domY_lower)), intersect(exp1_domY_lower, union(exp2_domY_higher, larson_domY_higher)), intersect(exp2_domY_lower, union(exp1_domY_higher, larson_domY_higher)), intersect(larson_domY_lower, union(exp1_domY_higher, exp2_domY_higher)))
print(length(exclude_musX))

#Get DE genes in both directions for each sex chromosome comparison, removing gene DE in opposite directions
exp1_musX_auto0<-c(exp1_musX_higher_auto, exp1_musX_lower_auto)
exp1_musX_auto<-exp1_musX_auto0[which(!(exp1_musX_auto0 %in% exclude_musX))]
exp1_musX_X0<-c(exp1_musX_higher_X, exp1_musX_lower_X)
exp1_musX_X<-exp1_musX_X0[which(!(exp1_musX_X0 %in% exclude_musX))]
exp1_domX_auto0<-c(exp1_domX_higher_auto, exp1_domX_lower_auto)
exp1_domX_auto<-exp1_domX_auto0[which(!(exp1_domX_auto0 %in% exclude_domX))]
exp1_domX_X0<-c(exp1_domX_higher_X, exp1_domX_lower_X)
exp1_domX_X<-exp1_domX_X0[which(!(exp1_domX_X0 %in% exclude_domX))]
exp1_musY_auto0<-c(exp1_musY_higher_auto, exp1_musY_lower_auto)
exp1_musY_auto<-exp1_musY_auto0[which(!(exp1_musY_auto0 %in% exclude_musY))]
exp1_musY_Y0<-c(exp1_musY_higher_Y, exp1_musY_lower_Y)
exp1_musY_Y<-exp1_musY_Y0[which(!(exp1_musY_Y0 %in% exclude_musY))]
exp1_domY_auto0<-c(exp1_domY_higher_auto, exp1_domY_lower_auto)
exp1_domY_auto<-exp1_domY_auto0[which(!(exp1_domY_auto0 %in% exclude_domY))]
exp1_domY_Y0<-c(exp1_domY_higher_Y, exp1_domY_lower_Y)
exp1_domY_Y<-exp1_domY_Y0[which(!(exp1_domY_Y0 %in% exclude_domY))]

exp2_musX_auto0<-c(exp2_musX_higher_auto, exp2_musX_lower_auto)
exp2_musX_auto<-exp2_musX_auto0[which(!(exp2_musX_auto0 %in% exclude_musX))]
exp2_musX_X0<-c(exp2_musX_higher_X, exp2_musX_lower_X)
exp2_musX_X<-exp2_musX_X0[which(!(exp2_musX_X0 %in% exclude_musX))]
exp2_domX_auto0<-c(exp2_domX_higher_auto, exp2_domX_lower_auto)
exp2_domX_auto<-exp2_domX_auto0[which(!(exp2_domX_auto0 %in% exclude_domX))]
exp2_domX_X0<-c(exp2_domX_higher_X, exp2_domX_lower_X)
exp2_domX_X<-exp2_domX_X0[which(!(exp2_domX_X0 %in% exclude_domX))]
exp2_musY_auto0<-c(exp2_musY_higher_auto, exp2_musY_lower_auto)
exp2_musY_auto<-exp2_musY_auto0[which(!(exp2_musY_auto0 %in% exclude_musY))]
exp2_musY_Y0<-c(exp2_musY_higher_Y, exp2_musY_lower_Y)
exp2_musY_Y<-exp2_musY_Y0[which(!(exp2_musY_Y0 %in% exclude_musY))]
exp2_domY_auto0<-c(exp2_domY_higher_auto, exp2_domY_lower_auto)
exp2_domY_auto<-exp2_domY_auto0[which(!(exp2_domY_auto0 %in% exclude_domY))]
exp2_domY_Y0<-c(exp2_domY_higher_Y, exp2_domY_lower_Y)
exp2_domY_Y<-exp2_domY_Y0[which(!(exp2_domY_Y0 %in% exclude_domY))]

larson_musX_auto0<-c(larson_musX_higher_auto, larson_musX_lower_auto)
larson_musX_auto<-larson_musX_auto0[which(!(larson_musX_auto0 %in% exclude_musX))]
larson_musX_X0<-c(larson_musX_higher_X, larson_musX_lower_X)
larson_musX_X<-larson_musX_X0[which(!(larson_musX_X0 %in% exclude_musX))]
larson_domX_auto0<-c(larson_domX_higher_auto, larson_domX_lower_auto)
larson_domX_auto<-larson_domX_auto0[which(!(larson_domX_auto0 %in% exclude_domX))]
larson_domX_X0<-c(larson_domX_higher_X, larson_domX_lower_X)
larson_domX_X<-larson_domX_X0[which(!(larson_domX_X0 %in% exclude_domX))]
larson_musY_auto0<-c(larson_musY_higher_auto, larson_musY_lower_auto)
larson_musY_auto<-larson_musY_auto0[which(!(larson_musY_auto0 %in% exclude_musY))]
larson_musY_Y0<-c(larson_musY_higher_Y, larson_musY_lower_Y)
larson_musY_Y<-larson_musY_Y0[which(!(larson_musY_Y0 %in% exclude_musY))]
larson_domY_auto0<-c(larson_domY_higher_auto, larson_domY_lower_auto)
larson_domY_auto<-larson_domY_auto0[which(!(larson_domY_auto0 %in% exclude_domY))]
larson_domY_Y0<-c(larson_domY_higher_Y, larson_domY_lower_Y)
larson_domY_Y<-larson_domY_Y0[which(!(larson_domY_Y0 %in% exclude_domY))]

print("Plotting (any direction)...")
pdf(paste("upset_plot.noDirectionality",cell_type,"pdf", sep="."), onefile=TRUE)
listInput_x_auto<-list(mus_domYVSmus=exp1_musX_auto, musXdomVSmusXdom_musY=exp2_musX_auto, musXdomVSmus=larson_musX_auto, dom_musYVSdom=exp1_domX_auto, domXmusVSdomXmus_domY=exp2_domX_auto, domXmusVSdom=larson_domX_auto)
upset(fromList(listInput_x_auto), sets=names(listInput_x_auto), keep.order=TRUE, mainbar.y.label="Number of DE genes - autos, same X") #order.by="freq"
listInput_y_auto<-list(dom_musYVSmus=exp1_musY_auto, domXmusVSmusXdom_musY=exp2_musY_auto, domXmusVSmus=larson_musY_auto, mus_domYVSdom=exp1_domY_auto, musXdomVSdomXmus_domY=exp2_domY_auto, musXdomVSdom=larson_domY_auto)
upset(fromList(listInput_y_auto), sets=names(listInput_y_auto), keep.order=TRUE, mainbar.y.label="Number of DE genes - autos,  same Y")
listInput_x_X<-list(mus_domYVSmus=exp1_musX_X, musXdomVSmusXdom_musY=exp2_musX_X, musXdomVSmus=larson_musX_X, dom_musYVSdom=exp1_domX_X, domXmusVSdomXmus_domY=exp2_domX_X, domXmusVSdom=larson_domX_X)
upset(fromList(listInput_x_X), sets=names(listInput_x_X), keep.order=TRUE, mainbar.y.label="Number of DE genes - X chr")
listInput_y_Y<-list(dom_musYVSmus=exp1_musY_Y, domXmusVSmusXdom_musY=exp2_musY_Y, domXmusVSmus=larson_musY_Y, mus_domYVSdom=exp1_domY_Y, musXdomVSdomXmus_domY=exp2_domY_Y, musXdomVSdom=larson_domY_Y)
upset(fromList(listInput_y_Y), sets=names(listInput_y_Y), keep.order=TRUE, mainbar.y.label="Number of DE genes - Y chr")
dev.off()

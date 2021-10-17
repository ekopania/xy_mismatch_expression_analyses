#PURPOSE: Generate volcano plots for expression datasets

#NOTES:
	#Figure out how to show X vs autos
	#Figure out how to show genes DE in Erica's dataset (especially musXdom vs mus on X chromosome)
	#Show both X versus autos and DE in Erica's dataset in same plot? Can do with a combination of colors and shape probably, if the program will let you do that, or can figure out my own code for doing that; volcano plots typically color by significance so this could be confusing but we might be able to get around this by just drawing a horizontal abline that is the significance threshold at whatever -log10(0.05) is
	#Volcano plot is -log10(p-value) by logFC; might be good to also do an MA or MD plot (logFC by expression level)

library(ggplot2)

comparisons<-c("mus+domY vs mus", "mus+domY vs dom", "dom+musY vs dom", "dom+musY vs mus")
comp_to_file<-c("CCPP_RSvsCCPPLY_RS", "CCPPLY_RSvsWWLL_RS", "WWLL_RSvsWWLLPY_RS", "CCPP_RSvsWWLLPY_RS")
names(comp_to_file)<-comparisons
pdf("volcano_plots.Yintro.pdf", onefile=TRUE)
for(c in comparisons){
	#Read in data
	mydata<-read.table(paste("topDEgenes.lrt.DE",comp_to_file[c],"txt", sep="."), header=TRUE)
	isX<-mydata$chrs=="X"
	final_df<-as.data.frame(cbind(mydata, isX))
	#plot
	if(c %in% c("mus+domY vs mus", "dom+musY vs dom", "dom+musY vs mus")){
		#using -logFC on x-axis because edgeR comparison was done in the opposite direction from what we want
		#using FDR-corrected P-value on y-axis
		p<-ggplot(mydata, aes(x=-logFC, y=-log10(FDR), color=isX)) + geom_point()
		p<-p + labs(title=paste("Volcano plot:", c), x="logFC", y="-log10(FDR-corrected P-value)")
		p<-p + theme(axis.text = element_text(size=18), axis.title = element_text(size=21), plot.title = element_text(size=36))
		p<-p + theme_minimal()
		print(p)
	} else{
		#using FDR-corrected P-value on y-axis
		p<-ggplot(mydata, aes(x=logFC, y=-log10(FDR), color=isX)) + geom_point()
		p<-p + labs(title=paste("Volcano plot:", c), x="logFC", y="-log10(FDR-corrected P-value)")
		p<-p + theme(axis.text = element_text(size=18), axis.title = element_text(size=21), plot.title = element_text(size=36))
		p<-p + theme_minimal()
		print(p)
	}
}
dev.off()

print("Done with 17_volcano_plots.r") 

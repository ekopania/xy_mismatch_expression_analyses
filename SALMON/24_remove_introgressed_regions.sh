#!/bin/bash
#PURPOSE: Remove putative introgressed regions from lists of DE genes
# Job name:
#SBATCH --job-name=faster_remove_introgressed_regions
#SBATCH --output=remove_introgressed_regions-%j.log
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
##SBATCH --mail-user=ekopania4@gmail.com # Where to send mail
#SBATCH --cpus-per-task=1 # Number of cores per MPI rank (ie number of threads, I think)
#SBATCH --nodes=1 #Number of nodes
#SBATCH --ntasks=1 # Number of MPI ranks (ie number of processes, I think)
#SBATCH --mem-per-cpu=8G #Not sure if I should mess with these...
#SBATCH --mem=0 #Not sure if I should mess with these...
# Partition:
## Since you want to run it on 72 cores, the partition good_cpu has nodes with 72 cores.
#SBATCH --partition=good_lab_reincarnation
##SBATCH -w, --nodelist=compute-0-4 # run on a specific node
#
## Command(s) to run:
yintro="LEWPY" #LEWPY or PWKLY

#Comparisons involving LEWPY only
ls topDEgenes.lrt.DE.CCPP_??vsWWLLPY_??.txt topDEgenes.lrt.DE.WWLL_??vsWWLLPY_??.txt topDEgenes.lrt.DE.PPLL_??vsPPLLPY_??.txt topDEgenes.lrt.DE.LLPP_??vsPPLLPY_??.txt | while read file; do
#ls topDEgenes.lrt.DE.CCPP_RSvsWWLLPY_RS.txt | while read file; do
#Comparisons involving PWKLY only
#ls topDEgenes.lrt.DE.CCPP_??vsCCPPLY_??.txt topDEgenes.lrt.DE.CCPPLY_??vsWWLL_??.txt topDEgenes.lrt.DE.LLPP_??vsLLPPLY_??.txt topDEgenes.lrt.DE.LLPPLY_??vsPPLL_??.txt | while read file; do 
#ls topDEgenes.lrt.DE.CCPP_RSvsCCPPLY_RS.txt | while read file; do
	name=$(echo "${file}" | cut -d "." -f 1-4)
	echo "${name}"
	rm temp2.bed
	#Get DEgenes file into bed format
	awk '{if (NR!=1) {print $1}}' "${file}" > temp_genes.txt
	grep -f temp_genes.txt /mnt/beegfs/ek112884/REFERENCE_DIR/Mus_musculus.GRCm38.102.gtf | grep -P "\tgene\t" | awk '{print $1 "\t" $4 "\t" $5 "\t" $10}' > temp.bed
	cat temp_genes.txt | while read gene; do
		if grep "${gene}" temp.bed >> temp2.bed; then
			echo "${gene} found in gtf"
		else
			echo -e "0\t0\t0\t0" >> temp2.bed
		fi
	done
	#awk '{print $1}' "${file}" | while read gene; do
	#	if grep "${gene}" /mnt/beegfs/ek112884/REFERENCE_DIR/Mus_musculus.GRCm38.102.gtf | grep -P "\tgene\t" | awk '{print $1 "\t" $4 "\t" $5}' >> temp.bed; then
	#		echo "${gene} found in gtf"
	#	else
	#		echo -e "0\t0\t0" >> temp.bed
	#	fi
		#grep "${gene}" /mnt/beegfs/ek112884/REFERENCE_DIR/Mus_musculus.GRCm38.102.gtf | grep -P "\tgene\t" | awk '{print $1 "\t" $4 "\t" $5}' >> temp.bed
	#done
	echo -e "chr\tstart\tstop\tgeneID" | cat - temp2.bed > temp3.bed
	paste temp3.bed "${file}" > temp_file.bed
	#Remove DE genes that fall into putatively introgressed regions
	awk '{print $2 "\t" $3 "\t" $4 "\t" $5}' putative_introgressed_regions.${yintro}.100kbWindows.removeBadRegions.bed > temp_introgressed.bed
	bedtools subtract -a temp_file.bed -b temp_introgressed.bed -A > "temp.removeIntrogressed.bed"
        awk '{print $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" $12 "\t" $13}' temp.removeIntrogressed.bed > temp.removeIntrogressed2.bed
        head -1 "${file}" | cat - temp.removeIntrogressed2.bed > "${name}.removeIntrogressed.bed"
done

echo "Done!"

#Flags explaination:
	#Bedtools subtract removes features from a that are found in b
	#-a and -b flags establish which file is "a" and which is "b"
	#By default, only removes part of a that is in b (ex: if a has 100 200 and b has 180 300, output will be 100 180)
	#-A flag removes entire feature from a if there is overlap with b (in above example, entire 100 200 would be removed)
	#Set -A flag in order to remove entire DE gene if it overlaps a putatively introgressed region

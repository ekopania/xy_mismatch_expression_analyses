#!/bin/bash
#PURPOSE: Remove putative introgressed regions from lists of DE genes
# Job name:
#SBATCH --job-name=remove_introgressed_regions
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
yintro="LEWPY"

#Comparisons involving LEWPY only
ls topDEgenes.lrt.DE.CCPP_RSvsWWLLPY_RS.txt topDEgenes.lrt.DE.WWLL_RSvsWWLLPY_RS.txt topDEgenes.lrt.DE.PPLL_RSvsPPLLPY_RS.txt topDEgenes.lrt.DE.LLPP_RSvsPPLLPY_RS.txt | while read file; do
#Comparisons involving PWKLY only
#ls topDEgenes.lrtDE.CCPP_RSvsCCPPLY_RS.txt topDEgenes.lrt.DE.CCPPLY_RSvsWWLL_RS.txt topDEgenes.lrt.DE.LLPP_RSvsLLPPLY_RS.txt topDEgenes.lrt.DE.LLPPLY_RSvsPPLL_RS.txt | while read file; do 
	name=$(echo "${file}" | cut -d "." -f 1-4)
	echo "${name}"
	rm temp.bed
	#Get DEgenes file into bed format
	awk '{print $1}' "${file}" | while read gene; do
		if grep "${gene}" /mnt/beegfs/ek112884/REFERENCE_DIR/Mus_musculus.GRCm38.102.gtf | grep -P "\tgene\t" | awk '{print $1 "\t" $4 "\t" $5}' >> temp.bed; then
			echo "${gene} found in gtf"
		else
			echo -e "0\t0\t0" >> temp.bed
		fi
		#grep "${gene}" /mnt/beegfs/ek112884/REFERENCE_DIR/Mus_musculus.GRCm38.102.gtf | grep -P "\tgene\t" | awk '{print $1 "\t" $4 "\t" $5}' >> temp.bed
	done
	echo -e "chr\tstart\tstop" | cat - temp.bed > temp2.bed
	paste temp2.bed "${file}" > temp_file.bed
	#Remove DE genes that fall into putatively introgressed regions
	bedtools subtract -a temp_file.bed -b putative_introgressed_regions.${yintro}.bed -A > "${name}.removeIntrogressed.${yintro}.bed"
done

echo "Done!"

#Flags explaination:
	#Bedtools subtract removes features from a that are found in b
	#-a and -b flags establish which file is "a" and which is "b"
	#By default, only removes part of a that is in b (ex: if a has 100 200 and b has 180 300, output will be 100 180)
	#-A flag removes entire feature from a if there is overlap with b (in above example, entire 100 200 would be removed)
	#Set -A flag in order to remove entire DE gene if it overlaps a putatively introgressed region

#!/bin/bash
#PURPOSE: Partition genome into windows and count variants in each window
#
# Job name:
#SBATCH --job-name=count_vars
#SBATCH --output=count_vars-%j.log
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
##SBATCH --mail-user=ekopania4@gmail.com # Where to send mail
#SBATCH --cpus-per-task=1 # Number of cores per MPI rank (ie number of threads, I think)
#SBATCH --nodes=1 #Number of nodes
#SBATCH --ntasks=1 # Number of MPI ranks (ie number of processes, I think)
#SBATCH --mem-per-cpu=40G #Not sure if I should mess with these...
#SBATCH --mem=96000 #Not sure if I should mess with these...
# Partition:
## Since you want to run it on 72 cores, the partition good_cpu has nodes with 72 cores.
#SBATCH --partition=good_lab_cpu
##SBATCH -w, --nodelist=compute-0-4 # run on a specific node
#
## Command(s) to run:
win_size=10 #Window size in kb
#chr=1
#genome="REFS/LEWES_EiJ.pseudogenome.forBedtools.genome"
#genome="REFS/WSB.newVariants.forBedtools.genome"
genome="REFS/PWK_PhJ.hapCallerRef.forBedtools.genome"
ref="PWK"
sample="PWKLY"
in_vcf="${sample}.to${ref}.filter.vcf"

#Make windows
if [ -f "windows.${sample}.to${ref}.${win_size}kb.bed" ]; then
	echo "Already made windows"
else
	echo "Making windows"
	bp=$(( $win_size*1000 ))
	bedtools makewindows -g ${genome} -w ${bp} > "windows.${sample}.to${ref}.${win_size}kb.bed"  #-s 500 add to command for sliding windows
fi

#Count SNPs and indels in each window
if [ -f "counts.${sample}.to${ref}.${win_size}kbWindows.txt" ]; then
	echo "Already counted"
else
	echo "Counting VCF overlap"
	bedtools coverage -a "windows.${sample}.to${ref}.${win_size}kb.bed" -b ${in_vcf} -counts > "counts.${sample}.to${ref}.${win_size}kbWindows.txt"
fi

#Plot number of variants in each window
Rscript 08_plot_vars.r "counts.${sample}.to${ref}.${win_size}kbWindows.txt"

echo "Done!"

#!/bin/bash
#PURPOSE: Merge sorted bams, following bwa mapping for separate sequencing runs
#
# Job name:
#SBATCH --job-name=merge_bams
#SBATCH --output=merge_bams-%j.log
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=emily.kopania@umconnect.umt.edu # Where to send mail
#SBATCH --cpus-per-task=8 # Number of cores per MPI rank (ie number of threads, I think)
#SBATCH --nodes=1 #Number of nodes
#SBATCH --ntasks=1 # Number of MPI ranks (ie number of processes, I think)
##SBATCH --mem-per-cpu=8G #Not sure if I should mess with these...
##SBATCH --mem=8G #Not sure if I should mess with these...
# Partition:
## Since you want to run it on 72 cores, the partition good_cpu has nodes with 72 cores.
#SBATCH --partition=good_lab_reincarnation
##SBATCH -w, --nodelist=compute-0-4 # run on a specific node
#
## Command(s) to run:
out_prefix="LEWES"

#Merge separate runs that correspond to same sample
samples=( $(awk '{print $2}' run_to_sample.txt | tail -n +2 | sort | uniq) )
for i in "${samples[@]}"; do
	echo "${i}"
	#Don't waste time on ones that already exist
	if [ -f "MAPPED/${i}.to${out_prefix}.bam" ]; then
		echo "bam exists"
	else
	        rm temp_bam_list.to${out_prefix}.txt
        	grep "${i}" run_to_sample.txt | while read line; do
                	accession=$(echo "${line}" | cut -f 1)
	                echo "${accession}"
                	ls MAPPED/${accession}*.to${out_prefix}.bam >> temp_bam_list.to${out_prefix}.txt
        	done
		#-r flag includes a readgroup from separate files (sequencing runs) - important for proper downstream processing like marking dups!!!
	        samtools merge -r -b temp_bam_list.to${out_prefix}.txt "MAPPED/${i}.to${out_prefix}.bam"
        	samtools index "MAPPED/${i}.to${out_prefix}.bam"
        fi
done
rm temp_bam_list.to${out_prefix}.txt

echo "Done!"

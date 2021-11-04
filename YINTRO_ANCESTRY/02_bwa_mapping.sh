#!/bin/bash
#PURPOSE: Map reads using bwa
#NOTE: Make sure you've created bwa index for reference
#
# Job name:
#SBATCH --job-name=bwa_mapping
#SBATCH --output=bwa_mapping_output-%j.log
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=emily.kopania@umconnect.umt.edu # Where to send mail
#SBATCH --cpus-per-task=8 # Number of cores per MPI rank (ie number of threads, I think)
#SBATCH --nodes=1 #Number of nodes
#SBATCH --ntasks=1 # Number of MPI ranks (ie number of processes, I think)
##SBATCH --mem-per-cpu=8G #Not sure if I should mess with these...
##SBATCH --mem=8G #Not sure if I should mess with these...
# Partition:
## Since you want to run it on 72 cores, the partition good_cpu has nodes with 72 cores.
#SBATCH --partition=good_lab_cpu
##SBATCH -w, --nodelist=compute-0-4 # run on a specific node
#
## Command(s) to run:
#ref_prefix="LEWES_EiJ.pseudogenome.fa"
#ref_prefix="WSB.newVariants.fa"
ref_prefix="PWK_PhJ.hapCallerRef.fa"
out_prefix="PWK"

#Align trimmed paired reads, sort, index
#ls /mnt/beegfs/ek112884/cnvs/DATA/FASTQS/WDIS_WGS/RESEQ/LLLL.trimmed_1P.fq.gz | while read file; do
#ls /mnt/beegfs/ek112884/cnvs/DATA/FASTQS/WDIS_WGS/LLLL*trimmed_1P.fq.gz | while read file; do
#ls /mnt/beegfs/ek112884/cnvs/DATA/FASTQS/WDIS_WGS/RESEQ/LEWPY.trimmed_1P.fq.gz | while read file; do
#ls /mnt/beegfs/ek112884/cnvs/DATA/FASTQS/WDIS_WGS/LEWPY*trimmed_1P.fq.gz | while read file; do
#ls /mnt/beegfs/ek112884/cnvs/DATA/FASTQS/WDIS_WGS/RESEQ/PWKLY.trimmed_1P.fq.gz | while read file; do
ls /mnt/beegfs/ek112884/cnvs/DATA/FASTQS/WDIS_WGS/PWKLY*trimmed_1P.fq.gz | while read file; do
	name=$(echo "${file}" | cut -d "." -f 1)
	echo "${name}"
	f2="${name}.trimmed_2P.fq.gz"
	echo "${f2}"
	short_name=$(echo "${name}" | cut -d "/" -f 9) #10
	#add "_RESEQ" to short_name in filenames if running from the resequenced dataset
	bwa mem -t 8 REFS/${ref_prefix} "${file}" "${f2}" | samtools sort -@8 -o "MAPPED/${short_name}.to${out_prefix}.bam" -
	samtools index "MAPPED/${short_name}.to${out_prefix}.bam"
done

#Align forward unpaired reads (after trimming) as if they are SE; sort and index
#ls /mnt/beegfs/ek112884/cnvs/DATA/FASTQS/WDIS_WGS/RESEQ/LLLL.trimmed_1U.fq.gz | while read file; do
#ls /mnt/beegfs/ek112884/cnvs/DATA/FASTQS/WDIS_WGS/LLLL*trimmed_1U.fq.gz | while read file; do
#ls /mnt/beegfs/ek112884/cnvs/DATA/FASTQS/WDIS_WGS/RESEQ/LEWPY.trimmed_1U.fq.gz | while read file; do
#ls /mnt/beegfs/ek112884/cnvs/DATA/FASTQS/WDIS_WGS/LEWPY*trimmed_1U.fq.gz | while read file; do
#ls /mnt/beegfs/ek112884/cnvs/DATA/FASTQS/WDIS_WGS/RESEQ/PWKLY.trimmed_1U.fq.gz | while read file; do
ls /mnt/beegfs/ek112884/cnvs/DATA/FASTQS/WDIS_WGS/PWKLY*trimmed_1U.fq.gz | while read file; do
	short_name=$(echo "${file}" | cut -d "/" -f 9 | cut -d "." -f 1) #10
        echo "${short_name}"
	#add "_RESEQ" to short_name in filenames if running from the resequenced dataset
	bwa mem -t 8 REFS/${ref_prefix} "${file}" | samtools sort -@8 -o "MAPPED/${short_name}.1U.to${out_prefix}.bam" -
	samtools index "MAPPED/${short_name}.1U.to${out_prefix}.bam"
done

echo "Done!"

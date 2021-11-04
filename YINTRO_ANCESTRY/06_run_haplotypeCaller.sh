#!/bin/bash
#PURPOSE: Call variants in Y-introgression samples using HaplotypeCaller and GenotypeGVCFs from GATK4
#
# Job name:
#SBATCH --job-name=call_variants
#SBATCH --output=call_variants-%j.log
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=emily.kopania@umconnect.umt.edu # Where to send mail
#SBATCH --cpus-per-task=64 # Number of cores per MPI rank (ie number of threads, I think)
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
source ~/software/anaconda/anaconda3/bin/activate
conda activate sambcfenv #bcftools and samtools 1.13; gatk4

#ref_fa="REFS/WSB.newVariants.fa"
#ref_fa="REFS/LEWES_EiJ.pseudogenome.fa"
ref_fa="REFS/PWK_PhJ.hapCallerRef.fa"
echo "Input reference fasta: ${ref_fa}"

#ls MAPPED/LLLL.toWSB.marked_duplicates.coordsort.bam | while read file; do
#ls MAPPED/LLLL.toLEWES.marked_duplicates.coordsort.bam | while read file; do
#ls MAPPED/LEWPY.toLEWES.marked_duplicates.coordsort.bam | while read file; do
ls MAPPED/LEWPY.toPWK.marked_duplicates.coordsort.bam | while read file; do
#ls MAPPED/PWKLY.toPWK.marked_duplicates.coordsort.bam | while read file; do
#ls MAPPED/PWKLY.toLEWES.marked_duplicates.coordsort.bam | while read file; do
	samp_name=$(echo "${file}" | cut -d "/" -f 2 | cut -d "." -f 1-2)
	if [ -f "${samp_name}.vcf.gz" ]; then
		echo "${samp_name}.vcf.gz exists; skipping..."
	else
		echo "Processing sample: ${samp_name}"
		echo "Calling variants..."
		#gatk --java-options "-Xmx4g" HaplotypeCaller -R ${ref_fa} -I ${file} -O "${samp_name}.g.vcf.gz" -ERC GVCF --native-pair-hmm-threads 24 #--sample-name "${samp_name}"
		#gatk --java-options "-Xmx4g" GenotypeGVCFs -R ${ref_fa} -V "${samp_name}.g.vcf.gz" -O "${samp_name}.vcf.gz"
		gatk --java-options "-Xmx4g" HaplotypeCaller -R ${ref_fa} -I ${file} -O "${samp_name}.vcf.gz" --native-pair-hmm-threads 64
		echo "Filtering..."
		bcftools filter -m+ -e 'MQ < 30.0 || FORMAT/DP < 3 || FORMAT/DP > 30 || ALT="*"' -s+ --IndelGap 5 -Ov -o "${samp_name}.filter.vcf" "${samp_name}.vcf.gz"
		bcftools view -i 'FILTER="PASS"' -Ov -o  "${samp_name}.hardFilter.vcf" "${samp_name}.filter.vcf"
	fi
done

echo "Done!" 

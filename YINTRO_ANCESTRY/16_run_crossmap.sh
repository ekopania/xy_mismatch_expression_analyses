#!/bin/bash
##PURPOSE: Run CrossMap to Liftover VCF coordinates to match GRCm38 mouse reference
#
# Job name:
#SBATCH --job-name=crossMap
#SBATCH --output=crossMap-%j.log
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=emily.kopania@umconnect.umt.edu # Where to send mail
#SBATCH --cpus-per-task=1 # Number of cores per MPI rank (ie number of threads, I think)
#SBATCH --nodes=1 #Number of nodes
#SBATCH --ntasks=1 # Number of MPI ranks (ie number of processes, I think)
#SBATCH --mem-per-cpu=40G #Not sure if I should mess with these...
#SBATCH --mem=96000 #Not sure if I should mess with these...
# Partition:
## Since you want to run it on 72 cores, the partition good_cpu has nodes with 72 cores.
#SBATCH --partition=good_lab_reincarnation
##SBATCH -w, --nodelist=compute-0-4 # run on a specific node
#
## Command(s) to run:
source ~/software/anaconda/anaconda3/bin/activate
conda activate ek_main_enviro

#ulimit -n 8024

in_vcf="PWKLY.toLEWES.hardFilter.vcf"
out_vcf="PWKLY.toLEWES.hardFilter.lifted_over.vcf"
chain_file="REFS/REF-to-LEWES_EiJ.formatted.chain"
ref="REFS/Mus_musculus.GRCm38.dna.primary_assembly.fa"

CrossMap.py vcf ${chain_file}  ${in_vcf}  ${ref}  ${out_vcf} --no-comp-alleles  

#--no-comp-alleles: If set, CrossMap does NOT check if the reference allele is different from the alt allele
#If this is not set, many variants will be lost b/c they are not variants compared to GRCm39
#This option is why I switched to CrossMap instead of picard LiftoverVcf

echo "Done!"

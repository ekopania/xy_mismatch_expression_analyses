#!/bin/bash
##PURPOSE: Run Picard's Liftover VCF to get coordinates to match GRCm38 mouse reference
#
# Job name:
#SBATCH --job-name=liftoverVCF
#SBATCH --output=liftoverVCF-%j.log
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

in_vcf="PWKLY.toPWK.hardFilter.vcf"
out_vcf="PWKLY.toPWK.hardFilter.lifted_over.vcf"
chain_file="REFS/REF-to-PWK_PhJ.formatted.chain"
reject_vcf="PWKLY.toPWK.hardFilter.liftOver_rejects.vcf"
ref="REFS/Mus_musculus.GRCm38.dna.primary_assembly.fa"

picard -Xmx40g LiftoverVcf I="${in_vcf}" O="${out_vcf}" CHAIN="${chain_file}" REJECT="${reject_vcf}" R="${ref}" LIFTOVER_MIN_MATCH=0 WARN_ON_MISSING_CONTIG=true RECOVER_SWAPPED_REF_ALT=true

echo "Done!"

#!/bin/bash
#PURPOSE: Use g2gtools to generate a chain file from VCF indels (relative to mouse reference GRCm38)
#https://g2gtools.readthedocs.io/en/latest/usage.html
# Job name:
#SBATCH --job-name=vcf2chain
#SBATCH --output=vcf2chain-%j.log
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=emily.kopania@umconnect.umt.edu # Where to send mail
#SBATCH --cpus-per-task=1 # Number of cores per MPI rank (ie number of threads, I think)
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
source ~/software/anaconda/anaconda3/bin/activate
conda activate g2gtools

REF=/mnt/beegfs/ek112884/cis_trans/MODFILES/Mus_musculus.GRCm38.dna.primary_assembly.fa
VCF_INDELS=/mnt/beegfs/ek112884/cis_trans/MODFILES/PWK_PhJ.mgp.v5.indels.dbSNP142.normed.vcf.gz
STRAIN=PWK_PhJ

echo "Generating chain file for ${STRAIN}..."

g2gtools vcf2chain -f ${REF} -i ${VCF_INDELS} -s ${STRAIN} -o REFS/REF-to-${STRAIN}.chain

echo "Done!"

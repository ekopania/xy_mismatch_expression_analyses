#!/bin/bash
#PURPOSE: g2gtools vcf2chain outputs chain file with a couple extra columns; this sed command just keeps the first 3, which is the official format for a UCSC chain file that picard tools needs to be happy
#https://genome.ucsc.edu/goldenPath/help/chain.html
# Job name:
#SBATCH --job-name=format_chain
#SBATCH --output=format_chain-%j.log
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
conda activate ek_main_enviro

in_chain=REFS/REF-to-PWK_PhJ.chain
out_chain=REFS/REF-to-PWK_PhJ.formatted.chain

sed -e 's/\t[A-Z].*$//' -e 's/MT$/22/' -e 's/X$/20/' -e 's/Y$/21/' "${in_chain}" > "${out_chain}"

#NOTE: chain ID has to be numeric; the last 3 -e commands conver MT, X, and Y to numeric ($ should make it the end of the line only, thus changing just the chain ID but not the chromosome names in other parts of the header)
#MT -> 22
#X -> 20
#Y -> 21

echo "Done!"

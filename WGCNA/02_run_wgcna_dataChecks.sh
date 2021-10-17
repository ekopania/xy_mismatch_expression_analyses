#!/bin/bash
#PURPOSE: Run WGCNA data checks
#
# Job name:
#SBATCH --job-name=WGCNA_dataChecks
#SBATCH --output=WGCNA_dataChecks.exp2_RS.log
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
##SBATCH --mail-user=ekopania4@gmail.com # Where to send mail
#SBATCH --cpus-per-task=1 # Number of cores per MPI rank (ie number of threads, I think)
#SBATCH --nodes=1 #Number of nodes
#SBATCH --ntasks=1 # Number of MPI ranks (ie number of processes, I think)
#SBATCH --mem-per-cpu=8G #Not sure if I should mess with these...
#SBATCH --mem=256G #Not sure if I should mess with these...
# Partition:
## Since you want to run it on 72 cores, the partition good_cpu has nodes with 72 cores.
#SBATCH --partition=good_lab_cpu
##SBATCH -w, --nodelist=compute-0-4 # run on a specific node
#
## Command(s) to run:

#All data (NOT recommended as experiments 1 and 2 were done on different sorters
#Rscript 03_WGCNA_dataChecks.r all

#Experiment 1 (XY mismatch on non-hybrid background)
#Rscript 03_WGCNA_dataChecks.r exp1
#Rscript 03_WGCNA_dataChecks.r exp1_LZ
#Rscript 03_WGCNA_dataChecks.r exp1_RS

#Experiment 2 (XY match on hybrid background)
#Rscript 03_WGCNA_dataChecks.r exp2
#Rscript 03_WGCNA_dataChecks.r exp2_LZ
Rscript 03_WGCNA_dataChecks.r exp2_RS

echo "Done!"

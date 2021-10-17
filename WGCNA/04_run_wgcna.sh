#!/bin/bash
#PURPOSE: Run WGCNA
#
# Job name:
#SBATCH --job-name=WGCNA
#SBATCH --output=WGCNA.exp2_RS.log
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
#signed
#Rscript 05_WGCNA.r all TRUE
#unsigned
#Rscript 05_WGCNA.r all FALSE

#Experiment 1 (XY mismatch on non-hybrid background)
#signed
#Rscript 05_WGCNA.r exp1 TRUE
#Rscript 05_WGCNA.r exp1_LZ TRUE
#Rscript 05_WGCNA.r exp1_RS TRUE
#unsigned
#Rscript 05_WGCNA.r exp1 FALSE

#Experiment 2 (XY match on hybrid background)
#signed
#Rscript 05_WGCNA.r exp2 TRUE
#Rscript 05_WGCNA.r exp2_LZ TRUE
Rscript 05_WGCNA.r exp2_RS TRUE
#unsigned
#Rscript 05_WGCNA.r exp2 FALSE

echo "Done!"

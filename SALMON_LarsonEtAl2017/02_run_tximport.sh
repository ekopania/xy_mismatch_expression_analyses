#!/bin/bash
#PURPOSE: Run tximport to convert salmon output (transcripts) to gene level counts
#
# Job name:
#SBATCH --job-name=salmon_tximport
#SBATCH --output=tximport-%j.log
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
##SBATCH --mail-user=ekopania4@gmail.com # Where to send mail
#SBATCH --cpus-per-task=1 # Number of cores per MPI rank (ie number of threads, I think)
#SBATCH --nodes=1 #Number of nodes
#SBATCH --ntasks=1 # Number of MPI ranks (ie number of processes, I think)
#SBATCH --mem-per-cpu=8G #Not sure if I should mess with these...
#SBATCH --mem=32G #Not sure if I should mess with these...
# Partition:
## Since you want to run it on 72 cores, the partition good_cpu has nodes with 72 cores.
#SBATCH --partition=good_lab_reincarnation
##SBATCH -w, --nodelist=compute-0-4 # run on a specific node
#
## Command(s) to run:
source ~/software/anaconda/anaconda3/bin/activate
conda activate r4
Rscript 03_tximport.r
echo "Done!"

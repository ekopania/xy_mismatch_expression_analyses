#!/bin/bash
#PURPOSE: Run salmon index to generate salmon index for salmon quantification
#
# Job name:
#SBATCH --job-name=salmon_index
#SBATCH --output=salmon_index.log
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
##SBATCH --mail-user=ekopania4@gmail.com # Where to send mail
#SBATCH --cpus-per-task=1 # Number of cores per MPI rank (ie number of threads, I think)
#SBATCH --nodes=1 #Number of nodes
#SBATCH --ntasks=1 # Number of MPI ranks (ie number of processes, I think)
#SBATCH --mem-per-cpu=8G #Not sure if I should mess with these...
#SBATCH --mem=0 #Not sure if I should mess with these...
# Partition:
## Since you want to run it on 72 cores, the partition good_cpu has nodes with 72 cores.
#SBATCH --partition=good_lab_reincarnation
##SBATCH -w, --nodelist=compute-0-4 # run on a specific node
#
## Command(s) to run:
echo "Building salmon index..."
salmon index -t /mnt/beegfs/ek112884/REFERENCE_DIR/Mus_musculus.GRCm38.cdna.all.fa.gz -i salmon_index.GRCm38.idx /mnt/beegfs/ek112884/REFERENCE_DIR/Mus_musculus.GRCm38.cdna.all.fa.gz
#Salmon tutorial uses an Ensembl reference transcriptome for indexing (ex: *.cdna.all.fa.gz)
echo "Done!"

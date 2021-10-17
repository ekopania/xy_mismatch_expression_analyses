#!/bin/bash
#PURPOSE: Merge and log transform normalized expression data for input into WGCNA; make table of "traits" (cross type and cell type)
#
# Job name:
#SBATCH --job-name=WGCNA_preprocess
#SBATCH --output=WGCNA_preprocess.log
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
##SBATCH --mail-user=ekopania4@gmail.com # Where to send mail
#SBATCH --cpus-per-task=1 # Number of cores per MPI rank (ie number of threads, I think)
#SBATCH --nodes=1 #Number of nodes
#SBATCH --ntasks=1 # Number of MPI ranks (ie number of processes, I think)
#SBATCH --mem-per-cpu=4G #Not sure if I should mess with these...
#SBATCH --mem=16G #Not sure if I should mess with these...
# Partition:
## Since you want to run it on 72 cores, the partition good_cpu has nodes with 72 cores.
#SBATCH --partition=good_lab_cpu
##SBATCH -w, --nodelist=compute-0-4 # run on a specific node
#
## Command(s) to run:

#paste /mnt/beegfs/ek112884/amplicon_expression_analysis/SALMON/fpkm_full_table.Yintro_exp1.LZ.95.txt <(cut -d"	" -f2- /mnt/beegfs/ek112884/amplicon_expression_analysis/SALMON/fpkm_full_table.Yintro_exp1.RS.95.txt) <(cut -d"	" -f2- /mnt/beegfs/ek112884/amplicon_expression_analysis/SALMON/fpkm_full_table.Yintro_exp2.LZ.95.txt) <(cut -d"	" -f2- /mnt/beegfs/ek112884/amplicon_expression_analysis/SALMON/fpkm_full_table.Yintro_exp2.RS.95.txt) > fpkm_full_table.txt

#join /mnt/beegfs/ek112884/amplicon_expression_analysis/SALMON/fpkm_full_table.Yintro_exp1.LZ.95.txt /mnt/beegfs/ek112884/amplicon_expression_analysis/SALMON/fpkm_full_table.Yintro_exp1.RS.95.txt > temp1.txt
#join temp1.txt /mnt/beegfs/ek112884/amplicon_expression_analysis/SALMON/fpkm_full_table.Yintro_exp2.LZ.95.txt > temp2.txt
#join temp2.txt /mnt/beegfs/ek112884/amplicon_expression_analysis/SALMON/fpkm_full_table.Yintro_exp2.RS.95.txt > fpkm_full_table.txt
#rm temp1.txt
#rm temp2.txt

Rscript 01_process_exp_data.r

echo "Done!"

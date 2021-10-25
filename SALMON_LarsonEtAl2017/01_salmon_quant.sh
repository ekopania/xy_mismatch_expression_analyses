#!/bin/bash
#PURPOSE: Run salmon quant for alignment-independent RNAseq quantification
#
# Job name:
#SBATCH --job-name=salmon_quant_IU
#SBATCH --output=salmon_quant.libTypeIU-%j.log
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
##SBATCH --mail-user=ekopania4@gmail.com # Where to send mail
#SBATCH --cpus-per-task=8 # Number of cores per MPI rank (ie number of threads, I think)
#SBATCH --nodes=1 #Number of nodes
#SBATCH --ntasks=1 # Number of MPI ranks (ie number of processes, I think)
#SBATCH --mem-per-cpu=32G #Not sure if I should mess with these...
#SBATCH --mem=0 #Not sure if I should mess with these...
# Partition:
## Since you want to run it on 72 cores, the partition good_cpu has nodes with 72 cores.
#SBATCH --partition=good_lab_cpu
##SBATCH -w, --nodelist=compute-0-4 # run on a specific node
#
## Command(s) to run:
source ~/software/anaconda/anaconda3/bin/activate
conda activate ek_main_enviro
echo "Running salmon quant..."

ls /mnt/beegfs/ek112884/amplicon_expression_analysis/DATA_LarsonEtAl2017/RAW_DATA/merged/*1P.fq.gz | while read file; do #Loop through all trimmed read 1 fastq files (unaligned sequence reads, post-trimmomatic) 
	sample=$(echo "${file}" | cut -d "/" -f 9 | cut -d "_" -f 1-4)
	echo "${sample}"
	long_name=$(echo "${file}" | cut -d "_" -f 1-8)
	file2="${long_name}_2P.fq.gz"
	echo "${file2}"
	#fileU="${long_name}_1U.fq.gz" #include unpaired forward reads
	#echo "${fileU}"
	#NOTE: --validateMappings deprecated in version 1.5 and higher; selective alignment mode (which this flag used to invoke) is now the default
        #NOTE: cannot combine PE and SE data (don't use -r flag with -1 and -2 flags)
        #PE reads, ISR library type, use pre-made mm10 index from salmon's refgenie site
        #http://refgenomes.databio.org/v3/assets/splash/0f10d83b1050c08dd53189986f60970b92a315aa7a16a6f1/salmon_sa_index?tag=default
        #salmon quant -p 8 -i /mnt/beegfs/ek112884/REFERENCE_DIR/refgenie_genome_folder/alias/mm10/salmon_sa_index/default/ -l ISR -1 "${file}" -2 "${file2}" -o "${sample}.salmon_output_libTypeISR_premadeIndex"
	#PE reads, ISF library type, use pre-made mm10 index from salmon's refgenie site
	#salmon quant -p 8 -i /mnt/beegfs/ek112884/REFERENCE_DIR/refgenie_genome_folder/alias/mm10/salmon_sa_index/default/ -l ISF -1 "${file}" -2 "${file2}" -o "${sample}.salmon_output_libTypeISF_premadeIndex"
	#PE reads, IU (unstranded) library type, use pre-made mm10 index from salmon's refgenie site
	salmon quant -p 8 -i /mnt/beegfs/ek112884/REFERENCE_DIR/refgenie_genome_folder/alias/mm10/salmon_sa_index/default/ -l IU -1 "${file}" -2 "${file2}" -o "${sample}.salmon_output_libTypeIU_premadeIndex"
	#PE reads - automatic library type detection
	#salmon quant -p 8 -i salmon_index.GRCm38.idx -l ISR -1 "${file}" -2 "${file2}" --validateMappings -o "${sample}.salmon_output_libTypeISR"
	#PE reads
	#-l flag is library type; ISF means inwards, stranded, read1 is forward
	#More info on salmon lib types: https://salmon.readthedocs.io/en/latest/library_type.html#fraglibtype
#	salmon quant -p 8 -i salmon_index.GRCm38.idx -l ISF -1 "${file}" -2 "${file2}" --validateMappings -o "${sample}.salmon_output_PE"
	#SE reads or 1U forward unpaired reads output from trimmomatic
	#Note libtype change to SF for stranded forward (no I for inward b/c no mate)
#	salmon quant -p 8 -i salmon_index.GRCm38.idx -l SF -r "${fileU}" --validateMappings -o "${sample}.salmon_output_SE"
	#Merge PE and SE quantifications; you need to do this separately by field for some reason
	#salmon quantmerge --quants "${sample}.salmon_output_PE" "${sample}.salmon_output_SE" --column=numreads -o "${sample}.salmon_output.quantmerge.numreads.txt" #Raw txpts quantification
	#salmon quantmerge --quants "${sample}.salmon_output_PE" "${sample}.salmon_output_SE" --column=len -o "${sample}.salmon_output.quantmerge.len.txt" #Txpt length
	#salmon quantmerge --quants "${sample}.salmon_output_PE" "${sample}.salmon_output_SE" --column=elen -o "${sample}.salmon_output.quantmerge.elen.txt" #txpt effective length
	#salmon quantmerge --quants "${sample}.salmon_output_PE" "${sample}.salmon_output_SE" -c --column=tpm -o "${sample}.salmon_output.quantmerge.tpm.txt" #TPM
done

echo "Done!"

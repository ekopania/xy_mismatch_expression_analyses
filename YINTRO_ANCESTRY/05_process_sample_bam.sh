#!/bin/bash
##PURPOSE: Process sample bam: Reorder to match ref if necessary, query sort, fix mates, coord sort, mark dups
#
# Job name:
#SBATCH --job-name=bam_processing
#SBATCH --output=bam_procesing-%j.log
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
ulimit -n 8024
#ls MAPPED/LLLL.toLEWES.reheader.bam | while read file; do
#ls MAPPED/LLLL.toWSB.reheader.bam | while read file; do
ls MAPPED/LEWPY.toLEWES.reheader.bam | while read file; do
#ls MAPPED/LEWPY.toPWK.reheader.bam | while read file; do
#ls MAPPED/PWKLY.toLEWES.reheader.bam | while read file; do
#ls MAPPED/PWKLY.toPWK.reheader.bam | while read file; do
        name=$(echo "${file}" | cut -d "." -f 1-2)
        echo "${name}"
	short_name=$(echo "${name}" | cut -d "/" -f 2)
	echo "${short_name}"
	#echo PROCESSED_BAMS/*/${short_name}.marked_duplicates.coordsort.bam
#	if [ -f PROCESSED_BAMS/WDIS_RESEQ_PROCESSED/${short_name}.marked_duplicates.coordsort.bam ]; then
#		echo "Bam exists; skipping..."
#	else
		echo "Processing bam..."
	#	picard -Xmx40g BuildBamIndex I="${file}"
	#	picard -Xmx40g ReorderSam I="${file}" O="${name}.reordered.bam" R="REFS/Mus_musculus.GRCm38.dna.primary_assembly.fa" ALLOW_INCOMPLETE_DICT_CONCORDANCE=true TMP_DIR=PICARD_TEMP_WDIS/ #Only run if sam header not in same chr/contig order as reference fasta
	#	picard -Xmx40g SortSam I="${name}.reordered.bam" O="${name}.querysorted.bam" SORT_ORDER=queryname TMP_DIR=PICARD_TEMP/
		picard -Xmx40g SortSam I="${file}" O="${name}.querysorted.bam" SORT_ORDER=queryname TMP_DIR=PICARD_TEMP/ #VALIDATION_STRINGENCY=SILENT
        	picard -Xmx40g FixMateInformation I="${name}.querysorted.bam" O="${name}.fixed_mate.bam" TMP_DIR=PICARD_TEMP/ #VALIDATION_STRINGENCY=SILENT
		picard -Xmx40g SortSam I="${name}.fixed_mate.bam" O="${name}.coordsort.bam" SORT_ORDER=coordinate TMP_DIR=PICARD_TEMP/ #VALIDATION_STRINGENCY=SILENT
		picard -Xmx40g MarkDuplicates I="${name}.coordsort.bam" O="${name}.marked_duplicates.coordsort.bam" M="${name}.dup_metrics.txt" MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=8000 TMP_DIR=PICARD_TEMP/ #VALIDATION_STRINGENCY=SILENT
	#	picard -Xmx40g MarkDuplicates I="${name}.fixed_mate.bam" O="${name}.marked_duplicates.TEST.bam" M="${name}.dup_metrics.TEST.txt" MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=8000 TMP_DIR=PICARD_TEMP_WDIS/ VALIDATION_STRINGENCY=SILENT
		#picard -Xmx40g SortSam I="${name}.marked_duplicates.TEST.bam" O="${name}.marked_duplicates.coordsort.TEST.bam" SORT_ORDER=coordinate TMP_DIR=PICARD_TEMP_WDIS/ VALIDATION_STRINGENCY=SILENT
        	#picard -Xmx40g CollectWgsMetrics I="${name}.marked_duplicates.coordsort.bam" O="${name}.wgs_metrics.txt" R="/mnt/beegfs/ek112884/REFERENCE_DIR/GRCm39/GCA_000001635.9_GRCm39_genomic.fna" TMP_DIR=PICARD_TEMP_WDIS/ VALIDATION_STRINGENCY=SILENT
	#        java -jar ~/software/picard.jar CollectAlignmentSummaryMetrics R="REFS/Mus_musculus.GRCm38.dna.primary_assembly.fa" I="${name}.marked_duplicates.coordsorted.bam" O="${name}.align_summary_metrics.txt" TMP_DIR=PICARD_WDIS/
	#fi
done

echo "Done!"

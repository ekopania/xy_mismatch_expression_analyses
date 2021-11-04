#NOT A REAL SCRIPT
#Just keeping track of what I did from the command line:

samtools view -H MAPPED/LLLL.toLEWES.bam > temp.header.sam
vim temp.header.sam
#Manually add these lines at the bottom:
#@RG     ID:LLLL_CSFP200005444-1a_H3NMWDSXY_L1.toLEWES   LB:lib1 PL:illumina     PU:unit1        SM:LLLL_CSFP200005444-1a_H3NMWDSXY_L1.toLEWES
#@RG     ID:LLLL_CSFP200005444-1a_H3LK3DSXY_L2.toLEWES   LB:lib1 PL:illumina     PU:unit1        SM:LLLL_CSFP200005444-1a_H3LK3DSXY_L2.toLEWES
#@RG     ID:LLLL_RESEQ.toLEWES   LB:lib1 PL:illumina     PU:unit1        SM:LLLL_RESEQ.toLEWES
#@RG     ID:LLLL_CSFP200005444-1a_H3NMWDSXY_L1.1U.toLEWES   LB:lib1 PL:illumina     PU:unit1        SM:LLLL_CSFP200005444-1a_H3NMWDSXY_L1.1U.toLEWES
#@RG     ID:LLLL_CSFP200005444-1a_H3LK3DSXY_L2.1U.toLEWES   LB:lib1 PL:illumina     PU:unit1        SM:LLLL_CSFP200005444-1a_H3LK3DSXY_L2.1U.toLEWES
#@RG     ID:LLLL_RESEQ.1U.toLEWES   LB:lib1 PL:illumina     PU:unit1        SM:LLLL_RESEQ.1U.toLEWES
samtools reheader temp.header.sam MAPPED/LLLL.toLEWES.bam > MAPPED/LLLL.toLEWES.reheader.bam
samtools index MAPPED/LLLL.toLEWES.reheader.bam
rm temp.header.sam

#Repeat for toWSB
samtools view -H MAPPED/LLLL.toWSB.bam >temp.header.toWSB.sam
vim temp.header.toWSB.sam
#Manually add these lines at the bottom:
#@RG	ID:LLLL_CSFP200005444-1a_H3NMWDSXY_L1.toWSB LB:lib1	PL:illumina	PU:unit1	SM:LLLL_CSFP200005444-1a_H3NMWDSXY_L1.toWSB
#@RG	ID:LLLL_CSFP200005444-1a_H3LK3DSXY_L2.toWSB	LB:lib1	PL:illumina	PU:unit1	SM:LLLL_CSFP200005444-1a_H3LK3DSXY_L2.toWSB
#@RG	ID:LLLL_RESEQ.toWSB	LB:lib1	PL:illumina	PU:unit1	SM:LLLL_RESEQ.toWSB
#@RG    ID:LLLL_CSFP200005444-1a_H3NMWDSXY_L1.1U.toWSB LB:lib1     PL:illumina     PU:unit1        SM:LLLL_CSFP200005444-1a_H3NMWDSXY_L1.1U.toWSB
#@RG    ID:LLLL_CSFP200005444-1a_H3LK3DSXY_L2.1U.toWSB     LB:lib1 PL:illumina     PU:unit1        SM:LLLL_CSFP200005444-1a_H3LK3DSXY_L2.1U.toWSB
#@RG    ID:LLLL_RESEQ.1U.toWSB     LB:lib1 PL:illumina     PU:unit1        SM:LLLL_RESEQ.1U.toWSB
samtools reheader temp.header.toWSB.sam MAPPED/LLLL.toWSB.bam > MAPPED/LLLL.toWSB.reheader.bam
samtools index MAPPED/LLLL.toWSB.reheader.bam
rm temp.header.toWSB.sam

#Readgroups For LEWPY
#I realized I was doing the sample tag wrong (SM); should be the same for all since it's all the same sample
#Repleace toLEWES with toPWK for mappings to PWK
#@RG    ID:LEWPY_RESEQ.toLEWES  LB:lib1     PL:illumina     PU:unit1        SM:LEWPY
#@RG    ID:LEWPY_CSFP200005448-1a_H3LK3DSXY_L2.toLEWES  LB:lib1     PL:illumina     PU:unit1        SM:LEWPY
#@RG    ID:LEWPY_RESEQ.1U.toLEWES  LB:lib1     PL:illumina     PU:unit1        SM:LEWPY
#@RG    ID:LEWPY_CSFP200005448-1a_H3LK3DSXY_L2.1U.toLEWES  LB:lib1     PL:illumina     PU:unit1        SM:LEWPY

#Readgroups for PWKLY
#Repleace toLEWES with toPWK for mappings to PWK
#@RG    ID:PWKLY_CSFP200005447-1a_H3LK3DSXY_L2.toLEWES	LB:lib1     PL:illumina     PU:unit1        SM:PWKLY
#@RG    ID:PWKLY_RESEQ.toLEWES	LB:lib1     PL:illumina     PU:unit1        SM:PWKLY
#@RG    ID:PWKLY_RESEQ.1U.toLEWES	LB:lib1     PL:illumina     PU:unit1        SM:PWKLY
#@RG    ID:PWKLY_CSFP200005447-1a_H3LK3DSXY_L2.1U.toLEWES	LB:lib1     PL:illumina     PU:unit1        SM:PWKLY

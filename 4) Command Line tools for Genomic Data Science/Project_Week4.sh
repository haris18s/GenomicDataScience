#!/usr/bin/env

#create a bowtie2 index, with prefix athal

bowtie2-build -f gencommand_proj4/athal_chr.fa indeces/athal

#q1)
samtools view tophat/8days/accepted_hits.bam |wc -l

#q2)
samtools view tophat/16days/accepted_hits.bam |wc -l

#q3-6)
#Check the align_summary.txt 
cat tophat/16days/align_summary.txt

#q7-8)
#spliced alignments--> check the accepted_hits.bam and find the cigar strings that contain "N"

#9-10)
#check unmapped


#q11 -20 use of cufflinks

#q11, 12, number of genes identified from cufflinks
cut -f 9 transcripts.gtf| cut  -d " " -f 2 | sort -u| wc -l


#q13,14
cut -f 9 transcripts.gtf| cut  -d " " -f 4 | sort -u| wc -l


#15,16
#all trnascripts -->\
cut -f3,9  transcripts.gtf|cut -f1,4 -d " "|grep transcript|wc -l

#find the genes that have multiple trnacripts .2 or .3 
cut -f3,9  transcripts.gtf|cut -f1,4 -d " "|grep transcript|grep '.2"'


#21-30 runn coffcompare 

#run cuffcompare in the cuffcomp dir, output goes in cufflinks dir
#17-20, Single and multi exon transcripts.
cat transcripts.gtf | cut -f 3,9 | cut -f 1| awk '{ if ($1=="transcript") count=0; count++; } count>1 {print "Transcript with multiple  exons"  }'


cat transcripts.gtf | awk 'BEGIN {transcriptCount=0; exonCount=0; currentTranscript="";} { if ($3=="exon") { exonCount++; } if ($3=="transcript") { if (exonCount > 1) transcriptCount++; exonCount=0; currentTranscript=$9; } } END { if (exonCount > 1) transcriptCount++; print "Transcripts with multiple exons:", transcriptCount; }'


awk -F'\t' '$3 == "exon" {transcript[$9]++} END {count=0; for (id in transcript) {if (transcript[id] == 1) count++} print "Total single-exon genes: " count}' transcripts.gtf



#q21-22
cat cuffcmp.transcripts.gtf.tmap |tail -n +2|wc -l

#q23-24
#grep AT4G20240 cuffcmp.transcripts.gtf.tmap

#q25-26 
cut -f 3 cuffcmp.transcripts.gtf.tmap |tail -n +2 | grep -c "c".


#q27-28
# 'j'  it represents a novel splce variant. Potentially novel isoform.
cut -f 3 cuffcmp.transcripts.gtf.tmap | grep  -c "j"



#q29-30
cut -f 3 cuffcmp.transcripts.gtf.tmap | tail -n +2| grep  -c "e"

#q31
cat merged.gtf |cut -d " " -f 2 |sort -u |wc -l

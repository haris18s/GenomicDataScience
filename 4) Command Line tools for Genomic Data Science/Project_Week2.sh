#!/usr/bin/env

#How many alignments does the set contain?
samtools flagstats Arabidopsis.bam


#How many alignments show the read’s mate unmapped?

#flag for unmapped mate, which is the fourth flag
samtools view -f 4 file.bam |cut -f1 |wc -l

#How many alignments contain a deletion (D)?
#find the cigar string check for D à deletion
samtools view atha.bam | cut -f6 | grep “D” 



#How many alignments show the read’s mate mapped to the same chromosome?
#find the column with mate chr, =, chr num or * if it is the same, different chr or non chr found, respectively that has mapped the mate alignment. 
samtools view atha.bam | cut -f7 | grep “=” |wc -l

#How many alignments are spliced?
#find the CIGAR string and check if it is a N ( skipped region – intron)
samtools view atha.bam | cut – f6 | grep “N”





#For a range Chr3:11,777,000-11,794,000”.
#Create a sorted alignment file with sort function from the samtools 
samtools sort athal_wu_0_A.bam athal_wu_0_A.sorted
#index this file
samtools index athal_wu_0_A.sorted.bam

#Extract alignment within this specific range, -b will create bam file
samtools view –b athal_wu_0_A.sorted.bam “Chr3:11777000-11794000”> athal_wu_0_A.region.bam


#Question 6 -..
#6. How many alignments does the set contain?
samtools view athal_wu.bam “Chr3:11777000-11794000"
#7. How many alignments show the read’s mate unmapped?
Samtools view atha_wu.bam “Chr3:11777000-11794000" | cut –f7 | grep “*” |wc –l

#8. How many alignments contain a deletion (D)?
Samtools view atha_wu.bam “Chr3:11777000-11794000" | cut –f6 | grep “D” |wc –l

#9. How many alignments show the read’s mate mapped to the same chromosome?
Samtools view atha_wu.bam “Chr3:11777000-11794000" | cut –f7 | grep “=” |wc –l

#10. How many alignments are spliced?
Samtools view atha_wu.bam “Chr3:11777000-11794000" | cut –f7 | grep “N” |wc –l

#11. How many sequences are in the genome file?
#Information found in the head of the bam file. Count the number of lines describing the sequences in the reference genome. The lines with SN.
samtools view –H athal_wu_0_A.bam | grep –c “SN:”

#12. What is the length of the first sequence in the genome file?
#It is the “LN” (denotes the length)after the “SN”.
samtools view –H athal_wu_0_A.bam | grep “SN:” | more

#13. What alignment tool was used?
#Denoted by the “PG”
samtools view –H athal_wu_0_A.bam | grep “^@PG”

#14. What is the read identifier (name) for the first alignment?
samtools view athal_wu_0_A.bam | head -1 | cut –f1

#15. What is the start position of this read’s mate on the genome? Give this as ‘chrom:pos’ if the read was mapped, or ‘*” if unmapped.
#Check the chr: and postion of the reads’mate.
samtools view athal_wu_0_A.bam | head -1 | cut –f1

#16. How many overlaps (each overlap is reported on one line) are reported?
#1st exatract  the alignment with the specific locus 
samtools view athal_wu_0_A.bam “Chr3:11777000-11794000">locus.bam
#Convert the athal_wu_0_A.bam into bed file
bedtools bamtobed –wo –a locus.bam -b athal_wu_0_A.bam


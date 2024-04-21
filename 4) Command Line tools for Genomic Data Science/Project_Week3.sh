#!/usr/bin/env
#1.How many sequences were in the genome?
grep - “>”  wu_0.v7.fas

#2-3. What was the name of the third/last  sequence in the genome file? Give the name only, without the “>” sign.
grep - “>”  wu_0.v7.fas | grep –n 3 --> the line 3, n indeixates the number



#4.How many index files did the operation create?
#create an index file
bowtie2-build –f wu__0_A.v7.fas indeces/wu_0


ls –1 /indeces |wc –l --> count number of lines

#5. What is the 3-character extension for the index files created?
#Bjust check the extension .bt2 from bowtie2.

#align the reads to the genome with the index file
bowtie2 –p 4 –x wu_A.v7.fas wu_0_A_wgs.fastq -S wu.bt2.sam

#convert to bam file
Samtools view –bT wu_0.v7.fas(reference_genome) wu.bt2.sam >wu_A.bt2.bam



#view the file
samtools view wu_A.bt2.bam

#6. How many reads were in the original fastq file?
#Each fastq recrd is represented by 4 lines in the wu_0)wgs.fastq. Count lines  divide by 4.

#7. How many matches (alignments) were reported for the original (full-match) setting?
#Exclude lines in the file containing unmapped reads.


bowtie2 –x wu_0/wu_0 –U wu_0_A_wgs.fastq –S wu_0.bt2.sam
bowtie2 –x wu_0/wu_0 –U wu_0_A_wgs.fastq –S wu_0_local.bt2.sam --local


cat wu_0.bt2.sam | grep -v "@"| cut -f 3 |grep -v "*" |wc -l.


#8. How many matches (alignments) were reported with the local-match setting? Exclude lines in the file containing unmapped reads.
cat wu_0_local.bt2.sam | grep -v "@"| cut -f 3 |grep -v "*" |wc -l.



#11.How many reads had multiple matches in the scenario in Question 7? You can find this in the bowtie2 summary; note that by default bowtie2 only reports the best match for each read.
#check bowtie summary


#12.How many reads had multiple matches in the scenario in Question 8? Use the format above. You can find this in the bowtie2 summary; note that by default bowtie2 only reports the best match for each read.
#check bowtie summary

#13. How many alignments contained insertions and/or deletions, in the scenario in Question 7?
cut -f6 |grep -v "*" |grep -E '[ID]' --color|wc -l


#14. How many alignments contained insertions and/or deletions, in the scenario in Question 8?
cut -f6 |grep -v "*" |grep -E '[ID]' --color|wc -l

#15. How many entries were reported for Chr3?

#Step 1: Start by converting the SAM file to BAM format as indicated, then sorting it:
samtools view –bT wu_0.v7.fas wu_0.bt2.sam > wu_A.bt2.bam
#then sorting it:
#It is critical to sort the BAM file before analyzing it with samtools.
samtools sort wu_A.bt2.bam -o wu_A_sorted.bt2.bam

#index the file
samtools index wu_A_sorted.bt2.bam

#run samtools mpileup to call variations, sammtools version=1.3
samtools mpileup -vu -f wu_0.v7.fas wu_A_sorted.bt2.bam >wu_variations.mpileup.vcf


#16. How many entries have ‘A’ as the corresponding genome letter?
cat wu_variations.mpileup.vcf |grep -v "#" |cut -f4|grep -P "^A$" |wc -l

#17. How many entries have exactly 20 supporting reads (read depth)?
cat wu_variations.mpileup.vcf |grep -v "#" |grep "DP=20"|wc -l

#18. How many entries represent indels?
cat wu_variations.mpileup.vcf |grep -v "#" |grep "INDEL"|wc -l


#19. How many entries are reported for position 175672 on Chr1?

cat wu_variations.mpileup.vcf |grep -v "#" |grep "175672"



#20. How many variants are called on Chr3?

#Step 1: First re-run ‘SAMtools mpileup’ with the BCF output option ‘-g’:
samtools mpileup -g -u -f wu_0.v7.fas wu_0_sorted.bt2.bam >wu_variations.mpileup.bcf

#then call variants using ‘BCFtools call’ with the multi-allelic caller (option ‘-m’), showing only variant sites (‘-v’) and presenting the output in uncompressed VCF format (‘-O v’), as instructed:
bcftools call --multiallelic-caller -v -O v wu_variations.mpileup.bcf >wu_variations_final.vcf

#Step 2: To answer the question, we count all reported variants that show ‘Chr3’ in column 1:

cat wu_variations_final.vcf  |grep –v “^#” |  cut –f1 | sort | uniq –c | grep “Chr3”



#21.How many variants represent an A->T SNP? If useful, you can use ‘grep –P’ to allow tabular spaces in the search term.

more wu_variations_final.vcf |grep -v "^#"|cut -f4,5|grep -P "^A\tT$"|wc -l

#22. How many entries are indels?
more wu_variations_final.vcf |grep -v "^#"|grep  "INDEL"|wc -l

#23.How many entries have precisely 20 supporting reads (read depth)?

cat wu_variations_final.vcf |grep -v "^#"|grep  'DP=20'

#24. What type of variant (i.e., SNP or INDEL) is called at position 11937923 on Chr3?
cat wu_variations_final.vcf |grep -v "^#"|grep  11937923


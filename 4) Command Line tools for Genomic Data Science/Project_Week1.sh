#!/usr/bin/env

#How many chromosomes are there in the genome?
grep “>” apple.genome |wc -l


#number of genes
cut -f1 apple.genes

#number of transcript variants (find the distinct number of distinct names in column2)
cut -f2 apple.genes| sort -u |wc-l

#or second method
cut -f2 apple.genes| uniq |wc-l



#How many genes have a single splice variant? Find the number of lines (genes)  that have a single line
cut -f1 apple.genes |sort |uniq -c |grep “ 1 “ |wc -l

#hoow many genes have a single splice variant?
cut –f1 apple.genes | sort | uniq –c | grep –c “ 1 “



#How may genes have 2 or more splice variants?
cut -f1 apple.genes |sort |uniq -c |grep “ 2 “ |wc -l

#How many genes are there on the ‘+’ strand?
cut -f1,4 apple.genes |uniq -c | grep + |wc -l

#How many genes are there on the ‘-’ strand?
cut -f1,4 apple.genes |uniq -c | grep + |wc -l

#How many genes are there on chromosome chr1? #relevant info in column 1,3
cut -f1,3 apple.genes|sort -u|cut -f2 |sort | uniq-c
#2nd way
cut -f1,3, apple.genes |sort -u |grep chr1 |wc -l



#How many genes are there on chromosome chr2? #relevant info in column 1,3
cut -f1,3, apple.genes |sort -u |grep chr2 |wc -l

#How many transcripts are there on chr1?

cut –f2,3 apple.genes | cut –f2 | sort | uniq –c (solution way)



#How many genes are in common between condition A and condition B?
cut -f1 apple.conditionA |sort -u >sortedA
cut -f1 apple.conditionB |sort -u >sortedB
comm -1 -2 sorted A sorted B

#How many genes are specific to condition B ?
comm -3 -1 sorted A sorted


#How many genes are in common to all three conditions?
cut -f1 apple.conditionC |sort -u >sortedC
comm -1 -2 sorted A sorted B >commAB
comm -1 -2 commAB sortedC


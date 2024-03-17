library(ALL)
library(biomaRt)
library(GenomicRanges)
library(airway)
library(AnnotationHub)
data(ALL)
ALL
#Question 1L  What is the mean expression across all features for sample 5 in the ALL dataset (from the ALL package)?
#get the actual data
sample_5 = mean(exprs(ALL)[,5])
mart <- useMart(host= 'feb2014.archive.ensembl.org', biomart = "ENSEMBL_MART_ENSEMBL")


#Question: Using this version of Ensembl, annotate each feature of the ALL dataset with the Ensembl gene id. 
#How many probesets (features) are annotated with more than one Ensembl gene id?


#check datasets
listMarts(mart)
#pick a dataset
ensemble_dataset = useDataset("hsapiens_gene_ensembl",mart)

#get the feature names of ALL dataset
values_probes = featureNames(ALL)

#check the ALL data, annotation based on hgu95av2 thats why this used in attribute and filter, 
#affy_hg_u95av2 this used in this version of ensembl

genes_overl = getBM(attributes =c("ensembl_gene_id", "affy_hg_u95av2"),
      filters = "affy_hg_u95av2", values = values_probes, mart = ensemble_dataset)


dup_probsets = sum(table(genes_overl[,2])>1)


#Question 3: How many probesets (Affymetrix IDs) are annotated with one or more genes on the autosomes (chromosomes 1 to 22).
#check the attributes that you can pass when search for the database -use chromosome name 
listAttributes(ensemble_dataset)

genes_overl_chr = getBM(attributes =c("ensembl_gene_id", "affy_hg_u95av2", "chromosome_name"),
                    filters = "affy_hg_u95av2", values = values_probes, mart = ensemble_dataset)

#filter the genes located in autosomes
genes_in_autosomal = genes_overl_chr[genes_overl_chr$chromosome_name %in% c(1:22),]

length(unique(genes_in_autosomal$affy_hg_u95av2))

#1st way: get with unique function probsets associated with genes in autosomes
num_probesets <- length(unique(genes_in_autosomal$affy_hg_u95av2))

#2nd way: subset the gene found in autosomes with probset column
num_probsets = sum(table(genes_in_autosomal[,2])>=1)


#Question4 - Clarification  4:Use the MsetEx dataset from the minfiData package. 
#Part of this question is to use the help system to figure out how to address the question.

#Question 4: What is the mean value of the Methylation channel across the features for sample “5723646052_R04C01”?
library(minfiData)

data(MsetEx.sub)
#get the assays names for the Data, in order to extract the data  
assayNames(MsetEx)
#extract data with assay function and the name of the assay
meth_data = assay(MsetEx, "Meth")
#mean for sample 5723646052_R04C01
mean(meth_data[,'5723646052_R04C01'])



#Question 5: Access the processed data from NCBI GEO Accession number GSE788. What is the mean expression level of sample GSM9024?
library(GEOquery)
eList = getGEO("GSE788")

#check the elist
names(eList)
length(eList)
#get the data
eData = eList[[1]]


#get the expression data for the sample GSM9024
data_sample = (exprs(eData)[,"GSM9024"])
mean(data_sample)


#Question 6: What is the average of the average length across the samples in the expriment?
#We are using the airway dataset from the airway package.
library(airway)
data(airway)
#get the medata for the samples
colData(airway)
#mean of avg length of the reads 
avg_length =colData(airway)$avgLength
mean(avg_length)


#Question 7: What is the number of Ensembl genes which have a count of 1 read or more in sample SRR1039512?
#We are using the airway dataset from the airway package. The features in this dataset are Ensembl genes.

#find the num of genes that has read >=1 for sample SRR1039512
num_of_reads_SRR= assay(airway, "counts")[,'SRR1039512']

genes_with_counts_gr_1 = length(num_of_reads_SRR[num_of_reads_SRR>=1])


#Question 8 The airway dataset contains more than 64k features. 
#How many of these features overlaps with transcripts on the autosomes (chromosomes 1-22),
#as represented by the TxDb.Hsapiens.UCSC.hg19.knownGene package?

#Clarification: A feature has to overlap the actual transcript, not the intron of a transcript. 
#So you will need to make sure that the transcript representation does not contain introns.
library(TxDb.Hsapiens.UCSC.hg19.knownGene)


rownames(assay(airway))
#get the tdxb data
txdb = TxDb.Hsapiens.UCSC.hg19.knownGene

#find all the airway features
airway_features_gr = rowRanges(airway)
#get the exon - and group them by gene
exons_tsc = exonsBy(txdb, "tx")

range(subsetByOverlaps(airway_features_gr, autosomal_exons))

#keep exons in autosomal chrosmomes
autosomes_exons = keepSeqlevels(exons_tsc, paste0("chr", c(1:22)), pruning.mode = "coarse")

#rename to be based on the 
new_style = mapSeqlevels(seqlevels(autosomes_exons), "NCBI")
autosomal_exons = renameSeqlevels(autosomes_exons,new_style)


overlapping_features <- subsetByOverlaps(airway_features_gr, autosomal_exons)
length(overlapping_features)


'#Question 9
The expression measures of the airway dataset are the number of reads mapping to each feature. 
In the previous question we have established that many of these features do not overlap autosomal transcripts 
from the TxDb.Hsapiens.UCSC.hg19.knownGene. But how many reads map to features which overlaps these transcripts?

#Question: For sample SRR1039508, how big a percentage (expressed as a number between 0 and 1) of the total reads
in the airway dataset for that sample, are part of a feature which overlaps an autosomal TxDb.Hsapiens.UCSC.hg19.knownGene 
transcript?'

data(airway)
#get the names for the features that overlap with transcripts in ghe autosomal
names_ovrl_trans =  names(overlapping_features)
#get the count data for the sample SRR1039508
data_sample = assay(airway, 'counts')[,1]
#subset the data counts of the sample
counts_overlapping = data_sample[names_ovrl_trans]

#get the total number of reads for the counts_overlaping with transcripts
total_reads_overlapping = sum(counts_overlapping)
#get the total number of reads for the sample SRR1039508
total_number_reads_SRR1039508 = sum(data_sample)

#get the percentage of the tota reads in the airway datset of the feautes that overlap
total_reads_overlapping/total_number_reads_SRR1039508



#Question 10:
exons_tsc = exonsBy(txdb, "tx")

#obtain H3K4me3 narrowPeaks from the E096 sample using Annotation hub
ahub = AnnotationHub()
#make a query
qhs  = query(ahub, c("H3K4me3", 'E096'))
data_narrow_methylat  = qhs[['AH30596']]

#keep autosomal chromosomes
data_narrow_autosomal = keepSeqlevels(data_narrow_methylat, paste0("chr", c(1:22)), pruning.mode = "coarse")
#rename the chrs of the data_narrow peak in order to compare with Granges from exons 
new_style = mapSeqlevels(seqlevels(data_narrow_autosomal), "NCBI")
data_narrow_autosomal = renameSeqlevels(data_narrow_autosomal, new_style)


#get the range from overlapping features (Q8 - the airway features overlapped with exons)
#ranges will combine all the intervals in each list object, NOTE: if you do not do range aftet getting range_overlappinfeatures you get smaller median
range_orverlapping_features = range(overlapping_features)
#get their promoters of this range
promoters_overlaping_range = promoters(range_orverlapping_features)

#subset the promoters with the data_narrow methylated (H3K4me3)
promoters_contain_narrowpeak = subsetByOverlaps(promoters_overlaping_range, data_narrow_autosomal)
#get the median number of counts for the promoters containing a H3K4me
median(data_sample[names(promoters_contain_narrowpeak)])

       
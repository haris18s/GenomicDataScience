#packages
library(yeastRNASeq)
library(ShortRead)
library(Biostrings)
library(leeBamViews)
library(GenomicRanges)
#Question 1:Question 1
'The yeastRNASeq experiment data package contains FASTQ files from an RNA seq experiment in yeast.
When the package is installed, you can access one of the FASTQ files by the path given by:'

#access the fastq filez
fastqFilePath <- system.file("reads", "wt_1_f.fastq.gz", package = "yeastRNASeq")
#read the fastq file
fastq_yeast_reads <- readFastq(fastqFilePath)
#get the reads
reads_yeast= sread(fastq_yeast_reads)

#get nucleotides for each position with consensusmatrix
consensus_yeast = consensusMatrix(reads_yeast)
#get fraction of reads that have A for the 5h position
A_in_5th =  consensus_yeast[,5][1]
total_reads = sum(consensusMatrix(a)[,5])
A_in_5th/total_reads


#Question 2: What is the average numeric quality value of the 5th base of these reads?
qualityReadsYeast = quality(fastq_yeast_reads)
#convert qualities to Phred scores and extract 5h position
score5thBase = as(qualityReadsYeast, 'matrix')[,5]
#get the mean
mean(score5thBase)



'Question 3: .
The leeBamViews  experiment data package contains aligned BAM files from an RNA seq experiment in yeast
(the same experiment as in Questions 1 and 2, but that is not pertinent to the question). 
You can access one of the BAM files by the path given by'

bamFilePath <- system.file("bam", "isowt5_13e.bam", package="leeBamViews")

bamFile <- BamFile(bamFilePath)

#ead bam file
aln <- scanBam(bamFile)
aln
aln <- aln[[1]]
names(aln)
aln$rname
#specify the region of interest
gr = GRanges(seqnames = "Scchr13", ranges = IRanges(start =800000 , end = 801000 ))
gr
params <- ScanBamParam(which = gr, what = scanBamWhat())
#return a character vecotr indicating which parts are available for extraction
scanBamWhat()
#read the bam file using the params
aln <- scanBam(bamFile, param=params)

#get the num of reads in each position, 
freq_table = table(aln[[1]]$pos)
#into a data frame
freq_df = as.data.frame(freq_table)
values_greater_than_1 <- freq_df$Freq[freq_df$Freq > 1]
#sum the num of reads
sum(values_greater_than_1)



'Question 4

This is a continuation of Question 3.

The package contains 8 BAM files in total, representing 8 different samples from 4 groups. A full list of file paths can be had as'

bpaths_q4 <- list.files(system.file("bam", package="leeBamViews"), pattern = "bam$", full=TRUE)
bpaths_q4
'An objective of the original paper was the discovery of novel transcribed regions in yeast. One such region is Scchr13:807762-808068.
Question: What is the average number of reads across the 8 samples falling in this interval?'

gr_q4 = GRanges(seqnames = "Scchr13", ranges = IRanges(start =807762 , end = 808068))
params_q4 <- ScanBamParam(which = gr_q4, what = scanBamWhat())

#make a loop to loop each file and calculate the av
num_samples = length(bpaths_q4)
total_reads_in_interval = 0
for (y in bpaths_q4) {
  file_path_perFile = y
  bamFile_perPath <- BamFile(file_path_perFile)
  aln_perBam =  scanBam(bamFile_perPath, param = params_q4)
  reads_perFile_interval = length(aln_perBam$`Scchr13:807762-808068`$seq)
  
  total_reads_in_interval = total_reads_in_interval + reads_perFile_interval 
  
}
avg_reads_inteval = total_reads_in_interval/num_samples



'In the lecture on the oligo package an ExpressionSet with 18 samples is constructed, 
representing normalized data from an Affymetrix gene expression microarray. 
The samples are divided into two groups given by the groupgroup variable.'

'Question 5: What is the average expression across samples in the control group for the “8149273” probeset 
(this is a character identifier, not a row number).'
library(oligo)


list.files("GSE38792")
list.files()

celfiles <- list.files('GSE38792/CEL', full = TRUE)
celfiles

#read the files
rawData <- read.celfiles(celfiles)
#file names
filename <- sampleNames(rawData)
sampleNames <- sub(".*_", "", filename)
sampleNames <- sub(".CEL.gz$", "", sampleNames)

#replace the sample nams of raw data
sampleNames(rawData) <- sampleNames

pData(rawData)$group <- ifelse(grepl("^OSA", sampleNames(rawData)),
                               "OSA", "Control")
#normalize the data
normData <- rma(rawData)

#exprs data for the 8149273 probeset
exprs_data <- exprs(normData)
mean(exprs_data['8149273',1:8])



'This is a continuation of Question 5.

Use the limma package to fit a two group comparison between the control group and the OSA group,
and borrow strength across the genes using eBayes()eBayes(). Include all 18 samples in the model fit.

Question 6: What is the absolute value of the log foldchange (logFClogFC) of the gene with the lowest P.valueP.value.'
library(limma)

design = model.matrix(~normData$group)
fit <- lmFit(normData, design)
fit <- eBayes(fit)

#get fold change and find the lowest p value, inf extract all the rows, t can be found evern without inf, toptable gives the top-ranked genes.
result_table = topTable(fit, number = Inf)
result_table_sorted = result_table[order(result_table$P.Value), ]
result_table_sorted[1,]

'#This is a continuation of Question 6.

Question 7: How many genes are differentially expressed between the two groups at an adj.P.valueadj.P.value cutoff of 0.05?'

topTable(fit, p.value = 0.05)





'Question 8

An example 450k dataset is contained in the minfiData

 package. This dataset contains 6 samples; 3 cancer and 3 normals. Cancer has been shown to be globally hypo-methylated (less methylated) compared to normal tissue of the same kind.

Take the RGsetEx dataset in this package and preprocess it with the preprocessFunnorm function. For each sample, compute the average Beta value (percent methylation) across so-called OpenSea loci.

Question: What is the mean difference in beta values between the 3 normal samples and the 3 cancer samples, across OpenSea CpGs?'


library(minfiData)
data(RGsetEx)

pData(RGsetEx)
#preprocess the data
procesRGset = preprocessFunnorm(RGsetEx)
procesRGset[getIslandStatus(procesRGset) == 'OpenSea']

mean(getBeta(procesRGset[getIslandStatus(procesRGset) == 'OpenSea'])[,c(1)])


sampleNames(procesRGset)
getIslandStatus(procesRGset) == 'OpenSea'


normal_samples = getBeta(procesRGset)[,c(1,2,5)]
cancer_samples = getBeta(procesRGset)[,c(3,4,6)]
open_seaIslands = getIslandStatus(procesRGset) == 'OpenSea'


diffrence =  mean(normal_samples[opean_seaIslands]) - mean(cancer_samples[opean_seaIslands])




'
Question 9
This is a continuation of Question 8.

The Caco2 cell line is a colon cancer cell line profiled by ENCODE. Obtain the narrowPeak DNase hyper sensitive sites computed by the analysis working group (AWG).

Question: How many of these DNase hypersensitive sites contain one or more CpGs on the 450k array?'
library(AnnotationHub)
#find the encode data for the Caco2 cell line 
ahub = AnnotationHub()
qhs = query(ahub, c('Caco2','DNase' ))
Caco_encode_meth = qhs[[1]]
#get the granes of 450K microarray
gr450K = granges(procesRGset)

subsetByOverlaps(Caco_encode_meth, gr450K)



'Question 10

The zebrafishRNASeq  package contains summarized data from an RNA-seq experiment in zebrafish in the form of a data.frame called zfGeneszfGenes. The experiment compared 3 control samples to 3 treatment samples.

Each row is a transcript; the data.frame contains 92 rows with spikein transcripts; these have a rowname starting with “ERCC”. Exclude these rows from the analysis.

Use DESeq2

 to perform a differential expression analysis between control and treatment. Do not discard (filter) genes and use the padjpadj results output as the p-value.

Question: How many features are differentially expressed between control and treatment (ie. padj <= 0.05padj <= 0.05)?'

library(DESeq2)
library(zebrafishRNASeq)
data(zfGenes)

#find and exclude the rows with spike in transcripts
zFgenes_withoutSpike = zfGenes[grep("^(?!ERCC)", rownames(zfGenes), perl = TRUE),]


as.data.frame(zFgenes_withoutSpike)
#create col data
col_data = data.frame(condition = c(rep("Ctrl",3), rep("Trt",3)), row.names = names(zFgenes_withoutSpike))
#create a container
dds = DESeqDataSetFromMatrix(countData = zFgenes_withoutSpike, 
                                        colData = col_data,
                                        design = ~ condition)
#run DESEq
dds <- DESeq(dds)

res_deseq <- results(dds)

res_df = as.data.frame(res_deseq)

filter_df1 = res_df[res_df$padj <= 0.05,]
col_data
dim(na.omit(filter_df2))






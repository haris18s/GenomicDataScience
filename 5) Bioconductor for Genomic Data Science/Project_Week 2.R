#Project Week 2 -----

#Question 1:What is the GC content of “chr22” in the “hg19” build of the human genome?
#Tip: The reference genome includes “N” bases; you will need to exclude those.
library(BSgenome)
library(AnnotationHub)
library("BSgenome.Hsapiens.UCSC.hg19.masked")
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(GenomicRanges)

available.genomes()


hg_19 = BSgenome.Hsapiens.UCSC.hg19.masked
hg_19_chr22 = hg_19$chr22

GC_chr22 = letterFrequency(hg_19_chr22, "GC", as.prob = TRUE) - letterFrequency(hg_19_chr22, "N", as.prob = TRUE)


#Question 2: Question: What is mean GC content of H3K27me3 “narrowPeak” regions from Epigenomics Roadmap from the H1 stem cell line on chr 22.

ahub = AnnotationHub()

qhs = query(ahub, c("H3K27me3", "H1", "E003"))
qhs

EPi_data_H3K27= qhs[[2]]

Epi_chr22_narrow = keepSeqlevels(EPi_data_H3K27, "chr22", pruning.mode = "coarse")



#get the nt sequence for each position in the human genome for chr 22
chr22_sequence = Views(hg_19, Epi_chr22_narrow)
average_GC_chr22 = mean(letterFrequency(chr22_sequence, "GC", as.prob = TRUE))

#Question 3: Question: What is the correlation between GC content and “signalValue” of these regions (on chr22)?
signal_EpiH3 = Epi_chr22_narrow$signalValue
average_signal_chr22 = mean(Epi_chr22_narrow$signalValue)
cor_GC_signal = cor(letterFrequency(chr22_sequence, "GC", as.prob = TRUE), Epi_chr22_narrow$signalValue)



#Question 4:Question: what is the correlation between the “signalValue” of the “narrowPeak” regions and the average “fc.signal” 
#across the same regions?
#Clarification: First compute the average “fc.signal” for across each region, for example using “Views”; this yields a single number of each region.
#Next correlate these numbers with the “signalValue” of ithe “narrowPeaks”.
qhs_fc = query(ahub, c("H3K27me3", "H1", "E003", "fc.signal"))
file_fc = qhs_fc[['AH32033']]

#load the data 
#bigwig_file <- "C:/Users/Haris/AppData/Local/R/cache/R/AnnotationHub/934268b6ca3_37473"

data_fc_ch22_E003 <- import(file_fc, format = "BigWig",
                  which = GRanges("chr22", IRanges(start=1, end = 10^8)),
                  as = "Rle")$chr22



#NOTE: for fc_signal each base has a signal
#find the regions that overlap between the "narrow" regions and fc_data set
overlap_regions_narrow_fc = Views(data_fc_ch22_E003, 
      start = start(Epi_chr22_narrow),
      end = end(Epi_chr22_narrow))


#find the mean in 
mean_fc_in_overlap =mean(overlap_regions_narrow_fc)

cor_mean_fc_signal = cor(mean_fc_in_overlap, signal_EpiH3)

#Question 5:Question: How many bases on chr22 have an fc.signal greater than or equal to 1?\
#rle object

#
signal_greater_1 = slice(data_fc_ch22_E003,1)
total_bases_great_1 = sum(width(signal_greater_1)) 




'
The H1 stem cell line is an embryonic stem cell line, a so-called pluripotent cell. 
Many epigenetic marks change upon differentiation. We will examine this. 
We choose the cell type with Roadmap ID “E055” which is foreskin fibroblast primary cells.
We will use the “fc.signal” for this cell type for the H3K27me3 mark, on chr22. 
We now have a signal track for E003 and a signal track for E055.
We want to identify regions of the genome which gain H3K27me3 upon differentiation. 
These are regions which have a higher signal in E055 than in E003. 
To do this properly, we would need to standardize (normalize) the signal across the two samples; we will ignore this for now.'

#Question 6: Identify the regions of the genome where the signal in E003 is 0.5 or lower and the signal in E055 is 2 or higher.

qhs_E055 = query(ahub, c("E055", "fc.signal", "H3K27me3"))
data_E055 = qhs_E055[["AH32470"]]
data_file_E055_chr22 = import(data_E055, format = "BigWig",
                        which = GRanges('chr22', IRanges(start=1, end = 10^8)),
                        as = "Rle")$chr22



E055_values_gr_2 = slice(data_file_E055_chr22,lower =2)
E003_less_th_05 = slice(data_fc_ch22_E003, upper=0.5)

v1 = as(E055_values_gr_2, "IRanges")
v2 = as(E003_less_th_05, "IRanges")
regions_overlap = sum((width(intersect(v1,v2))))




#Question 7 :What is the average observed-to-expected ratio of CpG dinucleotides for CpG Islands on chromosome 22?
Cpg_islands_query =  query(ahub, c("CpG", "hg19"))
Cpg_islands = Cpg_islands_query[[1]]

#Cpg for chr 22
Cpg_chr22 = keepSeqlevels(Cpg_islands, "chr22", pruning.mode = "coarse")



seqs_Cpg = Views(hg_19, Cpg_chr22)
observed_dinucl = dinucleotideFrequency(seqs_Cpg)


#loop to calculate ration per CpG island
sum_ratio = 0
for (i in seq_along(seqs_Cpg)) { # seq_along generates sequences of incs
  individual_cpg <- seqs_Cpg[i]
  observed_CG_per_CpG = dinucleotideFrequency(individual_cpg)[7]
  
  G_in_Cpg = (letterFrequency(individual_cpg,"G"))
  C_in_Cpg = (letterFrequency(individual_cpg,"C"))
  expected_per_cpg = G_in_Cpg * C_in_Cpg / width(individual_cpg)
  
  ratio_per_CpG = observed_CG_per_CpG / expected_per_cpg
  
  sum_ratio = sum_ratio + ratio_per_CpG
}


final_ratio_CpGs = sum_ratio/length(seqs_Cpg)


#Question 8:
#A TATA box is a DNA element of the form “TATAAA”. Around 25% of genes should have a TATA box in their promoter. 
#We will examine this statement.

#Question: How many TATA boxes are there on chr 22 of build hg19 of the human genome?
#Clarification: You need to remember to search both forward and reverse strands.

TATA = DNAString("TATAAA")
reverse_comp_TATa = reverseComplement(TATA)


match_forward = countPattern(TATA, hg_19_chr22)
match_reverse = countPattern(reverse_comp_TATa, hg_19_chr22)

total_matches = match_forward + match_reverse



#Question 9 :Question: How many promoters of transcripts on chromosome 22 containing a coding sequence, 
#contains a TATA box on the same strand as the transcript?

#Clarification: Use the TxDb.Hsapiens.UCSC.hg19.knownGene package to define transcripts and coding sequence. 
#Here, we defined a promoter to be 900bp upstream and 100bp downstream of the transcription start site.

txdb = TxDb.Hsapiens.UCSC.hg19.knownGene 
txdb

gr = GRanges(seqnames = 'chr22',  ranges = IRanges(start = start(hg_19_chr22), end = end(hg_19_chr22)))

#get the cds per transcript
cds_by_tx = unlist(cdsBy(txdb, by = "tx"))
#keep for chr22 and get their names
names_cds_by_tx_chr22 = names(keepSeqlevels(cds_by_tsc, "chr22", pruning.mode = "coarse"))

#take all the transcript and their names to subset 
all_tx <- transcripts(txdb, columns = "tx_id")

#subset the found transcripts from all trnascripts
subset_tx <- all_tx[all_tx$tx_id %in% names_cds_by_tx_chr22, ]
#find their promoters
subset_promoters <- promoters(subset_tx, upstream = 900, downstream = 100)
#find the matches of TATA box in chr22 of whole genome 
matches_chr22 = keepSeqlevels(vmatchPattern(TATA, hg_19), "chr22", pruning.mode = "coarse")

#overlap the promoters of transcripts and the matches found in chr22 that contain TATA

prom_contain_TATA = subsetByOverlaps(subset_promoters, matches_chr22)

#very nice function not uses
matches_gr <- as((matches_chr22_chr22), "GRanges")



#Question 10:
"
It is possible for two promoters from different transcripts to overlap, 
in which case the regulatory features inside the overlap might affect both transcripts. This happens frequently in bacteria.

Question: How many bases on chr22 are part of more than one promoter of a coding sequence?

Clarification: Use the TxDb.Hsapiens.UCSC.hg19.knownGene package to define transcripts and coding sequence.
Here, we define a promoter to be 900bp upstream and 100bp downstream of the transcription start site. In this case, 
ignore strand in the analysis."


#check the coverage of bases of the promoters
cov_bases_prom = coverage(subset_promoters)
#take the coverage for chr 22
cov_chr22 = cov_bases_prom[22]
#check how many bases have value >1 and sum 

bases_overlap = sum(cov_chr22>1)



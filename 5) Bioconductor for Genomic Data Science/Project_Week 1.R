#Project Week 1 -----

#Question 1:Use the AnnotationHub package to obtain data on "CpG Islands" in the human genome.

#Question: How many islands exists on the autosomes?

library(AnnotationHub)

ah = AnnotationHub()
ah$genome
#subset
ah = subset(ah, species == "Homo sapiens")

qhs = query(ah, "CpG Islands")
qhs

gr1 = qhs[[1]]
gr1

#keep only autosomal chroms- not working
keepStandardChromosomes(gr1, pruning.mode = "coarse")

#keepSeqlevels function working, first creating a vector of chrs
nums = c(1:22)
chromosomes = paste0("chr", nums)
CPG_islands  = keepSeqlevels(gr1, chromosomes, pruning.mode = "coarse")


#Question 2: How many islands exist on chromosome 4.
keepSeqlevels(gr1, "chr4", pruning.mode = "coarse")



#Question 3:  Obtain H3K4me3 for the H1 cell line from Epigenmics Roadmaps, keep regions only mapped to 1-22
#Question: How many bases does these regions cover?
qhs2 = query(ah, c("H3K4me3", "H1", "E003"))

qhs2


Epi_H3_E003 = qhs2[[2]]

#for autosomal chromosomes
Epi_H3K4_E003_chr1_22 = keepSeqlevels(Epi_H3_E003, chromosomes, pruning.mode = "coarse")
#check the chromosomes
seqlevels(Epi_H3K4_E003_chr1_22)


#bases in H3K4 modification
bases_in_H3K4 = sum(width(reduce(Epi_H3K4_E003_chr1_22)))

#reduce function reduces the num of columns, width gives the ranges for each genomic region and sum the positions
#number of bases in CPG
bases_in_CPG = sum(width(reduce(CPG_islands)))




#Question 4: mean signal value for histone modificaiton H3K27me3, H1 cell line

qhs3 = query(ah, c("H3K27me3", "H1", "E003"))
Epi_H327 = qhs3[[2]]
Epi_H327_chr1_22 = keepSeqlevels(Epi_H327, chromosomes, pruning.mode = "coarse")

mean(Epi_H327_chr1_22$signalValue)


#Question 5: Find Bivelant regions from these above histone marks
#we will say a region is bivalent if it is enriched in both H3K4me3 and H3K27me3

bivelant_regions = intersect(Epi_H3K4_E003_chr1_22, Epi_H327_chr1_22, ignore.strand =TRUE)
bases_in_bivelant = sum(width(bivelant_regions))


#Question 6: We will examine the extent to which bivalent regions overlap CpG Islands.
#Question 6: how big a fraction (expressed as a number between 0 and 1) of the bivalent regions, overlap one or more CpG Islands

num_biv_ov_CPG = sum(countOverlaps(bivelant_regions,CPG_islands))

fraction_biv_ov_CPG = num_biv_ov_CPG/ length(bivelant_regions)


"Question 7: Question: How big a fraction (expressed as a number between 0 and 1) of the bases which are part of CpG Islands,
are also bivalent marked."

#find common regions between CPG and bivelant with interest

gr_CPG_biv = intersect(CPG_islands, bivelant_regions, ignore.strand = TRUE)
bases_in_CPG_biv = sum(width(reduce(gr_CPG_biv)))

per_bases_in_CPG_biv = bases_in_CPG_biv/bases_in_CPG


"Question 8: How many bases are bivalently marked within 10kb of CpG Islands?"
#Consider resize function, to change the coordinates of the IRanges, change by adding 10kbp to width(of each irange)
?resize
res



KBp_CPG = resize(CPG_islands, width = width(CPG_islands+ 10000), ignore.strand = FALSE, fix = "center")

gr_CPG_biv_10kbp = intersect(KBp_CPG, bivelant_regions)

bases_biv_10KBp = sum(width(reduce(gr_CPG_biv_10kbp)))


"Question 9: How big a fraction (expressed as a number between 0 and 1) of the human genome is contained in a CpG Islands"
#numerato the part that is beig considered
#denominator: the whole or the total number
genome_size = sum(seqlengths(CPG_islands))

size_CPG_islands = sum(width(CPG_islands))

fraction_genome_in_CPGs = size_CPG_islands/genome_size


"Question 10 : Compute an odds-ratio for the overlap of bivalent marks with CpG islands."


odds_matrix  = matrix(data = 0, nrow = 2, ncol=2)
colnames(odds_matrix) = c("in", "out")
rownames(odds_matrix) = c("in", "out")
odds_matrix



#in -in the bases that are in both 
genome = sum(seqlengths(CPG_islands))

odds_matrix[1,1] = sum(width(intersect(CPG_islands,bivelant_regions, ignore.strand = TRUE)))


#in- out bases in CPS but not in bivelants
odds_matrix[1,2] =sum(width(setdiff(CPG_islands,bivelant_regions , ignore.strand =TRUE)))
#out -in bases not in CPGs but in bivelant
odds_matrix[2,1] = sum(width(setdiff(bivelant_regions, CPG_islands, ignore.strand =TRUE)))


#out-out is the human genome
odds_matrix[2,2] =  genome_size - sum(odds_matrix)
odds_ratio = odds_matrix[1,1] * odds_matrix[2,2] / (odds_matrix[2,1] * odds_matrix[1,2] )

odds_ratio


#!/usr/bin/env python

"""kmer_index.py: A k-mer index for indexing a text."""


import bisect
from bm_preproc import BoyerMoore


def naive_22m(pattern, text):
    """Similar to naive algorithm but it takes up to 2 mismatches acceptable"""

    occurances_2mm = []
    for i in range(len(text) -len(pattern) +1):
        count_mismatch = 0
        for j in range(len(pattern)):
            if not text[i+j] == pattern[j]:
                count_mismatch +=1
        if count_mismatch <3:
            occurances_2mm.append(i)
    return occurances_2mm


def read_reference_genome(filename):
    """ Description: This function reads a reference genome from a file and returns it as a string."""
    genome = ''
    with open(filename) as fh:
        for line in fh:
            if not line.startswith(">"):
                genome +=line.rstrip()
    return genome


class Index(object):
    """ This class holds a substring index for a given text t. It creates an index from all substrings of t of a specified length k"""

    def __init__(self, t, k):
        """ Create index from all substrings of t of length k """
        self.k = k  # k-mer length (k)
        self.index = []
        for i in range(len(t) - k + 1):  # for each k-mer
            self.index.append((t[i:i+k], i))  # add (k-mer, offset) pair
        self.index.sort()  # alphabetize by k-mer
    def query(self, p):
        """ Return index hits for first k-mer of p """
         # query with first k-mer
        i = bisect.bisect_left(self.index, (p, -1))  # binary search
        hits = []
        while i < len(self.index):  # collect matching index entries
            if self.index[i][0] != p:
                break
            hits.append(self.index[i][1])
            i += 1
        return hits


def partition_match_index(p, t, n):
    """THis functions is similar to approximate matching of pigenonhole principle
    DifferenceL now use pigeonole principle to match each partition with  list of index function
     question 4
     Question 5: index hits, how many times the segment matched the kme from index list"""
    # n = maximial number of mismatches that allowed

    segment_length = int(round(len(p) / (n + 1)))
    all_matches = set()  # fill the set with indeces that you found matches
    # loop into segments in p
    index_hits = 0
    for i in range(n + 1):
        start = i * segment_length
        end = min((i + 1) * segment_length,
                  len(p))  # min make sure that does not pass the end of the string, bcs might not be evengly divisible
        # make the preprocessing for good suffix rules and good characer rule
        index = Index(t,8)
        matches = index.query(p[start:end])
        #the number of times that a segment match the list of kmers
        index_hits += len(matches)
        # step on these positions to make sure that rest of p matches t with no more than n mismastches
        for m in matches:
            # if either this is true p run of the beginning or pass the end of t, so skip with continue
            if m < start or m - start + len(p) > len(t):
                continue

            mismatches = 0
            # this part compare the part of p  before the segment that we already compared
            for j in range(0, start):
                if not p[j] == t[m - start + j]:
                    mismatches += 1
                    if mismatches > n:
                        break

            # now compare the suffix after that segment
            for j in range(end, len(p)):
                if not p[j] == t[m - start + j]:
                    mismatches += 1
                    if mismatches > n:
                        break

            if mismatches <= n:
                # all matches is a st so if we have 4 partitions each partition will match 4 times, by having a set only 1 offset time wll be added
                all_matches.add(m - start)  # add the beginning, so m - start
    return list(all_matches),index_hits


"""question 6 - - Subsequence Indexing, For homework2 in coursera python for Algorithms for dna sequencing"""

"""example:https://nbviewer.org/github/BenLangmead/ads1-hw-examples/blob/master/hw2_query_subseq_index.ipynb"""

import bisect

class SubseqIndex(object):
    """ Holds a subsequence index for a text T """

    def __init__(self, t, k, ival):
        """ Create index from all subsequences consisting of k characters
            spaced ival positions apart.  E.g., SubseqIndex("ATAT", 2, 2)
            extracts ("AA", 0) and ("TT", 1). """
        self.k = k  # num characters per subsequence extracted
        self.ival = ival  # space between them; 1=adjacent, 2=every other, etc
        self.index = []
        self.span = 1 + ival * (k - 1)
        for i in range(len(t) - self.span + 1):  # for each subseq
            self.index.append((t[i:i + self.span:ival], i))  # add (subseq, offset)
        self.index.sort()  # alphabetize by subseq

    def query(self, p):
        """ Return index hits for first subseq of p """
        subseq = p[:self.span:self.ival]  # query with first subseq
        i = bisect.bisect_left(self.index, (subseq, -1))  # binary search
        hits = []
        while i < len(self.index):  # collect matching index entries
            if self.index[i][0] != subseq:
                break
            hits.append(self.index[i][1])
            i += 1
        return hits


def query_subseq(p,subseq_ind):
    """ Calculates the numberr of hits for given pattern(p) using subsequence index
    Input :
    - p (str): The pattern to search for.
    - subseq_ind (SubseqIndex): The subsequence index object used for querying.
    Returns:
        int: The total number of hits found for the pattern p in the subsequence index.
    """
    num_hits = [len(subseq_ind.query(p[i:])) for i in range(3)]

    return sum(num_hits)


if __name__ == "__main__":
    #For question4 and 5
    # Read reference genome from file
    genome = read_reference_genome("chr1.GRCh38.excerpt.fasta")
    pattern = "GGCGCGGTGGCTCACGCCTGTAAT"
    # index

    # Perform approximate matching using partitioning approach
    all_matches, index_hits = partition_match_index(pattern, genome, 2)
    print(len(all_matches))
    print(index_hits)

    # Perform approximate matching using naive approach
    print(f"naive{len(naive_22m(pattern, genome))} ")



    #Next code For question 6

    t = 'to-morrow and to-morrow and to-morrow creeps in this petty pace'
    p = 'to-morrow and to-morrow '


    pattern = "GGCGCGGTGGCTCACGCCTGTAAT"
    genome =read_reference_genome("chr1.GRCh38.excerpt.fasta")

    subseq_ind = SubseqIndex(genome, 8, 3)



    print(query_subseq(pattern,subseq_ind))


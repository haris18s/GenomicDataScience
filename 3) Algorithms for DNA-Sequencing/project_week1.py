"""
Project week 1 : Algorithms for DNA sequencing
Course: Genomic data science.



"""
def naive(pattern,text):
    """This function searches for occurances of a specific pattern within a given text with the native exact algorithm
    It compares a given pattern with all substrings of the text of the same length

    Input: PAttern(str): the pattern to search
    -Text (str) the text in which occurances are searched.
    Returns: a list containing the indiceas of all occurances."""

    occurances = []

    for i in range(len(text)-len(pattern) +1):
        match = True
        for j in range(len(pattern)):
            if not text[i+j] == pattern[j]:
                match = False
                break
        if match:
            occurances.append(i)
    return occurances


def reverse_complement(seq):
    """Teturns the reverse complement of a sequence
    Input: - seq (str) : The DNA sequence to find the rev.complement.
    Returns :St, The reverse commplement of the input DNA sequence.

    """

    complement = {"A": "T", "C":"G", "G":"C", "T":"A", "N":"N"}
    reverse_seq = ""
    for base in seq:
        reverse_seq = complement[base] +reverse_seq
    return reverse_seq

def naiver_reverse_compl(pattern,text):
    """This function combines the naive and reverse compl, to find the occurances of a pattern (of a reverse. compl) to
    a text
    Input: The pattern for which to find reve.complement occurance in the text
    text (str) The text in which to search for reverse complement occurances
     Returns:
    list: A list containing the indices of all occurrences of the reverse complement pattern within the text
    """
    reversed_pattern = reverse_complement(pattern)
    print(reversed_pattern)
    matches_reverse = naive(reversed_pattern, text)

    return matches_reverse



def read_fastQ(filename):
    """This function reads seqss and qualiity score from a FASTQ
    Input: - filename(str)- The name of the fastq file to read.

    Returns: a tuple containing two lists:one containing the seqs, and one the qualities.
    """
    sequences = []
    qualities =[]
    with open(filename) as fh:
        while True:
            fh.readline()
            seq = fh.readline()
            fh.readline()
            qual = fh.readline()
            if len(seq) == 0:
                break
            sequences.append(seq)
            qualities.append(qual)
    return sequences,qualities



def read_reference_genome(filename):
    """Function to read a refenrece genome from a FASTA file
    Input: -filename(str): The name of the FASTA file containging the reference genome
    Returns: str: The concatenated genome sequence read from the FASTA file."""

    genome = ''
    with open(filename) as fh:
        for line in fh:
            if not line.startswith(">"):
                genome +=line.rstrip()
    return genome

def naive_22m(pattern, text):
    """THis function extends the naive  exact matching algorithm to allow fofr 2 mismatches
    Input: pattern (str): The pattern to search within the text.
    Text:: the text in which to search for occurances

    Returns:  A List containing the indices of all occurances of the pattern within the text."""

    occurances_2mm = []
    for i in range(len(text) -len(pattern) +1):
        count_mismatch = 0
        for j in range(len(pattern)):
            if not text[i+j] == pattern[j]:
                count_mismatch +=1
        if count_mismatch <3:
            occurances_2mm.append(i)
    return occurances_2mm

def phred_33toQ(qual):
    """This function calculates the quality of a single characters (phrd+33) and returns numerical quality
    Parameters:

    Input-qual THe phred encocded quality score character to conver to a numerical quality
    Returns"
    int: the numerical quality score corresponding to the input Phrd+33 encoded quality.
    """
    return ord(qual)-33

if __name__ == "__main__":

    #Read the reference genome
    reference_genome = read_reference_genome("lambda_virus.fa")


    #test how many times this AGGT exists in lamda
    occurances_test = naive("AGTCGA", reference_genome,)

    #Test for its reverse complement within lamda genome
    reverse_occurance =naiver_reverse_compl("AGTCGA",reference_genome)

    #with up to 22 mismatchs
    occurances_2mm = naive_22m("AGGAGGTT",reference_genome)
    print(occurances_2mm)


    #Read fastq seqs and qual
    seq_man, qual_human = read_fastQ("ERR037900_1.first1000.fastq")

    hist_per_pos = [0]*101

    #Calculate the average quality per position and create a histogram
    for qual in qual_human:
        for i in range(len(qual)):
            score = phred_33toQ(qual[i])
            hist_per_pos[i] +=score


    for i in range(len(hist_per_pos)):
        hist_per_pos[i] /=100


    #Plot the histogram for qual per position
    import matplotlib.pyplot as plt
    plt.bar(range(len(hist_per_pos)), hist_per_pos)
    plt.show()
    test = phred_33toQ("#")

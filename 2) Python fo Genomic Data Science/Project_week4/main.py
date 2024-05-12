"""Project for the Genomic data science/ Python for genomic data science"""

working_file = "dna2.fasta"



def parse_fasta(filename):

    """
    It takes as an input a fasta file.
    Returns a dictionary with key the identifier and as a value the sequence
    """

    with open(filename, "r") as file:
        contents = file.readlines()
    dic_contents = {}
    for line in contents:
        line = line.rstrip()
        if line.startswith(">"):
            id = line[17:31]
            dic_contents[id] = ""
        else:
            dic_contents[id] = dic_contents[id] +line

    return dic_contents


def count_ids(filename):
    """Takes a fasta file as an input
     Returns the numbers of ids in the fasta file"""
    with open(filename, "r") as file:
        contents = file.readlines()

    num_ids = 0
    for line in contents:
        line = line.rstrip()
        if line.startswith(">"):
            num_ids +=1

    return num_ids

def find_lengths(dic_contents):
    """It takes a dictionary as an input, with key the identifier and the value the sequence

    Returns a different dic with id as a key and length o  a sequence as a value.
    """
    dic_lengths = {}
    for id, seq in dic_contents.items():
        length_seq = len(seq.strip())
        dic_lengths[id] = length_seq
    return dic_lengths


def find_longest(dic_contents):
    """It takes a dictionary as an input, with key the identifier and the value the sequence

    Returns a dic with the a id as key and the longest seqence, if there ar more than one sequence with same
    max length, all going to be shown.
    """
    longest = 0
    dic_longest = {}
    long_seq = ""
    for id, seq in dic_contents.items():
        if len(seq)>longest:
            longest = len(seq)
            long_seq = id
        elif len(seq) == longest:
            dic_longest[id] = len(seq)
    dic_longest[long_seq] = longest
    return dic_longest

def find_shortest(dic_contents):
    """It takes a dictionary as an input, with key the identifier and the value the sequence

       Returns a dic with the a id as key and the shortest seqence, if there ar more than one sequence with same
       shortest length, all going to be shown.
       """

    shortest = 10**9
    dic_shortest = {}
    shortest_seq = ""
    for id, seq in dic_contents.items():
        if len(seq) < shortest:
            shortest = len(seq)
            shortest_seq = id
        elif len(seq) == shortest:
            dic_shortest[id] = len(seq)
    dic_shortest[shortest_seq] = shortest
    return dic_shortest


def reverse_complement(sequence):
    """It takes as an input a sequence and returns the reverse complement sequence
    Input: a sting of a sequence
    Returns: string of the sequence, reverse complement
    """
    complementary_dict = {"A":"T", "G":"C", "T":"A", "C":"G"}

    complement_seq = [complementary_dict[base] for base in sequence]

    return "".join(complement_seq)


def prepare_reading_frames(sequence):
    """It creates all possible reading frames
    Inputs: takes a DNA string
    Returns: a list with 6 possible reading frames in triples"""
    #only forward reaing frames
    reading_frames = [sequence[i:] for i in range(3)]

    #reverse reading frames
    #first find rverse complement
    reverse_sequence = reverse_complement(sequence)
    return reading_frames



def find_ORFs(selected_frame, list_reading_frames):
    """This function it finds all the ORFs in a sequence

    Input: It takes as input a list of reeading frames, only for forward strand
    Returns a list of ORFs for a the given sequence"""
    ORFS = []
    stop_codons = ["TAG", "TAA", "TGA"]

    #chosen reading frame
    if selected_frame =="all":
        selected_frame =None
    else:
        selected_frame = int(selected_frame)

    for frame in range(len(list_reading_frames)) if selected_frame is None else [selected_frame]:
        seq = list_reading_frames[frame]
        for startmatch in range(0,len(seq),3):
            if "ATG" == seq[startmatch:startmatch+3]:
                sequence = "ATG"
                for endmatch in range(startmatch+3,len(seq),3):
                    codon = seq[endmatch:endmatch+3]
                    sequence +=codon
                    if codon in stop_codons:
                        if len(sequence)>6:
                            if sequence not in ORFS:
                                ORFS.append(sequence)
                        break
    return ORFS


def find_longest_ORF(ORFS):
    """It calculates the longest potential reading frame
    Input:It takes a list of ORFS
    Returns the longest orf from this list of ORFs.
    """
    try:
        longest = max(ORFS,key=len)
        return longest, len(longest)
    except ValueError:
        pass
        #print("empty sequence")




def find_orfs3(list_reading_frames):
    ORFS = []
    start_codon = "ATG"
    stop_codons = ["TAA", "TAG", "TGA"]
    for seq in list_reading_frames:
        start_indices = [i for i in range(len(seq)) if seq.startswith("ATG", i)]
        for start in start_indices:
            sequence = "ATG"
            for i in range(start + 3, len(seq), 3):
                codon = seq[i:i + 3]
                sequence += codon
                if codon in stop_codons:
                    if len(sequence) > 6:  # Minimum ORF length check
                        if sequence not in ORFS:
                            ORFS.append(sequence)
                    break
    return ORFS




def print_orfs(selected_frame,dict_positions, dict_longest_orf_per_seq):
    """This function prints the data in the console, longest orf with positions of the longest orf

    Input: Takes 2 dictionaries, one dic (key:id, value: longest orf of the id)
        Second dic: [dict_positions], :key:id, value:position of orf)
    """
    id_longest_orf_within_seqs= max(dict_longest_orf_per_seq, key = dict_longest_orf_per_seq.get)
    length_longest_orf = dict_longest_orf_per_seq[id_longest_orf_within_seqs]
    position_orf = dict_positions[id_longest_orf_within_seqs]
    message = (f"The longest orf at frame {selected_frame} is {id_longest_orf_within_seqs} with length {length_longest_orf}"
          f"and starting position {position_orf}")
    return message



def find_repeats(filename):
    """This function it finds the most frequent repeat in all sequences of the file

    Input: fileanem and it asks the length of the repeats that want to examine for
    Ouput: The most frequence repeat with the number of occurances.
    """
    n = int(input("Please type the length of the repeat you want to examine: "))
    occur = {}
    file_cont = parse_fasta(filename)

    for id,seq in file_cont.items():
        for i in range(len(seq)-n+1):
            substring = seq[i:i+n]
            if substring in occur:
                occur[substring] +=1
            else:
                occur[substring] =1

    repeat_numbe_of_occur = max(occur,key=occur.get)
    number_of_occuran = occur[repeat_numbe_of_occur]
    message = f"The repeat that occured most frequently is {repeat_numbe_of_occur} with number of occurance " \
              f"{number_of_occuran}"

    return message
def main():
    dic_fasta = parse_fasta(working_file)
    dic_longest = find_longest(dic_fasta)
    dic_shortest = find_shortest(dic_fasta)

    number_of_ids = count_ids(working_file)
    print(f"The number of sequence in the file are {number_of_ids}")

    longest_seq_in_file = find_longest(dic_fasta)
    print(f"THe longest sequence in the file is {longest_seq_in_file}")

    shortest_seq_in_file = find_shortest(dic_fasta)
    print(f"THe longest sequence in the file is {shortest_seq_in_file}")

    """initiate a dic to append each longest ORF for each time of loop, key longest orf of sequence/ value the position that 
    orf stars"""
    dict_longest_orfs = {}
    dict_positions_orf = {}
    user_input_frame = input("Please type the frame (as a number) that you want to examine for ORFS "
                             "or press (all) if you want to check for all reading frames:")

    for id, seq in dic_fasta.items():
        reading_frames = prepare_reading_frames(seq)
        try:
            # select the reading frame that you want to look for
            ORFs1 = find_ORFs(selected_frame=user_input_frame, list_reading_frames=reading_frames)
            orfs_method_1 = find_longest_ORF(ORFs1)
            longest_seq_orf1 = orfs_method_1[0]
            len_longest_seq_method1 = orfs_method_1[1]
            dict_longest_orfs[id] = len_longest_seq_method1
            dict_positions_orf[id] = seq.find(longest_seq_orf1)
        except TypeError:
            continue

    print_message = print_orfs(dict_longest_orf_per_seq=dict_longest_orfs, dict_positions=dict_positions_orf,
                               selected_frame=1)
    print(print_message)
    repeats = find_repeats(working_file)

    # #second method find orf
    # try:
    #     orfs3 = find_orfs3(reading_frames)
    #     orfs_method_1 = find_longest_ORF(orfs3)
    #
    #     longest_seq_orf2 = orfs_method_1[0]
    #     len_longest_seq_method2 = orfs_method_1[1]
    #     print(f"Method 2 find orfs.{id}, {len(orfs3)}, longest ORF: {len_longest_seq_method2} at position {seq.find(longest_seq_orf2)}")
    # except TypeError:
    #     print("Method1. empty sequence ")
if __name__ == "__main__":
    main()







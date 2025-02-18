#!/usr/bin/env python 

# Motif Marker by Wesley Gomersall
#   This program takes fasta and motif.txt inputs and generates a png image
#   with the same name as the input fasta which visualizes motif sequences 
#   as they appear in the fasta-format sequences.

import argparse
import cairo
import re

def get_args():
    parser = argparse.ArgumentParser(description="DNA motifs visualized to-scale as they appear in sequences given in fasta file.")
    parser.add_argument("-f", "--fasta", help="File path to fasta file. Entries of this file must lower-case introns and upper-case exons.", required=True, type=str)
    parser.add_argument("-m", "--motif", help="File path to motif file. One motif per line.", required=True, type=str)
    return parser.parse_args()

IUPAC_REGEX = {'A': '[Aa]', 'a': '[Aa]',
               'C': '[Cc]', 'c': '[Cc]',
               'T': '[Tt]', 't': '[Tt]', 
               'G': '[Gg]', 'g': '[Gg]',
               'U': '[UuTt]', 'u': '[UuTt]',
               'W': '[AaTt]', 'w': '[AaTt]',
               'S': '[CcGg]', 's': '[CcGg]',
               'M': '[AaCc]', 'm': '[AaCc]',
               'K': '[GgTt]', 'k': '[GgTt]',
               'R': '[AaGg]', 'r': '[AaGg]',
               'Y': '[CcTt]', 'y': '[CcTt]',
               'B': '[CcGgTt]', 'b': '[CcGgTt]',
               'D': '[AaGgTt]', 'd': '[AaGgTt]',
               'H': '[AaCcTt]', 'h': '[AaCcTt]',
               'V': '[AaCcGg]', 'v': '[AaCcGg]',
               'N': '[AaCcTtUuGg]', 'n': '[AaCcTtUuGg]'}

class FastaRecord:
    '''Input header line (newline stripped) and entire sequence from each fasta entry.

    Attributes: 
        header (str):       Fasta header
        sequence (str):     Fasta sequence 
                            (intron lower-case, exon upper-case,
                            exactly 0 or 1 exons with length > 1)
        exonstart (int):    1-based index of first base in exon
        exonend (int):      1-based index of last base in exon
        length (int):       Length of fasta sequence

    Method(s): 
        show():             Prints fasta header and sequence.
    '''
    def __init__(self, head, seq):
        self.header: str = head # remove (reverse-complement) here and add as an attribute? 
        self.sequence: str = seq

        exonlocation: list = [dex for dex, char in enumerate(self.sequence) if char.isupper()]
        
        if exonlocation != []: # check for one and only one exon in self.sequence with len > 1
            for i in range(len(exonlocation) - 1): 
                assert exonlocation[i] + 1 == exonlocation[i+1] # single continuous exon
            self.exonstart: int = exonlocation[0] + 1 
            self.exonend: int = exonlocation[-1] + 1
        else: 
            self.exonstart: int = None
            self.exonend: int = None

        self.length: int = len(seq) 

    def show(self): 
        print(self.header) 
        print(self.sequence)

class Motif:
    '''Take in string for each motif given in input file. 

    Attributes: 
        name (str):         Same string as found in input
        regex (str):        Regular expression based on IUPAC codes
        length (int):       Length of motif sequence

    Method(s): 
        convert_regex():    Creates regex from name based on 
                            IUPAC codes.
        locate():           Returns a list of 1-based indices of
                            the first character of self.regex regular 
                            expression in a given sequence.
    '''
    def __init__(self, given_motif): 
        self.name: str = given_motif.strip()
        self.regex: str = self.convert_to_regex()
        self.length: int = len(self.name) 

    def convert_to_regex(self) -> str:
        reg: str = "" 
        for char in self.name: 
            reg = reg + IUPAC_REGEX[char]
        return reg

    def locate(self, record: FastaRecord) -> list:
        '''Returns a list of all occurances of self.regex in sequence.
        Locations in sequence are from 1-based index'''
        locations = []
        # use lookahead to get overlapping matches, loop through matches
        for match in re.finditer(rf"(?=({self.regex}))", record.sequence): 
            locations.append(match.start() + 1)
        return locations

def generate_motif_list(filepath: str) -> list:
    '''Create list of motifs (entries class Motif).'''
    motifs: list = []
    with open(filepath, 'r') as mot_file: 
        for mot in mot_file:
            motifs.append(Motif(mot))
    return motifs

def parse_fasta(filepath: str, list_of_motifs: list) -> list: 
    '''Loop through fasta file, for each record call function find_motifs()
    Return a list of lookup tables from find_motifs(). 
    Each element of the list corresponds to a fasta entry.'''
    
    coordinates: list = []

    with open(filepath, 'r') as fasta_in:
        first = True
        for line in fasta_in:
            if line.startswith('>') or not line:

                if not first: # create FastaRecord and generate lookup table
                    entry = FastaRecord(header, sequence)
                    temp_dict: dict = find_motifs(entry, list_of_motifs)
                    coordinates.append(temp_dict)

                header = line.strip()
                sequence = "" 
                first = False
            else: 
                sequence = sequence + line.strip()

        else: # create FastaRecord and generate lookup table
            entry = FastaRecord(header, sequence)
            temp_dict: dict = find_motifs(entry, list_of_motifs)
            coordinates.append(temp_dict)

    return coordinates

def find_motifs(record: FastaRecord, motifs: list):
    '''Do something with record and return something else
    Each lookup table in this list is of the form:
        {header: [exonstart, exonend, length],
         motif1: [occurance1, occurance2, ..., occuranceN],
         motif2: [occurance1, occurance2, ..., occuranceN],
         ...
         motifM: [occurance1, occurance2, ..., occuranceN]} 
    these lists do not necessarily have the same N'''

    coord_dict: dict = dict()

    coord_dict[record.header] = [record.exonstart, record.exonend, record.length]
    for motif in motifs: 
        coord_dict[motif.name] = motif.locate(record) 
    return coord_dict

def get_colors(num_motifs: int) -> list: # WIP
    '''Generate list of color hex codes based on how many motifs need to be plotted.'''
    pass

def draw(filename: str, coordinates: list, colors: list): # WIP
    '''Iterate through list of dictionaries. For each, normalize coordinates and cairo draws '''

    # first, determine longest fasta record (header: [2] has this information) 
    # use longest to normalize the lengths of all subsequent coordinates
    # maybe this needs to be its own function? 

    # draw for each record
    # how to deal with perfectly overlapping motifs? 

    # draw motif color key

    pass

if __name__ == "__main__":
    args = get_args()

    motif_list = generate_motif_list(args.motif)
    num_motifs = len(motif_list) 

    list_for_draw = parse_fasta(args.fasta, motif_list ) 

    num_records = len(list_for_draw)

    colors = get_colors(num_motifs)

    draw(filename, list_for_draw, colors) 

    ###########
    # testing #
    ###########

    print(num_records)

    for item in list_for_draw: 
        print(item) 

    for i, motif in enumerate(motif_list): 
        print(motif.name)
        print(motif.regex)
        
    tr = FastaRecord("test header", "GGGATCGATCGATCGGGCCCCGGGCGGGG")
    tr2 = FastaRecord("test header", "GGGGG")
    tm = Motif("GGG")

    print(f"locating {tm.name} (regex: {tm.regex}) in {tr.sequence}: {tm.locate(tr)}")
    print(f"locating {tm.name} (regex: {tm.regex}) in {tr2.sequence}: {tm.locate(tr2)}")

    print(f"length of sequence {tr.sequence}: {tr.length}")
    print(f"length of sequence {tr2.sequence}: {tr2.length}")
    print(f"length of motif {tm.name}: {tm.length}")

    tr3 = FastaRecord("test header", "atcgATCGccc") 
    print(f"{tr3.sequence}")
    print(f"sequence length: {tr3.length}")
    print(f"exon starts at: {tr3.exonstart}, exon ends at: {tr3.exonend}")

#!/usr/bin/env python 

# Motif Marker by Wesley Gomersall
#   This program takes fasta and motif.txt inputs and generates a png image
#   with the same name as the input fasta which visualizes motif sequences 
#   as they appear in the fasta-format sequences.

import argparse
import cairo
import os
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
    '''Create list of motifs (entries class Motif).
    Input(s): 
        filepath (str):     String storing file name of input motif file.

    Output(s): 
        list:               List of Motif objects.
    '''
    motifs: list = []
    with open(filepath, 'r') as mot_file: 
        for mot in mot_file:
            motifs.append(Motif(mot))
    return motifs

def parse_fasta(filepath: str, list_of_motifs: list) -> list: 
    '''Loop through fasta file, for each record call function find_motifs()
    Return a list of lookup tables from find_motifs(). 
    Each element of the list corresponds to a fasta entry.

    Input(s): 
        filepath (str):     String storing file name of input fasta file.

    Output(s): 
        list:               List of dictionaries output by find_motifs().
    '''

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

def find_motifs(rec: FastaRecord, motifs: list) -> dict:
    '''Parse fasta record, generate dictionary with info about fasta exon and length, 
    and all motif locations in fasta record.

    Input(s): 
        rec (FastaRecord):  Fasta record to parse for motif 
                            sequences.
        motifs (list):      List of Motifs from generate_motif_list().

    Output(s): 
        dict:               {header: [exonstart, exonend, length],
                             motif1: [occurance1, occurance2, ...],
                             motif2: [occurance1, occurance2, ...],
                             ...
                             motifM: [occurance1, occurance2, ...]} 
    '''

    coord_dict: dict = dict()

    coord_dict[rec.header] = [rec.exonstart, rec.exonend, rec.length]
    for motif in motifs: 
        coord_dict[motif.name] = motif.locate(rec) 
    return coord_dict

def get_colors(motifs: dict) -> dict: 
    '''Generate dictionary of color codes for pycairo based on how many motifs need to be plotted.

    Input(s): 
        motifs (dict):     Output from generate_motif_list(). 

    Output(s): 
        dict:               {motif1: [1, 0, 0, 0]
                             motif2: [0, 1, 0, 0]
                             ...
                             motifM: [w, x, y, z]}
    '''

    colors: dict = dict()
    r = 0; g = 0; b = 0

    for i, m in enumerate(motifs): 
        match i % 3:
            case 0: 
                r += 1
            case 1:
                g += 1
            case 2: 
                b += 1
        if r == g and r == b and r != 0: 
                r = 0; g = 0
        color_list = [r, g, b, 0]
        
        colors.update({m.name: color_list})
    return colors

def draw(filename: str, coordinates: list, motifs:list, colors: list, png: bool = True): 
    '''Iterate through list of dictionaries. For each, normalize coordinates and cairo draws 

    Input(s): 
        filename (str):     String storing output file base name.
        coordinates (list): List of dictionaries output from parse_fasta()/find_motifs().
                            Each dictionary contains necessary coordinates for 
                            plotting introns, exon, and motifs for each fasta record.
        motifs (list):      List of Motif objects obtained from generate_motif_list().
        colors (dict):      Color dictonary obtained from get_colors().
        png (bool):         True/False whether to print to png (True), else print pdf.

    Output(s): 
        null():             Creates .png (default) or .pdf image of plotted motifs 
                            for each fasta record in the same figure. 
                            Base name of file given by filename input.
    '''
    
    # get motif names and thus the number of motifs to plot
    motif_names: set = set() 
    max_motif_len: int = 0 
    for m in motifs: # motifs is list of Motif objects
        motif_names.add(m.name)
        if len(m.name) > max_motif_len:
            max_motif_len = len(m.name)
    nummotifs: int = len(motif_names) 

    # get number of records and the length of the longest record for normalizing plots
    numrecords: int = len(coordinates) 
    longest_record: int = 0 
    for d in coordinates: # coordinates is a list of dictionaries output from find_motifs/parse_fasta
        for key in d.keys():
            if key not in motif_names and d[key][2] > longest_record: 
                longest_record = d[key][2]
                    
    # size of margins and record plots
    record_width: int = 1000
    record_height: int = 100 
    plot_height: float = record_height * .75
    key_height: int = 25
    margin_hor_left: int = 25
    margin_hor_right: int = 25
    margin_ver: int = 25
    title_height: int = 50
    title_size: int = 25
    label_size: int = 10

    # size of entire plot
    panel_height: int = 2 * margin_ver + title_height + record_height * numrecords + key_height * nummotifs
    panel_width: int = margin_hor_right + margin_hor_left + record_width

    if png: 
        surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, panel_width, panel_height) # for png
    else: 
        surface = cairo.PDFSurface(filename + ".pdf", panel_width, panel_height) # for pdf
    context = cairo.Context(surface) 

    # white background
    context.set_source_rgba(1, 1, 1, 1)
    context.rectangle(0, 0, panel_width, panel_height)
    context.fill()

    context.select_font_face("Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
    
    # draw plot title
    context.set_source_rgba(0, 0, 0, 1)
    context.set_font_size(title_size) 
    context.move_to(margin_hor_left, margin_ver + 0.3 * title_height)
    context.show_text("Motif Plot") 

    # draw each record
    for i, record in enumerate(coordinates): 

        x0: int = margin_hor_left
        y0: int = margin_ver + title_height + ((i + 1) - 0.5) * record_height 

        for key in record.keys():
            if key not in motif_names:
                seq_name: str = key
                scale: float = 1 / longest_record # scale = length of record / longest record
                tot_len: int = record[key][2]
                exon_begin = record[key][0]
                exon_end = record[key][1]

        # draw sequence name
        seq_name_y: int = margin_ver + title_height + 5 + i * record_height
        context.set_source_rgba(0, 0, 0, 1) # color black
        context.set_font_size(label_size) 
        context.move_to(margin_hor_left, seq_name_y)
        context.show_text(seq_name.strip('>').split()[0]) 
        
        # draw introns (total length of sequence) 
        intron2_end: float = margin_hor_left + record_width * tot_len * scale 

        context.set_source_rgba(0, 0, 0, 1) # color black
        context.set_line_width(2) # default line width = 2
        context.move_to(x0, y0) # (x,y), (0,0) is the top left of the canvas, (width, height) is bottom right
        context.line_to(intron2_end, y0)
        context.stroke()

        # draw exon
        exon_begin_x = margin_hor_left + record_width * exon_begin * scale 
        exon_end_x = margin_hor_left + record_width * exon_end * scale 
        exon_begin_y = y0 - 0.5 * plot_height

        context.set_source_rgba(0, 0, 0, 1)
        context.rectangle(exon_begin_x, exon_begin_y, exon_end_x - exon_begin_x, plot_height) # (x0, y0, width, height)
        context.fill()
        
        # draw ticks 
        if tot_len > 200:
            mark: int = 100 # tickmark every N bases
            tick_height: int = 10
            tickmarks: list = []
            # for t in range(1, math.ceil(tot_len / mark)): 
            for t in range(1, tot_len // mark + 1): 
                tickmarks.append(margin_hor_left + record_width * (t * mark) * scale) # NOT CORRECT

            for t in tickmarks: 
                context.move_to(t, y0 - tick_height) # (x,y), (0,0) is the top left of the canvas, (width, height) is bottom right
                context.line_to(t, y0 + tick_height)
                context.stroke()

        # draw motif(s) staggered on y axis
        motif_begin_y = y0 - 0.5 * plot_height

        for m, motif in enumerate(motifs): 
            motif_begin_y_m = motif_begin_y + m * plot_height / nummotifs
            motif_end_y_m = motif_begin_y_m + plot_height / nummotifs

            if motif.name in record.keys(): 
                for loc in record[motif.name]: 
                    color = colors[motif.name]
                    # make line width based on length of motif
                    motif_width: float = 4 * len(motif.name) / max_motif_len
                    context.set_line_width(motif_width)
                    context.set_source_rgb(color[0], color[1], color[2]) 
                    context.move_to(margin_hor_left + record_width * loc * scale, motif_begin_y_m) 
                    context.line_to(margin_hor_left + record_width * loc * scale, motif_end_y_m) 
                    context.stroke()

    # draw motif-color key
    for m, motif in enumerate(motifs): 
        key_begin_y: int = margin_ver + title_height + record_height * numrecords + 0.5 * key_height
        
        # Motif color dashes
        context.set_line_width(2)
        color = colors[motif.name]
        context.set_source_rgb(color[0], color[1], color[2]) 
        context.move_to(margin_hor_left, key_begin_y + m * key_height) 
        context.line_to(margin_hor_left + 25, key_begin_y + m * key_height) 
        context.stroke()

        # Motif names (text)
        context.set_source_rgba(0, 0, 0, 1) 
        context.move_to(margin_hor_left + 30, key_begin_y + m * key_height) 
        context.set_font_size(label_size)
        context.show_text(motif.name) 

    if png: 
        surface.write_to_png(filename + ".png")
    else:
        surface.finish() 

if __name__ == "__main__":
    args = get_args()

    motif_list = generate_motif_list(args.motif)

    list_for_draw = parse_fasta(args.fasta, motif_list ) 

    colors = get_colors(motif_list)

    basename = os.path.splitext(os.path.basename(args.fasta))[0]

    draw(basename, list_for_draw, motif_list, colors, True) 

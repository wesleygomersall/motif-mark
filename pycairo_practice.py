#!/usr/bin/env python

# Pycairo practice plotting by Wesley Gomersall

# Working code to generate a line and a rectangle, 
# not at the origin, using pycairo
# ... and a little more ;-)

import cairo
import re

# Locations of exon and motif
loc = {"Fasta header": [40, 80, 180], # [exonstart, exonend, totlength]
       "Motif_1": [50, 55, 70, 90, 130], # [motifloc1, motifloc2, ...]
       "Motif_2": [10, 20, 90],
       "Motif_3": [66]}

colors = {"Motif_1": [1, 0, 0], 
          "Motif_2": [0, 1, 0],
          "Motif_3": [0, 0, 1]}

longest = loc["Fasta header"][2] # pixel length of longest seq, there is only one in this example

# Ciaro 

width = 2 * 25 + longest; height = 100
surface = cairo.PDFSurface("WKG_example.pdf", width, height)
context = cairo.Context(surface) 

# Draw a line (introns) 

# Draw the fasta header as a title here 

# For multiple fasta entries, this length will be normalized using length of the longest seq 
scale_factor = 1 # this will change based on relative length of sequence

context.set_source_rgba(0, 0, 0, 1) # color black
context.move_to(25, 50) # (x,y), (0,0) is the top left of the canvas, (width, height) is bottom right
context.line_to(25 + scale_factor * longest, 50)
context.stroke()

# Draw Rectange (exon) 
exon = loc["Fasta header"] 
ex_start = 25 + (exon[0] * scale_factor) 
ex_end = 25 + (exon[1] * scale_factor) 

context.set_source_rgba(0, 0, 0, 1)
context.rectangle(ex_start, 40, ex_end - ex_start, 20) # (x0, y0, width, height)
context.fill()

# Remove the exon key from dic for the next section
loc.pop("Fasta header")

# Draw a few motifs
for key, locations in loc.items():
    for place in locations: 
        color = colors[key] # need to make a unique color from each key, this is from color dic
        context.set_source_rgb(color[0], color[1], color[2]) 
        context.move_to(25 + place * scale_factor, 40) 
        context.line_to(25 + place * scale_factor, 60)
        context.stroke()

        # Draw an element of the diagram key
        # color square + motif name

surface.finish()

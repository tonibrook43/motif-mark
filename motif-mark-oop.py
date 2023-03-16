#!/usr/bin/env python

#########
# Import #
#########

import argparse
import re
import cairo
from Bio import SeqIO
from collections import defaultdict

#########
# args #
#########

def get_args():
    parser = argparse.ArgumentParser(description="A program to convert cDNAgene_sequences into proteingene_sequences by searching for codons and translating into amino acids")
    parser.add_argument("-f", "--file", help="designates absolute file path to fasta file", required=True)
    parser.add_argument("-o", "--outfile", help="designates absolute file path to the output file", required=True)
    parser.add_argument("-m", "--motif", help="designates absolute file path to fasta file", required=True)
    
    return parser.parse_args()
args=get_args()


# to run script: ./motif-mark-oop.py -f /Users/tonibrooks/bioinfo/OOP/oneline_fasta.fa -m /Users/tonibrooks/bioinfo/OOP/Fig_1_motifs.txt -o /Users/tonibrooks/bioinfo/OOP/motif_outfile.txt

##################
# Pycairo set-up #
##################

def pycairo_set():
    '''Setting the parameters for pycairo surface and context'''
    #width and height of the image
    width= 1000
    height = 1000

    #surface variable with the specified width and height
    surface = cairo.ImageSurface(cairo.FORMAT_RGB24, width, height) #cairo.FORMAT_RGB24 did not work but cairo.FORMAT_ARGB32 did so yay

    #this creates a new context; we then use the methods of the context to 'draw' on the surface i.e., drawing lines and shapes
    context = cairo.Context(surface)

    context.set_source_rgb(255,255,255)
    context.paint()
    
    return surface, context

surface, context = pycairo_set()

#################################
# CONSTANTS AND CLASS VARIABLES #
#################################

motifs_file = args.motif
fasta_file = args.file
MOTIF_COLOR_LIST = [(1.0, 0.0, 0,0), (0.0, 0.0, 1.0), (0.0, 1.0, 0.0), (0.5, 0.0, 0.5)] #RBG colors as a list for the four motifs
MOTIF_COLOR_DICT = {} #empty dictionary that will have motif string (line) as a key and color as a value

motif_start = 0
motif_stop = 0
motif_dict={}
line=""

exon_start = 0
exon_end =0

#############
# functions #
#############

def draw_legend(context):
    context.move_to(800, 200) #set the starting position of the text
    context.select_font_face("Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL) #set the font
    context.set_font_size(16) #set the font size
    context.set_source_rgba(0, 0, 0) #set the color of the text

    #this will write the text for each item in the legend
    context.set_source_rgba(0, 0, 0) #set the color to red
    context.show_text("Legend:")
    context.move_to(800, 220)
    context.set_source_rgba(1, 0, 0) #set the color to red
    context.show_text("Exons")
    context.move_to(800, 240)
    context.set_source_rgba(0, 0, 1) #set the color to blue
    context.show_text("Introns")
    context.move_to(800, 260)
    context.set_source_rgba(0, 1, 0) #set the color to green
    context.show_text("Motifs")


###########
# classes #
###########


class Gene:
    def __init__(self, context, gene_start, gene_end, y_start_end, gene_item, header):
        self.context = context
        self.gene_start = gene_start #this is how you save x_start
        self.gene_end = gene_end
        self.y_start_end = y_start_end
        self.gene_item = gene_item
        self.header = str(header)

    def __repr__(self):
        return f"Motif({self.gene_start}, {self.gene_end}, {self.y_start_end}, {self.gene_item}, {self.header})"

    def draw_gene(self):

        #Width of 2 for a fine line
        self.context.set_line_width(2)

        self.context.set_source_rgba(0, 0, 0)
        self.context.stroke()

        #how far away from the margin on the x-axis
        LEFT_MARGIN_X = 75

        #how far away from the margin on the y-axis
        GENE_HEADER_Y_MARGIN = 75

        #start position
        self.context.move_to(self.gene_start+LEFT_MARGIN_X , self.y_start_end)        #(x,y)
        
        # end position of line from start position
        self.context.line_to(self.gene_end+LEFT_MARGIN_X, self.y_start_end)
        self.context.stroke()


        self.context.set_font_size(15)
        self.context.select_font_face("Arial",
                        cairo.FONT_SLANT_NORMAL,
                        cairo.FONT_WEIGHT_NORMAL)
        
        self.context.move_to(LEFT_MARGIN_X, self.y_start_end-GENE_HEADER_Y_MARGIN)

        self.context.show_text(self.header)


class Exon:
    def __init__(self, context, exon_start: int, exon_end: int, exon_start_end, sequence):
        self.context = context
        self.exon_start = exon_start
        self.exon_end = exon_end
        self.exon_start_end = exon_start_end
        self.sequence = sequence

    
    def draw_exon(self, context):
    
        # Width of 15 for a thicker line
        self.context.set_line_width(15)
        self.context.set_source_rgba(0.0,0.0,0.0) #exons will be in a dif color than intron
        self.context.fill() #fill the exon with color


        #how far away from the margin on the x-axis
        LEFT_MARGIN_X = 0

        #start position
        self.context.move_to(self.exon_start+LEFT_MARGIN_X, self.exon_start_end)        #(x,y)
        
        # end position of line from start position
        self.context.line_to(self.exon_end+LEFT_MARGIN_X, self.exon_start_end)
        self.context.stroke()


class Motif:
    def __init__(self, context, motif_start, motif_stop, y_start_end, line) -> None:
        self.motif_start = motif_start
        self.motif_stop = motif_stop
        self.line =line #motif string
        self.context = context
        self.y_start_end = y_start_end
        self.color = MOTIF_COLOR_DICT[line]

    def __repr__(self):
        return f"Motif({self.motif_start}, {self.motif_stop}, {self.line}, {self.y_start_end})"
    
    def draw_motif(self, context):
        
        # Width of 7 for a fine line and colors based on motif
        self.context.set_line_width(7)
        self.context.set_source_rgba(*self.color)
        
        #how far away from the margin on the x-axis
        LEFT_MARGIN_X = 0

        #start position
        self.context.move_to(self.motif_start+LEFT_MARGIN_X , self.y_start_end)        #(x,y)
        
        # end position of line from start position
        self.context.line_to(self.motif_stop+LEFT_MARGIN_X, self.y_start_end)
        self.context.stroke()

class MotifIUPAC:
    def __init__(self, motifs_file):
        self.iupac = {"A":'[A]',"C":'[C]',"G":'[G]',"T":'[T]',"U":'[U]',"R":'[AG]',"Y":'[CT]',"S":'[GC]',"W":'[AT]',"K":'[GT]',"M":'[AC]',"B":'[CGT]',"D":'[AGT]',"H":'[ACT]',"V":'[ACG]',"N":'[A-Za-z]',
             "a":'[A]',"c":'[C]',"g":'[G]',"t":'[T]',"u":'[U]',"r":'[AG]',"y":'[CT]',"s":'[GC]',"w":'[AT]',"k":'[GT]',"m":'[AC]',"b":'[CGT]',"d":'[AGT]',"h":'[ACT]',"v":'[ACG]',"n":'[A-Za-z]'}
        self.motifs_file = motifs_file
        self.regexed_motifs = []
        self.character_counts = []
   
    def parse_motifs(self):
        # motif_dict={}
        with open(self.motifs_file, 'r') as motifs_file:
            character_count =0
            for i, line in enumerate(motifs_file): 
                motif = "" 
                line = line.strip('\n').upper()
                for character in line:  #parsing through characters in each motif
                    if character.capitalize() in self.iupac:
                        motif += self.iupac[character.capitalize()]   #if the character is in self.iupac, add it to the motif string
                        character_count +=1
                    else:
                        raise ValueError('Error: motif contains characters that cannot be translated :/') #raise an error if the motif cannot be translated or doesnt contain the necessary characters
                motif_dict[line] = motif
                MOTIF_COLOR_DICT[line] = MOTIF_COLOR_LIST[i]
            return motif_dict

#########
# main  #
#########

#GENE
y_start_end= 0
gene_item = 0

#MOTIF
#this is where we are building the motif into an object
found_motifs = MotifIUPAC(args.motif)
found_motifs.parse_motifs()
motif_match = []
motifs = []
motif_dict = found_motifs.parse_motifs()

#EXON
exon_start_end = 0
exon_position_list = []

with open(args.file, "r") as fa:
    for record in SeqIO.parse(fa, "fasta"):
        header = record.description
        sequence = str(record.seq)

        ## GENE START AND END + DRAW GENE OBJECTS ##

        y_start_end += 200 #how far apart I want to draw my genes on the figure
        gene_item += 1 #this is to keep track of number of gene records      
        gene_start = 1 #this should be where the first nucleotide position is
        gene_end = len(sequence) #this is the end of the record sequence
        gene = Gene(context, gene_start, gene_end, y_start_end, gene_item, header) #call back class and all the objects
        gene.draw_gene()

        ## EXON START AND END + DRAW EXON OBJECTS ##

        exon_start_end += 200
        for exon_match in re.finditer(r'[A-Z]+', sequence):
            exon_position = (exon_match.start(), exon_match.end())
            exon_position_list.append(exon_position)
            exon_start, exon_end = exon_position
            exon = Exon(context, exon_start, exon_end, exon_start_end, sequence)
            drawn_exon = exon.draw_exon(context)

        ## MOTIF START AND END + DRAW MOTIF OBJECTS ##

        for line, motif in motif_dict.items():
            matches = re.finditer(motif, str(sequence), re.IGNORECASE)
            for match in matches:
                motif_match.append(match)
                motif_start = match.start()
                motif_stop = match.end()
                motif = {'motif': line, 'start': motif_start, 'stop': motif_stop}
                motifs.append(motif)
            motif_draw = Motif(context, motif_start, motif_stop, y_start_end, line) #call back class and all the objects
            drawn_motif = motif_draw.draw_motif(context)


#### OUTPUT FINAL IMAGE TO PNG FILE ###
draw_legend(context)
context.save()
surface.write_to_png('Figure_1.png')
surface.finish()

# motif-mark

Developed by Toni Brooks at the University of Oregon

This is a Python program that searches for known motifs in a DNA sequence and marks the occurrences of these motifs in the sequence. It uses object-oriented programming (OOP) principles to implement the solution and uses pycairo to draw the gene, exons, and motifs.

Requirements:
  Python 3.x
  pycairo
  argparse
  
The program requires two inputs:
  A FASTA file containing DNA sequences (sequences up to 1000 bases)
  A motif file containing motifs (up to 10 bases each, one motif per line in a text file)
  

The program is capable of handling motifs with ambiguous nucleotides and multiple sequences (up to 10 in the data provided) and multiple motifs (up to 5 in the data provided). This porgram can also handle overlapping motifs and denote introns/exons, and all features (motifs, introns, exons) are to scale.

The program outputs:
  A Single, well-labeled figure in PNG format.

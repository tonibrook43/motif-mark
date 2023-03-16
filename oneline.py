#!/usr/bin/env python

#########
# Import #
#########

import argparse


#########
# args #
#########

def get_args():
    parser = argparse.ArgumentParser(description="A program to convert multi-line fasta into one line fasta")
    parser.add_argument("-f", "--file", help="designates absolute file path to fasta file", required=True)
    parser.add_argument("-o", "--outfile", help="designates absolute file path to the output file", required=False)    
    return parser.parse_args()
args=get_args()

# to run script: ./oneline.py -f /Users/tonibrooks/bioinfo/OOP/Figure_1.fasta


def oneline_fasta():
    '''Makes FASTA sequences on one line.'''
    # make dict with headers as keys and sequences as values
    sequence_dictionary = {}
    header_line = ""
    with open(args.file, 'r') as fasta:
        for line in fasta:
            line = line.strip('\n')
            # only get header lines
            if line[0] == '>':
                header_line = line
            else:
                if header_line not in sequence_dictionary:
                    sequence_dictionary[header_line] = line
                else:
                    sequence_dictionary[header_line] += line
    # write out to file
    oneline_fasta = open('oneline_fasta.fa', 'w')
    for keys,vals in sequence_dictionary.items():
        oneline_fasta.write(str(keys) + '\n' + str(vals) + '\n')

    oneline_fasta.close()

oneline_fasta()
  
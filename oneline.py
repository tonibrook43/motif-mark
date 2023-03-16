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
    return parser.parse_args()
args=get_args()

def oneline_fasta():
    '''Makes FASTA sequences on one line.'''
    sequence_dictionary = {}
    header_line = ""
    with open(args.file, 'r') as fasta:
        for line in fasta:
            line = line.strip('\n')
            if line[0] == '>':
                header_line = line
            else:
                if header_line not in sequence_dictionary:
                    sequence_dictionary[header_line] = line
                else:
                    sequence_dictionary[header_line] += line
    # write out to file
    oneline_fasta = open('oneline_fasta.fa', 'w')
    for keys, values in sequence_dictionary.items():
        oneline_fasta.write(str(keys) + '\n' + str(values) + '\n')

    oneline_fasta.close()

oneline_fasta()
  

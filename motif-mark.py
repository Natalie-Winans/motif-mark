#!/usr/bin/env python

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Motif-finding tool
# Author: Natalie Winans
# Last Updated:
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import argparse
import cairo
import math

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# USER INPUT
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# def get_args():
#     parser = argparse.ArgumentParser()
#     parser.add_argument('-f', '--fasta', help = "input fasta file containing genes of interest")
#     parser.add_argument('-m', '--motifs', help = "input file containing motifs of interest")

# args = get_args()
# fasta = args.fasta
# motifs = args.motifs

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FUNCTIONS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



base_dict = {
    'U' : 'T',
    'W' : ('A', 'T'), 
    'S' : ('C', 'G'), 
    'M' : ('A', 'C'), 
    'K' : ('G', 'T'), 
    'R' : ('A', 'G'), 
    'Y' : ('C', 'T'),
    'B' : ('C', 'G', 'T'),
    'D' : ('A', 'G', 'T'),
    'H' : ('A', 'C', 'T'),
    'V' : ('A', 'C', 'G'),
    'N' : ('A', 'C', 'G', 'T'),
    'Z' : 0
}
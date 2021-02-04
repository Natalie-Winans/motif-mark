#!/usr/bin/env python

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Motif-finding tool
# Author: Natalie Winans
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import argparse
from os import getgid
import cairo
import math
import re
from itertools import groupby

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

def degenerate_bases():
    """Creates dictionary in which keys are IUPAC degenerate base symbols and
    values are regex statements representing all possible values in a DNA or 
    RNA sequence"""
    base_dict = {
        'A' : '[Aa]',
        'C' : '[Cc]',
        'G' : '[Gg]',
        'T' : '[TtUu]',
        'U' : '[UuTt]',
        'W' : '[AaTtUu]', 
        'S' : '[CcGg]', 
        'M' : '[AaCc]', 
        'K' : '[GgTtUu]', 
        'R' : '[AaGg]', 
        'Y' : '[CcTtUu]',
        'B' : '[CcGgTtUu]',
        'D' : '[AaGgTtUu]',
        'H' : '[AaCcTtUu]',
        'V' : '[AaCcGg]',
        'N' : '[AaCcGgTtUu]',
        'Z' : '[]'
    }
    return base_dict



def convert_motifs(motif_file):
    """Takes text file of motifs, converts to regex statements based on IUPAC conventions, returns list of converted motifs ready to be searched with regex"""
    base_dict = degenerate_bases()
    re_motifs = []
    with open(motif_file, 'r') as fh:
        for motif in fh:
            new_motif = ""
            motif = motif.upper().strip()
            for base in motif:
                new_motif += base_dict[base] 
            re_motifs.append(new_motif)
        return re_motifs



def parse_fasta(fasta_file):
    """Takes fasta file and yields header and seq string objects,
    adapted from biopython"""
    header, seq = None, []
    for line in fasta_file:
        line = line.rstrip()
        if line.startswith(">"):
            if header: 
                yield(header, ''.join(seq))
            header, seq = line, []
        else:
            seq.append(line)
    if header: 
        yield(header, ''.join(seq))


def get_positions(fasta):
    """Takes fasta file, extracts positional information from header lines 
    and exons as well as exon length"""
    # with open(fasta) as fh:
    start_pos = []
    end_pos = []
    seq_length = []
    exon_start = []
    exon_end = []
    exon_length = []
    for header, seq in parse_fasta(fasta): 
        # get sequence positional information
        pos = re.search('\d+-\d+', header).group(0)
        start_pos.append(pos.split('-')[0])
        end_pos.append(pos.split('-')[1])
        seq_length.append(len(seq))

        # # get exon positional information
        exon = re.search('[A-Z]+', seq)
        exon_start.append(exon.span()[0])
        exon_end.append(exon.span()[1])
        #exon_seq = exon.group()
        exon_length.append(len(exon.group()))

    return start_pos, end_pos, seq_length, exon_start, exon_end, exon_length
    
# def find_motifs(motif_file):
#     re_motifs = convert_motifs(motif_file)
#     motif1 = re.search(re_motifs[1])

#     return re_motifs, motif1

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MAIN CODE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


motif_file = 'Fig_1_motifs.txt'
print(convert_motifs(motif_file))

fasta = 'Figure_1.fasta'

with open(fasta, 'r') as fasta:
    start_pos, end_pos, seq_length, \
        exon_start, exon_end, exon_length \
            = get_positions(fasta)
    print(seq_length, '\n', exon_length)



# set surface dimensions based on number of genes and max length
width = int(max(seq_length)) + 100
height = len(start_pos) * 150 + 150
surface = cairo.SVGSurface('plot.svg', width, height)
ctx = cairo.Context(surface)
ctx.rectangle(0, 0, width, height)
ctx.set_source_rgba(1, 1, 1, 1)
ctx.fill()

line_y = 150
#draw read maps with exons and motifs
for i in range(len(seq_length)):
    ctx.set_line_width(5)
    ctx.set_source_rgba(1, 0, 0, 1)
    ctx.move_to(50, line_y)
    ctx.line_to(50 + exon_start[i], line_y)
    ctx.stroke()

    ctx.set_line_width(20)
    ctx.set_source_rgba(0, 1, 0, 1)
    ctx.move_to(50 + exon_start[i], line_y)
    ctx.line_to(50 + exon_end[i], line_y)
    ctx.stroke()

    ctx.set_line_width(5)
    ctx.set_source_rgba(1, 0, 0, 1)
    ctx.move_to(50 + exon_end[i], line_y)
    ctx.line_to(50 + seq_length[i], line_y)
    ctx.stroke()

#annotate map with positions
# ctx.move_to(50, line_y + 20)
# ctx.set_source_rgb(0, 0, 0)
# ctx.set_font_size(12)
# ctx.show_text(start_pos)

# ctx.move_to(seq_length, line_y + 20)
# ctx.set_source_rgb(0, 0, 0)
# ctx.set_font_size(12)
# ctx.show_text(end_pos)

    line_y += 150

surface.write_to_png('plot.png')

















# with open(fasta) as fh:
#     for header, seq in parse_fasta(fh): 
#         # get sequence positional information
#         pos = re.search('\d+-\d+', header).group(0)
#         start_pos = pos.split('-')[0]
#         end_pos = pos.split('-')[1]
#         seq_length = len(seq)

#         # get exon positional information
#         exon = re.search('[A-Z]+', seq)
#         exon_start = exon.span()[0]
#         exon_end = exon.span()[1]
#         exon_seq = exon.group()
#         exon_length = len(exon.group())



        
#         #draw read map with exon and motifs
#         ctx.set_line_width(5)
#         ctx.set_source_rgba(1, 0, 0, 1)
#         ctx.move_to(50, line_y)
#         ctx.line_to(50 + exon_start, line_y)
#         ctx.stroke()

#         ctx.set_line_width(20)
#         ctx.set_source_rgba(0, 1, 0, 1)
#         ctx.move_to(50 + exon_start, line_y)
#         ctx.line_to(50 + exon_end, line_y)
#         ctx.stroke()

#         ctx.set_line_width(5)
#         ctx.set_source_rgba(1, 0, 0, 1)
#         ctx.move_to(50 + exon_end, line_y)
#         ctx.line_to(50 + seq_length, line_y)
#         ctx.stroke()

#         #annotate map with positions
#         ctx.move_to(50, line_y + 20)
#         ctx.set_source_rgb(0, 0, 0)
#         ctx.set_font_size(12)
#         ctx.show_text(start_pos)

#         ctx.move_to(seq_length, line_y + 20)
#         ctx.set_source_rgb(0, 0, 0)
#         ctx.set_font_size(12)
#         ctx.show_text(end_pos)

#         line_y += 150

#     surface.write_to_png('plot.png')

















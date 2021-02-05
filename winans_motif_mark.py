#!/usr/bin/env python

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Motif-finding tool
# Author: Natalie Winans
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import argparse
from os import getgid, remove
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
    re_motifs = {}
    with open(motif_file, 'r') as fh:
        for motif in fh:
            new_motif = ""
            motif = motif.upper().strip()
            for base in motif:
                new_motif += base_dict[base] 
            re_motifs[motif] = new_motif
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
    """Takes fasta file, extracts headers, and gene and exon positional information"""
    # with open(fasta) as fh:
    headers = []
    seqs = []
    start_pos = []
    end_pos = []
    seq_length = []
    exon_start = []
    exon_end = []
    exon_length = []
    for header, seq in parse_fasta(fasta): 
        seqs.append(seq)
        headers.append(header)
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

    return headers, seqs, start_pos, end_pos, seq_length, exon_start, exon_end, exon_length
    
# def find_motifs(motif_file):


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MAIN CODE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


motif_file = 'Fig_1_motifs.txt'
fasta = 'Figure_1.fasta'

with open(fasta, 'r') as fasta:
    headers, seqs, start_pos, end_pos, seq_length, \
        exon_start, exon_end, exon_length \
            = get_positions(fasta)
    #print(headers, '\n', seq_length, '\n', exon_length)
    #print(seqs)

re_motifs = convert_motifs(motif_file)
all_motif_spans = {}
for i, seq in enumerate(seqs): 
    seq_motif_spans = {}
    for j, motif in enumerate(re_motifs):
        motif_match= re.finditer(re_motifs[motif], seq)
        #motifs0 = [m.group(0) for m in motif_j]
        seq_motif_spans[motif] = [m.span() for m in motif_match]
        #print(motif)
    
    all_motif_spans["seq{0}".format(i)] = seq_motif_spans
    # print(seq_motif_spans)

    
# set surface dimensions based on number of genes and max length
width = int(max(seq_length)) + 250
height = len(start_pos) * 150 + 150
surface = cairo.SVGSurface('plot.svg', width, height)
ctx = cairo.Context(surface)
ctx.rectangle(0, 0, width, height)
ctx.set_source_rgba(1, 1, 1, 1)
ctx.fill()

colors = [[.9, 0, 0, .8], \
    [1, .8, 0, .8], \
    [0, .8, .2, .8], \
    [0, .2, 1, .8], \
    [.6, .2, .9, .8]]

line_y = 150
#draw read maps with exons and motifs
for i, dict in enumerate(all_motif_spans):

    #label with gene name
    ctx.move_to(50, line_y - 35)
    ctx.set_source_rgb(0, 0, 0)
    ctx.set_font_size(16)
    ctx.show_text(headers[i])

    ctx.set_line_width(5)
    ctx.set_source_rgb(0, 0, 0)
    ctx.move_to(50, line_y)
    ctx.line_to(50 + exon_start[i], line_y)
    ctx.stroke()

    #draw exon
    ctx.set_line_width(20)
    ctx.set_source_rgb(.4, .4, .4)
    ctx.move_to(50 + exon_start[i], line_y)
    ctx.line_to(50 + exon_end[i], line_y)
    ctx.stroke()

    ctx.set_line_width(5)
    ctx.set_source_rgb(0, 0, 0)
    ctx.move_to(50 + exon_end[i], line_y)
    ctx.line_to(50 + seq_length[i], line_y)
    ctx.stroke()

    #draw motifs
    for j, motif in enumerate(all_motif_spans[dict]):
        spans = all_motif_spans[dict][motif]
        ctx.set_source_rgba(*colors[j])
        for span in spans:
            ctx.move_to(50 + span[0], line_y)
            ctx.set_line_width(30)
            ctx.line_to(50 + span[1], line_y)
            ctx.stroke()

    line_y += 150


#make legend
ctx.move_to(max(seq_length) + 50 ,100)
ctx.set_source_rgb(0, 0, 0)
ctx.set_font_size(18)
ctx.show_text("Motifs")

motifs = list(re_motifs.keys())
for i in range(len(seqs)):
    ctx.move_to(max(seq_length) + 50, 100 + 30 * (i + 1))
    ctx.set_source_rgba(*colors[i])
    ctx.set_line_width(30)
    ctx.line_to(max(seq_length) + 60, 100 + 30 * (i + 1))
    ctx.stroke()
    
    ctx.move_to(max(seq_length) + 80, 100 + 30 * (i + 1.25))
    ctx.set_source_rgb(0, 0, 0)
    ctx.set_font_size(16)
    ctx.show_text(motifs[i])






#annotate map with positions
# ctx.move_to(50, line_y + 20)
# ctx.set_source_rgb(0, 0, 0)
# ctx.set_font_size(12)
# ctx.show_text(start_pos)

# ctx.move_to(seq_length, line_y + 20)
# ctx.set_source_rgb(0, 0, 0)
# ctx.set_font_size(12)
# ctx.show_text(end_pos)



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

















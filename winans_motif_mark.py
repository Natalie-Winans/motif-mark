#!/usr/bin/env python

# Motif Mark: Tool for visualizing protein-binding motifs
# Author: Natalie Winans            


import argparse
import os
import re
import cairo


def get_args():
    parser = argparse.ArgumentParser(
        prog="MotifMark",
        description="Tool for visualization of motifs on sequences",
        add_help=True)
    parser.add_argument(
        '-f', '--fasta', help="input fasta file containing sequences of interest")
    parser.add_argument(
        '-m', '--motifs', help="input file containing motifs of interest")
    return parser.parse_args()


def extract_fasta_name(fasta_file):
    """Extract input filename without extension for figure title and output files."""
    return os.path.basename(fasta_file).split('.')[0]


def degenerate_bases():
    """Create dictionary in which keys are IUPAC degenerate base symbols and
    values are regex statements representing all possible values in a DNA or 
    RNA sequence."""
    base_dict = {
        'A': '[Aa]',
        'C': '[Cc]',
        'G': '[Gg]',
        'T': '[TtUu]',
        'U': '[UuTt]',
        'W': '[AaTtUu]',
        'S': '[CcGg]',
        'M': '[AaCc]',
        'K': '[GgTtUu]',
        'R': '[AaGg]',
        'Y': '[CcTtUu]',
        'B': '[CcGgTtUu]',
        'D': '[AaGgTtUu]',
        'H': '[AaCcTtUu]',
        'V': '[AaCcGg]',
        'N': '[AaCcGgTtUu]',
        'Z': '[]'
    }
    return base_dict


def convert_motifs(motif_file):
    """Take text file of motifs, convert to regex statements based on IUPAC conventions, 
    return list of converted motifs ready to be searched with regex."""
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
    """Take fasta file and yield header and seq string objects.
    (adapted from biopython)"""
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
    """Take fasta file, extract and return headers, 
    sequence and exon positional information"""
    headers = []
    seqs = []
    seq_start = []
    seq_length = []
    exon_start = []
    exon_end = []
    for header, seq in parse_fasta(fasta):
        seqs.append(seq)
        headers.append(header)
        # get sequence positional information
        pos = re.search('\d+-\d+', header).group(0)
        seq_start.append(pos.split('-')[0])
        seq_length.append(len(seq))
        # get exon positional information
        exon = re.search('[A-Z]+', seq)
        exon_start.append(exon.span()[0])
        exon_end.append(exon.span()[1])
    return headers, seqs, seq_start, \
        seq_length, exon_start, exon_end


def motif_spans(motif_file, seqs):
    """
    Take motif file and list of gene/pre-mRNA sequences,
    return dictionary of dictonaries of lists of tuples.

    Outer dictionary:
        keys = gene/pre-mRNA sequence indexes
        values = inner dictionaries:
            keys = motifs
            values = lists of tuples representing span of 
                     motif in the sequence
    """
    re_motifs = convert_motifs(motif_file)
    all_motif_spans = {}
    for i, seq in enumerate(seqs):
        seq_motif_spans = {}
        for j, motif in enumerate(re_motifs):
            motif_match = re.finditer(re_motifs[motif], seq)
            seq_motif_spans[motif] = [m.span() for m in motif_match]

        all_motif_spans["seq{0}".format(i)] = seq_motif_spans

    return all_motif_spans


def palette(num_motifs):
    """Take number of motifs (up to 6) and return 
    colorblind-friendly color palette."""
    if num_motifs == 1:
        pal = [.13, .53, .2, .8]
    if num_motifs == 2:
        pal = [[.13, .53, .2, .8],
               [.67, .2, .47, .8]]
    if num_motifs == 3:
        pal = [[.27, .47, .67, .8],
               [.13, .53, .2, .8],
               [.67, .2, .47, .8]]
    if num_motifs == 4:
        pal = [[.4, .8, .93, .8],
               [.13, .53, .2, .8],
               [.8, .73, .27, .8],
               [.67, .2, .47, .8]]
    if num_motifs == 5:
        pal = [[.27, .47, .67, .8],
               [.4, .8, .93, .8],
               [.13, .53, .2, .8],
               [.8, .73, .27, .8],
               [.67, .2, .47, .8]]
    if num_motifs == 6:
        pal = [[.27, .47, .67, .8],
               [.4, .8, .93, .8],
               [.13, .53, .2, .8],
               [.8, .73, .27, .8],
               [.93, .4, .47, .8],
               [.67, .2, .47, .8]]
    return pal


def main():
    args = get_args()
    fasta = args.fasta
    motif_file = args.motifs

    fasta_name = extract_fasta_name(fasta)

    with open(fasta, 'r') as fasta:
        headers, seqs, seq_start, \
            seq_length, exon_start, exon_end \
            = get_positions(fasta)

    all_motif_spans = motif_spans(motif_file, seqs)

    # set surface dimensions based on number of sequences and max length
    width = int(max(seq_length)) + 250
    height = len(seq_start) * 150 + 100
    surface = cairo.SVGSurface('%s.svg' % fasta_name, width, height)
    ctx = cairo.Context(surface)
    ctx.rectangle(0, 0, width, height)
    ctx.set_source_rgba(1, 1, 1, 1)
    ctx.fill()

    # set motif colors
    pal = palette(len(all_motif_spans))

    # write figure title from fasta name
    ctx.move_to(50, 60)
    ctx.set_source_rgb(0, 0, 0)
    ctx.set_font_size(24)
    ctx.show_text("FASTA file: %s" % fasta_name)

    # draw sequence maps with exons and motifs
    line_y = 150
    for i, dict in enumerate(all_motif_spans):

        # label with sequence name
        ctx.move_to(50, line_y - 35)
        ctx.set_source_rgb(0, 0, 0)
        ctx.set_font_size(16)
        ctx.show_text(headers[i].split(' ')[0])

        # draw 5' intron
        ctx.set_line_width(5)
        ctx.set_source_rgb(0, 0, 0)
        ctx.move_to(50, line_y)
        ctx.line_to(50 + exon_start[i], line_y)
        ctx.stroke()

        # draw exon
        ctx.set_line_width(35)
        ctx.set_source_rgb(.7, .7, .7)
        ctx.move_to(50 + exon_start[i], line_y)
        ctx.line_to(50 + exon_end[i], line_y)
        ctx.stroke()

        # draw 3' intron
        ctx.set_line_width(5)
        ctx.set_source_rgb(0, 0, 0)
        ctx.move_to(50 + exon_end[i], line_y)
        ctx.line_to(50 + seq_length[i], line_y)
        ctx.stroke()

        # label 5' and 3' ends
        ctx.move_to(35, line_y + 5)
        ctx.set_source_rgb(0, 0, 0)
        ctx.set_font_size(12)
        ctx.show_text("5'")
        ctx.move_to(55 + seq_length[i], line_y + 5)
        ctx.set_source_rgb(0, 0, 0)
        ctx.set_font_size(12)
        ctx.show_text("3'")


        # draw motifs
        for j, motif in enumerate(all_motif_spans[dict]):
            spans = all_motif_spans[dict][motif]
            ctx.set_source_rgba(*pal[j])
            for span in spans:
                ctx.move_to(50 + span[0], line_y)
                ctx.set_line_width(35)
                ctx.line_to(50 + span[1], line_y)
                ctx.stroke()

        line_y += 150

    # make legend
    leg_y = 115
    ctx.move_to(max(seq_length) + 50, leg_y)
    ctx.set_source_rgb(.7, .7, .7)
    ctx.set_line_width(30)
    ctx.line_to(max(seq_length) + 80, leg_y)
    ctx.stroke()
    ctx.move_to(max(seq_length) + 90, leg_y + 7)
    ctx.set_source_rgb(0, 0, 0)
    ctx.set_font_size(16)
    ctx.show_text("Exon")

    motifs = list(convert_motifs(motif_file).keys())
    for i in range(len(motifs)):
        ctx.move_to(max(seq_length) + 50, leg_y + 5 + 30*(i+1) + i*5)
        ctx.set_source_rgba(*pal[i])
        ctx.set_line_width(30)
        ctx.line_to(max(seq_length) + 80, leg_y + 5 + 30*(i+1) + i*5)
        ctx.stroke()

        ctx.move_to(max(seq_length) + 90, leg_y + 5 + 30*(i+1.25) + i*5)
        ctx.set_source_rgb(0, 0, 0)
        ctx.set_font_size(16)
        ctx.show_text(motifs[i])

    surface.write_to_png('%s.png' % fasta_name)


if __name__ == "__main__":
    main()

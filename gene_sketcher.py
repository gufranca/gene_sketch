#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Author: Gustavo Starvaggi Franca
Program name: gene_sketcher.py
Date: 2014-05-01
Last date modified: 2014-05-01
License: GPL

1. What it does:
   This program uses Genome Tools library to draw genes from GFF files. It also
   adds feature tracks in the sketch if the file is passed as an option.

2. Input:
   GFF3 file, STYLE file according to Genome Tools specifications, OUTPUT file,
   and FEATURES file (to be drawn as custom tracks). The last one is optional.

3. Output:
   A diagram of a gene and feature tracks.

4. Usage:
   python gene_sketcher.py <GFF_FILE> <STYLE_FILE> <OUTPUT_FILE> [<TRACKS_FILE>]
   or try python gene_sketcher.py --help

"""

import argparse
from gt.core import *
from gt.extended import *
from gt.annotationsketch import *
from gt.annotationsketch.custom_track import CustomTrack
from gt.core.gtrange import Range


class CustomTrackInsertions(CustomTrack):
    def __init__(self, sidelength, data):
        super(CustomTrackInsertions, self).__init__()
        self.sidelength = sidelength
        self.data = data
    
    def get_height(self):
        return 20
    
    def get_title(self):
        return "Insertion site"
    
    def render(self, graphics, ypos, rng, style, error):
        height = (self.sidelength*math.sqrt(3))/2
        margins = graphics.get_xmargins()
        t_color = Color(0.27, 0.741, 0.835, 1.0)
        for pos, desc in self.data.iteritems():
            drawpos = margins + (float(pos)-rng.start)/(rng.end-rng.start+1) \
                      * (graphics.get_image_width()-2*margins)
            graphics.draw_line(drawpos-self.sidelength/2, ypos + height,
                             drawpos, ypos,                                     
                             t_color, 1)
            graphics.draw_line(drawpos, ypos,                                   
                             drawpos+self.sidelength/2, ypos + height,          
                             t_color, 1)
            graphics.draw_line(drawpos-self.sidelength/2, ypos + height,        
                             drawpos+self.sidelength/2, ypos + height,          
                             t_color, 1)
            graphics.draw_text_centered(drawpos, ypos + height + 13, str(desc))
        return 0


def draw_tracks(TRACK_FILE, GFF_FILE, diagram):
    """
    Draw custom tracks if the file is passed as an argument.
    File must be in following format: Gene\ttrack_name\tposition
    The second argument is the diagram object, where the tracks will be
    inserted.

    """

    from track_check import read_track_file, overlap_recursive

    # get gene name from gff file
    gff = open(GFF_FILE, 'r')
    lines = gff.readlines()
    gene_name = lines[2].split(";")[1].split('=')[1].strip()

    track_sorted = read_track_file(TRACK_FILE)
    if gene_name in track_sorted:
        gene_track_list = track_sorted[gene_name]
        
        # check if there is more than one track for the gene, then deal with
        # overlapping tracks. If not, just copy gene_track_list.
        if len(gene_track_list) > 1:
            tracks = overlap_recursive(gene_track_list)
        else:
            tracks = [gene_track_list]
        
        # this is ugly, but we have to reorder the elements in tuple,
        # position must be the first and the track name is the second
        tracks_reorder = []
        for track in tracks:
            hold = []
            for t in track:
                hold.append((t[1], t[0]))
            tracks_reorder.append(hold)
    
        # add custom tracks for each track (each on each line)
        for track in tracks_reorder:
            # Ex: ctt = CustomTrackInsertions(15, {126243176:"foo",
            #                                      126255180:"bar"})
            ctt = CustomTrackInsertions(15, dict(track))
            diagram.add_custom_track(ctt)
    return diagram


def main():
    """
    Get arguments, build objects and call functions to draw diagrams.

    """

    parser = argparse.ArgumentParser(description="""Draw gene sketches and 
    custom tracks in the diagram""")

    parser.add_argument('gff_file', help="GFF3 file for a gene.")
    parser.add_argument('style_file', help="Genome Tools style file.")
    parser.add_argument('output_file', help="Output file")
    parser.add_argument('-t', dest='tracks_file', help="""Features file
                         containing tracks to be drawn.""")
    args = parser.parse_args()

    #--------------------------------#
    # Draw genes using Genome tools! #
    #--------------------------------#
    
    style = Style()
    style.load_file(args.style_file)

    # read gff file
    in_stream = GFF3InStream(args.gff_file)
    # add introns between exons:
    add_introns_stream = AddIntronsStream(in_stream)
    # create feature index (one per gene, so just one in our case.)
    feature_index = FeatureIndexMemory()
    # populate the feature index with our data stream
    feature_stream = FeatureStream(add_introns_stream, feature_index)

    # XXX this needs to be run before querying the feature index!
    gene_node = feature_stream.next_tree()

    seqid = feature_index.get_first_seqid()
    feature_range = feature_index.get_range_for_seqid(seqid)

    final_diagram = Diagram.from_index(feature_index, seqid,
                                       feature_range, style)

    if args.tracks_file:
        final_diagram = draw_tracks(args.tracks_file, args.gff_file,
                                    final_diagram)

    layout = Layout(final_diagram, 680, style)
    height = layout.get_height()
    canvas = CanvasCairoFile(style, 680, height)
    layout.sketch(canvas)
    
    # TODO! Choose the output format (png, pdf, svg)
    canvas.to_file(args.output_file)


if __name__ == "__main__":
    main()
    

#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Author: Gustavo Starvaggi Franca
Program name: track_check.py
Date: 2014-04-26
Last date modified: 2014-04-26
License: GPL

1. What it does:
   This module was designed to get a set of features for each gene, read 
   these features and check if they overlap each other. This is done
   "recursively" to get a final output of a list of lists of non overlapping
   elements. This module was primarily used to draw miRNA features in host
   genes and have to be imported in your script of genome tools.

2. Input:
   A tab separated file containing Gene Feature Position.

3. Output:
   Depends on the functions used. The final output is a list of lists with non-
   overlapping elements. This list must be handled in your script using genome
   tools library.

4. Usage:
   As a module: import track_check *
   As a script: python track_check.py <features_list.txt> > <output>


"""

import sys
from operator import itemgetter

def read_track_file(track_file):
    """
    Read a tab delimited file containing track information.
    Ex: GENE track_name position
    Returns -> A dictionary with position sorted tracks
    track_sorted["GENE"] = [('track1', pos1), ('track2', pos2)]

    """

    track_info = {}
    track_sorted = {}
    with open(track_file, 'r') as f:
        for line in f:
            data = line.strip().split("\t")
            gene, track, track_pos = data[0], data[1], int(data[2])
            if gene in track_info:
                track_info[gene].append((track, track_pos))
            else:
                track_info[gene] = [(track, track_pos)]
        for gene, tracks in track_info.items():
            track_sorted[gene] = sorted(tracks, key=itemgetter(1))        
        return track_sorted


def check_overlap(gene_track_list, minimum=5000):
    """
    Check if the track elements that overlap with previous one.
    Returns -> Two lists of overlaped elements and non-overlaped elements.

    """
    
    assert len(gene_track_list) > 1, """Track list must contain more than one
    element. Use a conditional in your track list before calling it."""
    
    prev_track = -100000
    overlap = []
    not_overlap = []
    for track in gene_track_list:
        # append to overlap list if the actual overlaps to the previous
        if track[1] - prev_track <= minimum:
            overlap.append((track[0], track[1]))
        # if not overlap, append to not overlap list
        else:
            not_overlap.append((track[0], track[1]))
        prev_track = track[1]    
    return overlap, not_overlap

            
def overlap_recursive(gene_track_list):
    """
    Calls check_overlap() for a list of tracks containing N elements until 
    there are no overlapping elements. This creates subtracks containing only
    non-overlapping tracks.
    Returns -> A final list of lists containing non-overlapping subtracks.
               [[(track1, pos1), (track2, pos2)], [(track3, pos3)...]...]
    """
    
    assert len(gene_track_list) > 1, """Track list must contain more than one
    element. Use a conditional in your tracks list before calling it."""
    
    overlap, not_overlap = check_overlap(gene_track_list)
    final_tracks = [not_overlap]
    # Call check_overlap() again only if overlap is not empty
    if overlap:
        while len(overlap) > 1:
            overlap, not_overlap = check_overlap(overlap)
            final_tracks.append(not_overlap)
        else:
            final_tracks.append(overlap)
    return final_tracks


if __name__ == "__main__":
    
    # just to test the functions
    f = sys.argv[1]
    a = read_track_file(f)
    gene_track_list = a["CHST13"]

    final = overlap_recursive(gene_track_list)
    print final



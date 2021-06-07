#!/usr/bin/env python3

"""Assembly based on Ovelap Layout Consensus"""

from Bio import SeqIO
from numpy import zeros
import itertools
import argparse
import time
import os

parser = argparse.ArgumentParser()
parser.add_argument('input_reads', type=str, help="Input reads to assemble (fasta format)")
parser.add_argument('output_contigs', type=str, help="Output contigs (fasta format)")
parser.add_argument('-cov','--coverage', type=int, required=False, help="Sequencing coverage", default=5 )
parser.add_argument('-minlap', '--minmial_overlap', type=int, required=False, help="Minimal overlap for the sequences", default=20)
args = parser.parse_args()


def fasta_format(filename):
	# check if name for fastafile is correct
	if not filename.endswith('.fasta'):
		filename = filename + '.fasta'
	return filename


def readin_fasta(reads_file):
	"""Get read sequences from fastafile"""
	read_indexed = dict()
	for ind, record in enumerate(SeqIO.parse(reads_file, 'fasta')):
		read_indexed[ind] = str(record.seq)
	return read_indexed


def overlap(a, b, k_min=3):
    """ Return length of longest suffix of 'a' matching
        a prefix of 'b' that is at least 'min_length'
        characters long.  If no such overlap exists,
        return 0. """
    start = 0  # start all the way at the left
    while True:
        start = a.find(b[:k_min], start)  # look for b's suffx in a
        if start == -1:  # no more occurrences to right
            return 0
        # found occurrence; check for full suffix/prefix match
        if b.startswith(a[start:]):
        	if start == 0: #if the sequences are identical they are not following reads!
        		return 0
        	return len(a)-start
        start += 1  # move just past previous match


def overlapScore(xc, yc):
    """Cost function assigning 0 to matching and 4 to to unmatching, 2 for gap in algnment"""
    minc, maxc = min(xc, yc), max(xc, yc)
    if xc == yc: return 0 # match
    if minc == 'A' and maxc == 'G': return 2 # transition
    elif minc == 'C' and maxc == 'T': return 2 # transition
    return 4


def approximate_overlap(x, y, s, k_min, minval=8):
	"""find global overlap - below the minval - threshold, for x,y -sequences, 
	s - score function, minlen - minimal suffix lengthm"""
	D = zeros((len(x)+1, len(y)+1), dtype=int)
	for j in range(1, len(y)+1):
		D[0, j] = 1000
	for i in range(1, len(x)+1):
	    D[i, 0] = 0
	for i in range(1, len(x)+1):
	    for j in range(1, len(y)+1):
	        D[i, j] = min(D[i-1, j-1] + s(x[i-1], y[j-1]), # diagonal
	                      D[i-1, j  ] + s(x[i-1], '-'),    # vertical
	                      D[i  , j-1] + s('-',    y[j-1])) # horizontal
	bestlap = 0
	for lap, scr in enumerate(reversed(D[-1][k_min:])):
		if scr < minval:
			minval = scr
			bestlap = len(D[0]) - (lap + 1)
			return bestlap
	return bestlap


def overlap_graph(reads, k_min, k_enough=30, min_approx=12):
    """ Find overlaps for sequence a above the set threshold"""
    overlaps = dict()
    start = time.time()
    print('Building graph start time:', start)
    for a, b in itertools.permutations(reads.items(), 2):
    	if a[1] == b[1]:
    		continue
    	olen = overlap(a[1], b[1], k_min)
    	if olen != 0:
    		if a[0] not in overlaps:
    			overlaps[a[0]] = dict()
    		overlaps[a[0]][b[0]] = olen
    print(f'Finished basic overlap finding in {(time.time() - start)/60} minutes.')
    for a in reads.items():
    	# if we find staight forward very good overlap we can continue to next one
    	if a[0] in overlaps.keys():
    		i, olap = find_longest_overlap(a[0], overlaps)
    		if olap > k_enough:
    			continue
    	for  b in reads.items():
    		if a[1] == b[1]: continue
    		olen = approximate_overlap(a[1], b[1], s=overlapScore, k_min=k_min, minval=min_approx)
    		if olen != 0:
    			if a[0] not in overlaps:
    				overlaps[a[0]] = dict()
    			overlaps[a[0]][b[0]] = olen
    		if olen > k_enough:
    			break
    print(f'Building graph finished in {(time.time() - start)/60} minutes.')
    return overlaps


def join_into_contig(a, b, overlap):
	joined_reads = a[:-overlap] + b
	return joined_reads


def find_longest_overlap(ind_a, overlaps):
    laps = overlaps[ind_a]
    if len(laps) > 0:
    	longest = sorted(laps.items(), key=lambda x: x[1])[-1]
    return longest


def create_contigs(reads, overlaps):
	to_delete = []
	for ind_a in overlaps.keys():
		if ind_a in reads.keys():
			ind_b, olap = find_longest_overlap(ind_a, overlaps)
			if ind_b in reads.keys():
				contig = join_into_contig(reads[ind_a], reads[ind_b], olap)
				to_delete.append(ind_a)
				del reads[ind_a]
				reads[ind_b] = contig
	for ind in to_delete:
		del overlaps[ind]
	return reads, overlaps


def remove_biased_reads(reads, init_len):
	""" remove reads that are not in any contig """
	only_contigs = reads.copy()
	for r in reads.items():
		if len(r[1]) <= init_len:
			del only_contigs[r[0]]
	return only_contigs


def write_contigs(reads, filename):
	i = 1
	with open(filename, "w") as contigs:
		for ind, contig in reads.items():
			ctg = ">contig" + str(i) + '\n' + str(contig) + '\n'
			i += 1
			print(ctg)
			contigs.write(ctg)


if __name__ == '__main__':
	if not os.path.isfile(args.input_reads) or not args.input_reads.endswith('.fasta'):
		raise IOError("Given input file does not exist or does not have fasta file extension")
	out_contigs = fasta_format(args.output_contigs)
	minolap = args.minmial_overlap
	start = time.time()
	reads = readin_fasta(args.input_reads)
	init_len = len(reads[0])
	start_all = time.time()
	overlaps = overlap_graph(reads, k_min=minolap, min_approx=10)
	time_gone = (time.time()-start_all)/60
	reads, overlaps = create_contigs(reads, overlaps)
	while time_gone < 55:
		reads, overlaps = create_contigs(reads, overlaps)
		time_gone = (time.time()-start_all)/60
	reads = remove_biased_reads(reads, init_len)
	print(f"Found {len(reads)} contigs")
	write_contigs(reads, out_contigs)

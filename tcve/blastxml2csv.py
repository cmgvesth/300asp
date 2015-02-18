#!/usr/bin/python

'''##################################################################
# Imports
##################################################################'''
import sys,os
sys.path.append(os.path.join(os.path.dirname(__file__), "../utils/"))
from aspmine_imports import *

def blastxml2csv(filename):
	source_exists("file", filename) # costum function, will exit if not exists
	records = SearchIO.parse(filename, 'blast-xml')
	allvalues = []

	for item in records:
		for hsp in item.hsps:

			#filename = (blast.split("/")[-1]).replace(".txt", "")
			#qorg = filename.split("_")[0]
			#sorg = filename.split("_")[-1]

			pident = (float(hsp.ident_num) / float(hsp.aln_span)) * 100
			values = (	filename, hsp.aln_span, hsp.alphabet, hsp.bitscore, hsp.bitscore_raw, hsp.evalue, hsp.gap_num, \
						hsp.hit_end, hsp.hit_frame, hsp.hit_id, str(hsp.hit_range),  \
						hsp.hit_span, hsp.hit_start, hsp.hit_strand, hsp.ident_num, hsp.is_fragmented, hsp.pos_num, hsp.query_end, \
						hsp.query_frame, hsp.query_id, str(hsp.query_range), hsp.query_span, hsp.query_start,hsp.query_strand, pident)

			valuenames = "filename", "aln_span", "alphabet", "bitscore", "bitscore_raw", "evalue", "gap_num",\
						"hit_end", "hit_frame", "hit_id", "hit_range",\
						"hit_span", "hit_start", "hit_strand", "ident_num", "is_fragmented", "pos_num", "query_end",\
						"query_frame", "query_id", "query_range", "query_span", "query_start", "query_strand", "pident"
			allvalues.append(values)

	cursor2csv(valuenames, allvalues, filename+".csv")	
	print "# INFO: blastxml2csv wrote to file %s" % filename+".csv"		

blastxml2csv("Aniger1_vs_Aniger3.txt")
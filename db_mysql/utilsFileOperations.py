#!/usr/bin/python

#------------------------------------------------------------------
# Imports
#------------------------------------------------------------------
import sys, getopt, argparse, re, glob, os
import gzip

#######################################################################
# Decompress TAB files
# Example use: records = gzip_tab_file(filepath)
# Returns a list of lines
#######################################################################

def gzip_tab_file (filepath):
	records = []
	try:
		decompressed_file = gzip.open(filepath, "rb")
		tmp_records = (decompressed_file.read()).split("\n")

	except:
		sys.exit("# ERROR utils gzip: not GZIP file" )

	for i in tmp_records:
		# Exclude comment and empty lines ------------------------------------
		if i.startswith("#") or i == '' or i.startswith("protein"): continue	
		else: records.append(i) 

	return records

#######################################################################
# Decompress fasta files
# Example use: 	records = gzip_fasta_file(filepath) 
# Returns list of FASTA sequence objects
#######################################################################
def gzip_fasta_file (filepath):
	try:
		decompressed_file = gzip.open(filepath, "rb")
		records = list(SeqIO.parse(decompressed_file , "fasta"))
	except:
		sys.exit("# ERROR utils FASTA gzip: not GZIP file" )
	return records
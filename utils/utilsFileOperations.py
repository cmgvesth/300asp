#!/usr/bin/python

#------------------------------------------------------------------
# Imports
#------------------------------------------------------------------
'''##################################################################
# Imports
##################################################################'''
import sys 
from aspmine_imports import *

	
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
		sys.exit("# ERROR utils TAB gzip: not GZIP file" )

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

	# If this is empty, there was nothing that matched a FASTA format in the file	 ------------------------------------
	if len(list(records)) < 1:	sys.exit("# ERROR: Not FASTA format or empty FASTA entry")
	
	return records

#######################################################################
# Verify that file exists, directory exists and is not empty
# Example use: source_exists("dir", source)
# Does not return any values but errors if somthing is wrong
#######################################################################
def source_exists(stype, source):
	allfiles = []
	if stype == "file":
		if os.path.isfile(source):
			print "# INFO: file exists: %s" % source 
		else:
			sys.exit("# ERROR: file does not exists: %s" % source  )

	if stype == "dir":
		# Make sure that the path exists
		if os.path.exists(source):	print "# INFO: path exists: %s" % source 
		else:	sys.exit("# ERROR: path does not exists: %s" % source   )
		
		# Loop through all files in directory and check if it is empty
		for root, subfold, files in os.walk(source):
			allfiles.append(files)
			
		if allfiles == []:
			sys.exit("# ERROR: directory is empty: %s" % source  )
	
				
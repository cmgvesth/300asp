#!/usr/bin/python

#------------------------------------------------------------------
# Imports
#------------------------------------------------------------------
import gzip
import sys, getopt, argparse, subprocess, os

#------------------------------------------------------------------
# Format argument helper
#------------------------------------------------------------------
class SmartFormatter(argparse.HelpFormatter):
	width = 100
	def _split_lines(self, text, width):
	# this is the RawTextHelpFormatter._split_lines
		if text.startswith('R|'):
			return text[2:].splitlines()  
		return argparse.HelpFormatter._split_lines(self, text, width)

#------------------------------------------------------------------
# Get command line arguments
#------------------------------------------------------------------
parser = argparse.ArgumentParser(formatter_class=SmartFormatter, 
								 usage='%(prog)s -source filename\n'
								 "Example: python %(prog)s -source Aspsac1_AssemblyScaffolds.fasta.gz\n")

parser.add_argument('-source' , required=True,help="FASTA file in GZ format")

args = parser.parse_args() # Store argument values
filename = args.source.split("/")[-1]

#------------------------------------------------------------------
# Try to open GZIP file
#------------------------------------------------------------------

try:
	decompressed_file = gzip.open(args.source, "rb")
	print "# Success decompress"
except:
	print "# Error: not GZIP file" 
	sys.exit()

antismash_path = subprocess.check_output(["which run_antismash"],shell=True)

fh = open("tmp.fasta", "w")
fh.write(decompressed_file.read())

decompressed_file.close()

#print filename+".antismash"

if not os.path.exists(filename+".antismash"):

# antismash_path.rstrip() removes trailing line break
	os.system("%s --logfile %s --statusfile %s --all-orfs --outputfolder %s --input-type nucl tmp.fasta" % (antismash_path.rstrip(), filename+".log" ,filename+".status" , filename+".antismash"))
else:
	print "# ERROR: file already exists, delete to rerun: "+filename+".antismash"
	sys.exit()
	
os.system("rm tmp.fasta")



#!/usr/bin/python

#------------------------------------------------------------------
# Imports
#------------------------------------------------------------------
import MySQLdb as mdb
import sys, getopt, argparse, re, glob, os
import gzip
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from mysql_classes import Fasta	# custom class
from mysql_classes import Gff	# custom class

#------------------------------------------------------------------
# Format argument helper
#------------------------------------------------------------------
class SmartFormatter(argparse.HelpFormatter):
	width = 100
	def _split_lines(self, text, width):
	#------------------------------------------------------------------
	# this is the RawTextHelpFormatter._split_lines
	#------------------------------------------------------------------
		if text.startswith('R|'):
			return text[2:].splitlines()  
		return argparse.HelpFormatter._split_lines(self, text, width)

#------------------------------------------------------------------
# Get command line arguments
#------------------------------------------------------------------

parser = argparse.ArgumentParser(formatter_class=SmartFormatter, 
								 usage='%(prog)s -action action -source filename -atype analysis_type\n'
								 "Example: python %(prog)s -dbname aspminedb -action test -atype antismash -source Aspsac1_AssemblyScaffolds.fasta.gz\n"
								 "Example: python %(prog)s -action load -atype antismash -source Aspsac1_AssemblyScaffolds.fasta.gz\n"
								 "Example: python %(prog)s -prefix new2 -action load -atype antismash -source Aspsac1_AssemblyScaffolds.fasta.gz\n"

#------------------------------------------------------------------
# Choose if data should be loaded or only tested
#------------------------------------------------------------------
parser.add_argument('-action' , required=True, default = "test", choices=['load', 'test'],
					help='R|Test or test and load data\n'
					"-test\t: Read files and verify formats\n"
					"-load\t: Read files and verify formats followed by storing data in MySQL tables\n")

#------------------------------------------------------------------
# At least one file filetype must be provided, or the option "all" which expects a JGI directory format
#------------------------------------------------------------------
parser.add_argument('-atype' , required=True, default = "antismash",choices=['antismash'],
					help="R|Select analysis type.\n"
					"-atype pf\t: Antismash annotation\n")

parser.add_argument("-source", required=True, help="R|File or directoryname, example: -file Aspbr1_AssemblyScaffolds.fasta.gz")

parser.add_argument("-dbname", required=True, default = "asminedb", help="R|Database name")

parser.add_argument("-prefix", required=False, help="R|Organism prefix, force new database entry")

args = parser.parse_args()

#------------------------------------------------------------------
# Check prefix validity
#------------------------------------------------------------------
if args.prefix and not re.search(r"^[0-9A-Fa-f]*$", args.prefix):
	print "# ERROR : prefix contains characters other than letters and numbers"
	sys.exit()
#------------------------------------------------------------------
# Print argument values to screen
#------------------------------------------------------------------
print "# Arguments accepted:\n# source\t:%s\n# analysis type\t:%s\n# action\t:%s\n# prefix\t:%s" % (args.source, args.atype, args.action, args.prefix)

#------------------------------------------------------------------
# Connect to specific DB
#------------------------------------------------------------------
try:
    db = mdb.connect("localhost","asp","1234",args.dbname)

except mdb.Error, e:
    print "# Error %d: %s" % (e.args[0],e.args[1])
    sys.exit(1)


if args.atype == "antismash":
	

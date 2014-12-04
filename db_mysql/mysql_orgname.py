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
"""
CODEBOOK INFO:

Program reads a TAB separated file with nice organism names and organism short code.
The codes are compared to the codes in the db and the real names are stored. 

Create source list:

ll | awk -F, 'BEGIN{FS=" "}{print $(NF)}' | grep "^Asp" | awk -F, 'BEGIN{FS="["}{print $1"\t["$2}' > organismnames.txt
head organismnames.txt
Aspergillus_acidus_v1.0_	[Aspfo1]/
Aspergillus_aculeatinus_CBS_121060_v1.0_	[Aspacu1]/
Aspergillus_aculeatus_ATCC16872_v1.1_	[Aspac1]/
Aspergillus_brasiliensis_v1.0_	[Aspbr1]/
Aspergillus_brunneoviolaceus_CBS_621.78_v1.0_	[Aspbru1]/
Aspergillus_carbonarius_ITEM_5010_v3_	[Aspca3]/
Aspergillus_clavatus_NRRL_1_from_AspGD_	[Aspcl1]/
........
"""
			
#------------------------------------------------------------------
# Format argument helper
#------------------------------------------------------------------
class SmartFormatter( argparse.HelpFormatter ):
	width = 100
	def _split_lines( self, text, width ):
	# this is the RawTextHelpFormatter._split_lines
		if text.startswith( 'R|' ):
			return text[2:].splitlines()  
		return argparse.HelpFormatter._split_lines( self, text, width )

#------------------------------------------------------------------
# Get command line arguments
#------------------------------------------------------------------
parser = argparse.ArgumentParser( formatter_class=SmartFormatter, 
								 usage='%(prog)s -action connect\n'
								 "Example: python %(prog)s -dbname aspminedb -source organismnames.txt\n" )

#------------------------------------------------------------------
# Choose if data should be loaded or only tested
#------------------------------------------------------------------

parser.add_argument( "-source",		required=True,	help="R|File name, example: -source organismnames.txt" )
parser.add_argument( "-dbname",		required=False, help="R|Database name", default = "aspminedb" )

args = parser.parse_args()

#------------------------------------------------------------------
# Print argument values to screen
#------------------------------------------------------------------
print "#--------------------------------------------------------------"
print "# ARGUMENTS :\n# source\t:%s\n# database\t:%s" % ( args.source,args.dbname )
print "#--------------------------------------------------------------"

#------------------------------------------------------------------
# Connect to specific DB
#------------------------------------------------------------------
try:
    db = mdb.connect( "localhost","asp","1234",str( args.dbname ) )
    db.close() 
except mdb.Error, e:
    sys.exit( "# ERROR %d: %s" % ( e.args[0],e.args[1] ) )

db = mdb.connect( "localhost","asp","1234",str( args.dbname ) )
cursor = db.cursor()

#------------------------------------------------------------------
# Get lines from file
#------------------------------------------------------------------
file_obj = open( args.source, "r" )
file_lines = file_obj.readlines()
file_obj.close()        

print "# INFO: number of organism in file %s" % len(file_lines)

# Check that only organisms expected are in the file. Also prevents strange files
# These problems SHOULD be prevented in the initial construction of the file, but if the list is manually made this is the first validation
for f in file_lines:	
	if not re.search("^Asp", f) and not re.search("^Nep", f) and not re.search("^Pen", f):
		sys.exit("# ERROR: File contains lines not starting with appropriate species name (Asp, Neo, Pen)")
	
	(org_name, org_code) = (f.split("\t")[0], f.split("\t")[1])
	org_code = org_code.replace("[", "").replace("]", "").replace("/", "").replace("\n", "")
	org_name = org_name.replace("_", " ").replace(" $", "")

	try:	
		cursor.execute( "UPDATE organism SET real_name = '%s' WHERE name = '%s';" % (org_name,org_code))
		print cursor._last_executed
		db.commit()	# Add changes to database
		
	except mdb.Error, e:
		sys.exit( "# ERROR organism name %s %d: %s" % ( args.source, e.args[0],e.args[1] ) )
	

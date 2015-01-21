#!/usr/bin/python

#------------------------------------------------------------------
# Imports
#------------------------------------------------------------------
import MySQLdb as mdb
import sys, getopt, argparse, re, glob, os, gzip
from datetime import datetime

''' bio python '''
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

''' math '''
import csv

''' local libs '''
sys.path.append('/home/tcve/github/utils/')
from utilsArgparse import * # custom functions



##############################################################################################################################################
##############################################################################################################################################
####################################################################### ARGUMENTS ############################################################
##############################################################################################################################################
##############################################################################################################################################
startTime = datetime.now() # record runtime

#------------------------------------------------------------------
''' Get command line arguments '''
#------------------------------------------------------------------
parser = CustomArgumentParser(formatter_class=SmartFormatter, 
								description="Load primary metabolism model data into local mysql database", 
								 usage='%(prog)s -out filename\n'
								 "Example: python %(prog)s -dbname aspminedb -filetype model -action test -source NigerModel2013.csv")
#parser.add_argument("-limit", "-l", required=False, default = 100, help="limit length selects, dafault 100, to select all set to 'all' ")
parser.add_argument("-out", "-o", required=False, default = "lengths.tab", help="Name of output file")
parser.add_argument("-dbname", "-d", required=False, default = "aspminedb", help="Database name")
parser.add_argument("-sec", required=False, help="Section name")
parser.add_argument("-plot", "-p", required=False, help="Create")

#parser.add_argument("-verbose", "-v", required=False, help="R|Increase output verbosity for bug tracking" , action='store_true')

args = parser.parse_args()
section = args.sec

#------------------------------------------------------------------
# Print argument values to screen
#------------------------------------------------------------------
print "#--------------------------------------------------------------"
print "# ARGUMENTS :\n# Database\t:%s\n# Outfile\t:%s" % (args.dbname, args.out)
print "#--------------------------------------------------------------"

""" TO DO: 
make sure that the R script is in the right location and make some simple grep to veryfy what is in it 
else - exit
"""
outfile = open(args.out, 'w') 

##############################################################################################################################################
##############################################################################################################################################
####################################################################### DATABASE  ############################################################
##############################################################################################################################################
##############################################################################################################################################

#------------------------------------------------------------------
# Connect to specific DB
#------------------------------------------------------------------

try:
    db = mdb.connect("localhost","asp","1234",str(args.dbname))

except mdb.Error, e:
    sys.exit("# ERROR %d: %s" % (e.args[0],e.args[1]))

try:
	cursor = db.cursor()

except mdb.Error, e: 
	sys.exit("# ERROR %d: %s" % (e.args[0],e.args[1]))

##############################################################################################################################################
##############################################################################################################################################
####################################################################### DATA      ############################################################
##############################################################################################################################################
##############################################################################################################################################

#------------------------------------------------------------------
# Create table protein lengths
#------------------------------------------------------------------
try:
	cursor.execute("drop table IF EXISTS lengths ;")
	cursor.execute("create table lengths as select org_id, prot_seqname, length(prot_seq) as len from proteins;")

except mdb.Error, e:
	sys.exit("# ERROR %d: %s" % (e.args[0],e.args[1]))

#------------------------------------------------------------------
# Get protein lengths
# Using panda sql to get data from the database into a R readable format
# psql.read_frame
#------------------------------------------------------------------
result = ""

try: 
	if section:
		cursor.execute("SELECT organism.org_id, section, name, real_name, len FROM lengths join organism using (org_id) where section = \'" + section + "\';")
	else:
		cursor.execute("SELECT organism.org_id, section, name, real_name, len FROM lengths join organism using (org_id);")

except mdb.Error, e:
	sys.exit("# ERROR %d: %s" % (e.args[0],e.args[1]))

result = cursor.fetchall()

if not type(result) is tuple:
	sys.exit("# ERROR : result of query is not a tuple but a %s" % type(result))

##############################################################################################################################################
##############################################################################################################################################
####################################################################### PLOT      ############################################################
##############################################################################################################################################
##############################################################################################################################################
f = open("tmp.csv",'wb')
writer = csv.writer(f, dialect = 'excel')
writer.writerows(result)
f.close

if args.plot:
	os.system("R CMD BATCH aspmine_analysis_lengths.R")







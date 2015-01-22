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
parser.add_argument("-out", "-o", required=False, default = "lengths.csv", help="Name of output file")
parser.add_argument("-dbname", "-d", required=False, default = "aspminedb", help="Database name")
parser.add_argument("-sec", required=False, help="Section name")
parser.add_argument("-plot", "-p", default = "n", required=False, help="Create plot", choices=['n', 'y'])
parser.add_argument("-tab", "-t", default = "n", required=False, help="Create table", choices=['n', 'y'])
parser.add_argument("-R", default="aspmine_analysis_lengths.R", required=False, help="Name of R script, default aspmine_analysis_lengths.R")

#parser.add_argument("-verbose", "-v", required=False, help="R|Increase output verbosity for bug tracking" , action='store_true')

args 	= parser.parse_args()
section = args.sec
plot 	= args.plot
table 	= args.tab
rscript = args.R
outfile = args.out

#------------------------------------------------------------------
# Print argument values to screen
#------------------------------------------------------------------
print "#--------------------------------------------------------------"
print "# ARGUMENTS :\n# Database\t:%s\n# Outfile\t:%s \n# Create Plot\t:%s \n# Create Table\t:%s \n# Rscript\t:%s" % (args.dbname, args.out, plot, table, rscript)
print "#--------------------------------------------------------------"

if not os.path.isfile(rscript):
	sys.exit("# ERROR: Rscript file does not exist, %s" % rscript)

if table=="n": 
	if not os.path.isfile(outfile):
		sys.exit("# ERROR: length file does not exist, %s , run with -tab option " % outfile)

#outfile = open(args.out, 'w') 

##############################################################################################################################################
##############################################################################################################################################
####################################################################### FUNCTIONS ############################################################
##############################################################################################################################################
##############################################################################################################################################

##############################################################################################################################################
####################################################################### DATABASE  ############################################################
##############################################################################################################################################
def make_table(outfile):
	print "# INFO: creating length table"

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

	f = open(outfile,'wb')
	writer = csv.writer(f, dialect = 'excel')
	writer.writerows(result)
	f.close()

##############################################################################################################################################
####################################################################### PLOT      ############################################################
##############################################################################################################################################
def make_plot(outfile):
	print "# INFO: running Rscript"
	os.system("R CMD BATCH '--args %s' aspmine_analysis_lengths.R test.out " % outfile)


if table=="y":
	make_table(outfile)
	if plot=="y":
		make_plot(outfile)

if table=="n" and plot=="y":
	make_plot(outfile)







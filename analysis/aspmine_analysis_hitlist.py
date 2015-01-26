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

import warnings
warnings.filterwarnings("ignore", "Unknown table.*")


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
								description="Get best hits from specific gene name or ID", 
								 usage='%(prog)s -dbname [database name] -out [output filename] -gene [JGI gene name or ID]\n'
								 "Example: python %(prog)s -dbname aspminedb -out hitlist.csv -gene An15g00560")
parser.add_argument("-out", "-o", required=False, default = "hitlist.csv", help="Name of output file, default hitlist.csv")
parser.add_argument("-tab", "-t", required=False, help="Create table new or use existing tmp")
parser.add_argument("-dbname", "-d", required=False, default = "aspminedb", help="Database name")
parser.add_argument("-gene", "-g", nargs = '*',  required=True, default=[None], action='store', help="List of gene IDs or names")
parser.add_argument("-strain", "-s", nargs = '*',  required=True, default=[None], action='store', help="List of organism JGI keys")

args 	= parser.parse_args()
genes 	= args.gene
strains = args.strain
outfile = args.out

#------------------------------------------------------------------
# Print argument values to screen
#------------------------------------------------------------------
print "#--------------------------------------------------------------"
print "# ARGUMENTS :\n# Database\t:%s\n# Outfile\t:%s \n# Genes\t:%s \n# Strains\t:%s" % (args.dbname, outfile, genes, strains)
print "#--------------------------------------------------------------"

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
####################################################################### FUNCTIONS ############################################################
##############################################################################################################################################
##############################################################################################################################################

##############################################################################################################################################
####################################################################### make_database  #######################################################
##############################################################################################################################################
def make_database(strains, genes):
	print "# INFO: creating subset table"

	#------------------------------------------------------------------
	# This table includes two way hits, as it gets both sseq and qseq data
	# These hits could be the same but in rare cases they might differ
	# Common case: A -> B AND B -> A # not nessesaraly with the same score
	# Rare case: A -> B BUT B ! -> ! A
	#------------------------------------------------------------------
	try:
		cursor.execute("drop table IF EXISTS tmp ;")
		cursor.execute("create table tmp as \
						select a.* from blast a \
						left join blast b \
						on b.blast_qseq_id=a.blast_sseq_id \
						and b.blast_sseq_id=a.blast_qseq_id \
						and b.blast_qseq_id<>b.blast_sseq_id \
						where ( " + gene_string + ") \
						and (" + strain_string + ") \
						and (b.blast_qseq_id is null or b.blast_qseq_id >= a.blast_qseq_id);")

	except mdb.Error, e:
		sys.exit("# ERROR %d: %s" % (e.args[0],e.args[1]))

	cursor.execute("SELECT count(*) FROM tmp")
	result = cursor.fetchall()

	if not type(result) is tuple:
		sys.exit("# ERROR : result of query is not a tuple but a %s" % type(result))

	#------------------------------------------------------------------
	# Get protein lengths
	# Using panda sql to get data from the database into a R readable format
	# psql.read_frame
	#------------------------------------------------------------------
##############################################################################################################################################
####################################################################### get_data  ############################################################
##############################################################################################################################################

def get_data (strains, genes):

	result = ""

	try: 
		cursor.execute("select blast_sseq_id,blast_qseq_id,blast_pident from tmp r \
						where blast_pident=(select max(cast( blast_pident as DECIMAL(5,2))) from tmp s \
						where s.blast_qseq_jg1=r.blast_qseq_jg1 and s.blast_sseq_jg1=r.blast_sseq_jg1  \
						group by s.blast_sseq_jg1,s.blast_qseq_jg1);")
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
##############################################################################################################################################
####################################################################### MAIN #################################################################
##############################################################################################################################################
##############################################################################################################################################

#------------------------------------------------------------------
# Create strain and gene query strings
#------------------------------------------------------------------
(strain_string, gene_string) = ("","")

for s in range(0,len(strains)):
	strain_string += " b.blast_qseq_jg1 = \'" + strains[s] + "\' or b.blast_sseq_jg1 = \'" + strains[s] + "\'"
	if s != len(strains)-1:
		strain_string += " or "


for i in range(0,len(genes)):
	gene_string += " a.blast_qseq_id like \'%" + genes[i] + "%\' or a.blast_sseq_id like \'%" + genes[i] + "%\'"
	if i != len(genes)-1:
		gene_string += " or "

#------------------------------------------------------------------
# Check if table tmp already contains the needed genesa and strains 
#------------------------------------------------------------------

#------------------------------------------------------------------
# create table if specified or if tmp does not contain the right genes and strains
#------------------------------------------------------------------
if table :
	make_database(strains, genes)

#------------------------------------------------------------------
# Get data from tmp table
#------------------------------------------------------------------
get_data(strains, genes)



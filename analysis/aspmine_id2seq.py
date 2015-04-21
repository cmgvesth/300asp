#!/usr/bin/python

'''##################################################################
# Imports
##################################################################'''
import sys,os
root = os.path.abspath(__file__).split("github")[0]
sys.path.append(root+"github/utils/")
from aspmine_imports import *

'''##################################################################
# ARGUMENTS and setup
##################################################################'''
startTime = datetime.now() # record runtime
parser = CustomArgumentParser(formatter_class=SmartFormatter, 
								description="Get best hits from specific gene name or ID", 
								 usage='%(prog)s -dbname [database name] -out [output filename] -gene_strain [series of gene ids and JGI organism codes] \n'
								 "Example: python %(prog)s -dbname aspminedb -out hitlist.csv -gene_strain An15g00560 Aspni_DSM_1 1123159 Aspni7")
parser.add_argument("-out", "-o", required=False, default = "hitlist.csv", help="Name of output file, default hitlist.csv")
parser.add_argument("-tab", "-t", required=False, action='store_true', help="Create table new or use existing tmp")
parser.add_argument("-dbname", "-d", required=False, default = "aspminedb", help="Database name")
parser.add_argument("-gene_strain", "-gs",  nargs = '*', required=True, default=[None], action='store', help="List of gene ID or name")

args 	= parser.parse_args()
gene_strain = args.gene_strain
outfile = args.out

#------------------------------------------------------------------
#------------------------------------------------------------------
sys.exit("# Script is currently under development!")
#------------------------------------------------------------------
#------------------------------------------------------------------

'''------------------------------------------------------------------
# Print argument values to screen
------------------------------------------------------------------'''
print "#--------------------------------------------------------------\n\
# ARGUMENTS\n\
# Description\t: Script return list of best bidirectional hit for a gene in a specific strain\n\
# Database\t: %s\n\
# Outfile\t: %s\n\
# Gene\t\t: %s\n\
# Strain\t: %s\n\
# Table\t\t: %s\n\
#--------------------------------------------------------------" % (args.dbname, outfile, gene, strain, table)

'''------------------------------------------------------------------
# Connect to specific DB
------------------------------------------------------------------'''
try:
    db = mdb.connect("localhost","asp","1234",str(args.dbname))
except mdb.Error, e:
    sys.exit("# ERROR %d: %s" % (e.args[0],e.args[1]))

try:
	cursor = db.cursor()
except mdb.Error, e: 
	sys.exit("# ERROR %d: %s" % (e.args[0],e.args[1]))

#!/usr/bin/python

'''##################################################################
# Imports
##################################################################'''
import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), "../utils/"))
import MySQLdb as mdb
import getopt, argparse, re, glob, os, gzip
from datetime import datetime
from utilsArgparse import * # custom functions


"""
import MySQLdb as mdb
db = mdb.connect(host="192.38.13.9", port= 3306 ,user="setd",passwd="1234",db="testasp")
passwd does not include a hyphen
"""

##############################################################################################################################################
##############################################################################################################################################
####################################################################### ARGUMENTS ############################################################
##############################################################################################################################################
##############################################################################################################################################
startTime = datetime.now() # record runtime


#------------------------------------------------------------------
# Get command line arguments
#------------------------------------------------------------------
parser = CustomArgumentParser(formatter_class=SmartFormatter, 
								description="Load data into local mysql database", 
								 usage='%(prog)s -dbname [databasename] \n'
								 "Example: python %(prog)s -dbname aspminedb\n"
								 "Example: python %(prog)s -dbname aspminedb")
#------------------------------------------------------------------
# Choose if data should be loaded or only tested
# At least one file filetype must be provided, or the option "all" which expects a JGI directory format
#------------------------------------------------------------------
parser.add_argument("-dbname", "-d", required=False, default = "aspminedb", help="R|Database name")
parser.add_argument("-host", required=False, default = "localhost", help="R|Database host IP")

args = parser.parse_args()
host = args.host
db = args.dbname

#------------------------------------------------------------------
# Print argument values to screen
#------------------------------------------------------------------
print "#--------------------------------------------------------------"
print "# ARGUMENTS :\n# Database\t:%s\n# Host\t:%s" % (args.dbname, args.host)
print "#--------------------------------------------------------------"

#------------------------------------------------------------------
# Connect to specific DB
#------------------------------------------------------------------

def database( db, host ) :
	try:
	    db = mdb.connect(host,"asp","1234",db)
	    print "# INFO: success connectiong to existing database %s" % dbname

	except mdb.Error, e:
		sys.exit("# ERROR %d: %s" % (e.args[0],e.args[1]))

database( str(db), str(host) )
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
								description="Test connection to mysql database", 
								 usage='%(prog)s -dbname [databasename] -user [mysql username] -p [Mysql password] -d [mysql databasename]\n'
								 "Example: python %(prog)s -host 192.38.13.9 -user asp -password 1234 -d aspminedb")
#------------------------------------------------------------------
# Choose if data should be loaded or only tested
# At least one file filetype must be provided, or the option "all" which expects a JGI directory format
# klara: 192.38.13.9
#------------------------------------------------------------------
parser.add_argument("-dbname", "-d", required=True, help="R|Database name")
parser.add_argument("-host", required=True,  help="R|Database host IP")
parser.add_argument("-user", required=True, help="R|Database username")
parser.add_argument("-password", "-p", required=True, help="R|Database password")

args = parser.parse_args()
host = args.host
dbname = args.dbname
password = args.password
user = args.user

#------------------------------------------------------------------
# Print argument values to screen
#------------------------------------------------------------------
print "#--------------------------------------------------------------"
print "# ARGUMENTS :\n\
# Database\t:%s\n\
# Host\t\t:%s\n\
# Password\t:%s\n\
# Username\t:%s" % (dbname, host, password, user)
print "#--------------------------------------------------------------"

#------------------------------------------------------------------
# Connect to specific DB
#------------------------------------------------------------------

def database( db, host ) :
	try:
	    db = mdb.connect(host,"asp","1234",dbname)

	    print "# INFO: success connecting to existing database %s" % dbname

	except mdb.Error, e:
		sys.exit("# ERROR %d: %s" % (e.args[0],e.args[1]))

database( str(dbname), str(host) )

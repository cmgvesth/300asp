#!/usr/bin/python

'''##################################################################
# Imports
##################################################################'''
import sys 
sys.path.append(os.path.join(os.path.dirname(__file__), "../utils/"))
from aspmine_imports import *


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
								description="Load primary metabolism model data into local mysql database", 
								 usage='%(prog)s -filetype filetype -action action -source filename\n'
								 "Example: python %(prog)s -dbname aspminedb -filetype model -action test -source NigerModel2013.csv")
#------------------------------------------------------------------
# Choose if data should be loaded or only tested
#------------------------------------------------------------------
parser.add_argument('-action' , "-a", required=True, default = "test", choices=['load', 'test'],
					help='R|Test or test and load data\n'
					"-test\t: Read files and verify formats\n"
					"-load\t: Read files and verify formats followed by storing data in MySQL tables\n")

#------------------------------------------------------------------
# At least one file filetype must be provided, or the option "all" which expects a JGI directory format
#------------------------------------------------------------------
parser.add_argument('-filetype' ,"-f",  required=True, default = "model", choices=["model"],
					help="R|Select filetype: model, CSV dump of excel document")

parser.add_argument("-source", "-s", required=True, help="R|File , example: -source NigerModel2013.csv")#, type=argparse.FileType('r'))
parser.add_argument("-dbname", "-d", required=False, default = "aspminedb", help="R|Database name")
parser.add_argument("-orgkey", "-o",required=False, help="R|Organism key/name")
parser.add_argument("-verbose", "-v", required=False, help="R|Increase output verbosity for bug tracking" , action='store_true')

args = parser.parse_args()

#------------------------------------------------------------------
# Print argument values to screen
#------------------------------------------------------------------
print "#--------------------------------------------------------------"
print "# ARGUMENTS :\n# source\t:%s\n# filetype\t:%s\n# action\t:%s\n# Database\t:%s\n# Organism key\t:%s" % (args.source, args.filetype, args.action, args.dbname, args.orgkey)
print "#--------------------------------------------------------------"

#------------------------------------------------------------------
# Test functions
#------------------------------------------------------------------
if args.verbose:
	test__get_org_id()
	test__fasta()

	#test__gff()
	#test__ipr()
	#test__kegg()
	#test__go()
	#test__sigp()
	#test__kog()

	sys.exit("# FINISHED testing functions")

#------------------------------------------------------------------
# Connect to specific DB
#------------------------------------------------------------------
try:
    db = mdb.connect("localhost","asp","1234",str(args.dbname))

except mdb.Error, e:
    sys.exit("# ERROR %d: %s" % (e.args[0],e.args[1]))

##############################################################################################################################################
##############################################################################################################################################
####################################################################### Main #################################################################
##############################################################################################################################################
##############################################################################################################################################
source = args.source
filetype = args.filetype
action = args.action
dbname = args.dbname

#-----------------------------------
# Read in data
# select first row as titles
# replace leading spaces
#-----------------------------------
records = []
try:
	fileobj = open(source, "r")
	tmp_records = (fileobj.read()).split("\n")
except:
	sys.exit("# ERROR utils TAB gzip: not GZIP file" )

for i in tmp_records:
	if tmp_records.index(i) == 0:
		title_line = i
	else:
		tmp_record = i.split("\t")
		records.append([t.lstrip(" ") for t in tmp_record]) 

#-----------------------------------
# Identify reaction title lines
# Get upper and lower level titles
#-----------------------------------
(upper_flag, lower_flag) = (0,0)

for r in records:
	#if r.count("") == 17 or r.count("") == 16:
	if 
		if upper_flag == 0:
			upper_level = r[1]
			upper_flag = 1

		if upper_flag == 1:
			lower_level = r[1]
			upper_flag == 0
		
		print upper_level, lower_level




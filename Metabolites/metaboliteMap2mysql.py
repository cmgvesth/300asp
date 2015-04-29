#!/usr/bin/python

'''------------------------------------------------------------------
# Imports
------------------------------------------------------------------'''
import sys,os
root = os.path.abspath(__file__).split("github")[0]
sys.path.append(root+"github/utils/")
from aspmine_imports import *

'''------------------------------------------------------------------
# Code uses local functions: 
------------------------------------------------------------------

CustomArgumentParser
DBconnect
executeQuery

------------------------------------------------------------------'''

'''------------------------------------------------------------------
# ARGUMENTS and setup
------------------------------------------------------------------'''
startTime = datetime.now() # record runtime
parser = CustomArgumentParser(formatter_class=SmartFormatter, usage='%(prog)s -infile [name of csv file]')
parser.add_argument("-dbname", "-d", required=False, default = "aspminedb", help="Database name")

""" INPUT FILES """
parser.add_argument("-infile", "-if", required=False, default = "Aspergillus section Nigri-metab_JCF.csv", help="CSV file from excel of genomes and metabolits")
# https://files.podio.com/152765850
# Aspergillus section Nigri-metab_JCF.csv

""" PARSE ARGUMENTS """
args 	= parser.parse_args()
infile	= args.infile
startTime = datetime.now() # record runtime

'''------------------------------------------------------------------
# Print argument values to screen
------------------------------------------------------------------'''
print "#--------------------------------------------------------------\n\
# ARGUMENTS\n\
# Database\t: %s\n\
# CSV file from excel of genomes and metabolits\t: %s\n\
#--------------------------------------------------------------" % (args.dbname, infile)
# Compile new antismash2organism -clean1\t: %s\n\
# Analysis, cluster VS organism - with best candidate, looping\t\t: %s\n\

csvfile = open(infile, "r")

orgnames = []
vectors = []
data = {}
csvlines = csv.reader(csvfile, delimiter=';')

met1 = next(csvlines)	# First row - metabolite upper level desc
met2 = next(csvlines)	# Second row - metabolite lower level desc

for line in csvlines:
	l = [ e for e in line if e == "" ] # Remove empty lines
	if not len(l) == len(line):	
		orgkey = line[0]
		orgnames.append( orgkey )
		data[ orgkey ] = line[1:]

csvfile.close()

'''------------------------------------------------------------------
# Connection to database
------------------------------------------------------------------'''
cursor = DBconnect("localhost", "asp", "1234", args.dbname) # costum function

'''------------------------------------------------------------------
# Connect orgnames to aspminedb organism table
------------------------------------------------------------------'''
result = {}

for org in orgnames:
	idflag = 0;
	words = org.split(" ")
	for word in words:
		if idflag == 0:
			query = "SELECT name, org_id, real_name from organism where real_name like \'%" + word + "%\'"
			cursor.execute(query)
			(columns, orginfo) = executeQuery(cursor, query)
			if len(orginfo) == 1:
				idflag = 1
				result[orginfo] = data[org][5:] 
		
'''------------------------------------------------------------------
# Create metaboliteMap table
------------------------------------------------------------------'''

startTime_metaboliteMap = datetime.now()
print "# INFO: Dropping and re-creating metaboliteMap"

cursor.execute("DROP TABLE IF EXISTS metaboliteMap")

query="CREATE TABLE metaboliteMap (\
		org_id varchar(100) NOT NULL,\
		name varchar(100) NOT NULL,\
		real_name varchar(100) NOT NULL,\
		metabolite1 varchar(100) NOT NULL,\
		metabolite2 varchar(100) NOT NULL,\
		status varchar(100),\
		PRIMARY KEY (real_name, metabolite1,metabolite2, status)\
		) ENGINE=MyISAM DEFAULT CHARSET=latin1;"
executeQuery(cursor, query)

# If not specified to create table test if the table actually does exist
if not	cursor.execute("Show tables LIKE 'metaboliteMap'"):
	sys.exit("# ERROR: table does not exist, metaboliteMap")

print "# INFO Runtime metaboliteMap: ", (datetime.now()-startTime_metaboliteMap)	

'''------------------------------------------------------------------
# Insert values into metaboliteMap table
------------------------------------------------------------------'''

for r in result:
	i = 0
	for m in result[r]:
		values = (r[0][0],r[0][1],r[0][2], m, met1[6+i], met2[6+i])
		values = list(values)
		values = [ "NULL"  if e=='' else e for e in values ]
		values = [ str(e).replace("xxxx", "NULL") for e in values]
		values = tuple(values)
		query = ("INSERT INTO metaboliteMap (name, org_id, real_name, status, metabolite1, metabolite2)\
				VALUES (" +("'%s'," * len(values)).rstrip(",") + ");" )
		executeQuery(cursor, query % values)
		i += 1



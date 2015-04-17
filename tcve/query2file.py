#!/usr/bin/python

'''##################################################################
# Imports
##################################################################'''
import sys,os
sys.path.append(os.path.join(os.path.dirname(__file__), "../utils/"))
from aspmine_imports import *

'''##################################################################
# ARGUMENTS and setup
##################################################################'''
startTime = datetime.now() # record runtime
parser = CustomArgumentParser(formatter_class=SmartFormatter, usage='%(prog)s -dbname [database name]')
parser.add_argument("-dbname", "-d", required=False, default = "aspminedb", help="Database name")
parser.add_argument("-qfile", "-qf", required=False, help="File containing mysql query")
parser.add_argument("-tqfile", "-tqf", required=False, help="File containing mysql test query")
parser.add_argument("-fqfile", "-fqf", required=False, help="File containing mysql query - result is sequence and should be converted to FASTA")
#parser.add_argument("-stats", "-s", required=False, action='store_true', help="Compile stats table")

args 		= parser.parse_args()
queryfile 	= args.qfile
testfile 	= args.tqfile
fastafile 	= args.fqfile
startTime 	= datetime.now() # record runtime

if fastafile and not queryfile:
	queryfile=fastafile
if not fastafile and not queryfile and not testfile:
	sys.exit("# ERROR: no arguments given")

'''------------------------------------------------------------------
# Print argument values to screen
------------------------------------------------------------------'''
print "#--------------------------------------------------------------\n\
# ARGUMENTS\t\t:\n\
# Database\t\t: %s\n\
# Query file\t\t: %s\n\
# Test query file\t: %s\n\
# FASTA query file\t: %s\n\
#--------------------------------------------------------------" % (args.dbname, queryfile, testfile, fastafile)

'''------------------------------------------------------------------
# Get queries from files
------------------------------------------------------------------'''
(query, tquery, fquery) = (None, None, None)

source_exists("file", queryfile)
qf = open(queryfile, 'r')
query =  qf.read()

if testfile:
	source_exists("file", testfile)
	tqf = open(testfile, 'r')
	tquery = tqf.read()

if fastafile:
	source_exists("file", fastafile)
	fqf = open(fastafile, 'r')
	fquery = fqf.read()

'''------------------------------------------------------------------
# Print queriesv to screen
------------------------------------------------------------------'''
print "\
# Query\t: \n%s\n\
# Test query\t: \n%s\n\
# FASTA query\t: \n%s\n" % (query, tquery, fquery)

'''------------------------------------------------------------------
# Connection to database
------------------------------------------------------------------'''
cursor = DBconnect("localhost", "asp", "1234", args.dbname) # costum function

'''------------------------------------------------------------------
# Execute main query
------------------------------------------------------------------'''
(columns, result)	 = executeQuery(cursor, query)
if columns and result:
	cursor2csv(columns, result, queryfile+"_main.csv") # costum function

'''------------------------------------------------------------------
# Execute test query
------------------------------------------------------------------'''
if tquery:
	(columns, result)	 = executeQuery(cursor, tquery)
	if columns and result:
		cursor2csv(columns, result, testfile+"_test.csv") # costum function

'''------------------------------------------------------------------
# Execute FASTA query
------------------------------------------------------------------'''
if fquery:
	(columns, result) = executeQuery(cursor, fquery) # costum function
	if columns and result:
		f = open(fastafile+".fsa",'wb')
		
		cursor2csv(columns, result, "tmp.csv")
		seqcsv = open("tmp.csv", 'r')
		seqs = csv.reader(seqcsv,delimiter=';')
		seqs.next() # skip first line (header)

		for row in seqs:
			seq = row[-1]
			comment = "_".join( row[0:len(row)-1] )
			comment = comment.replace(" ", "_")
			f.write(">" + comment + "\n")
			f.write(seq + "\n")
	f.close()		

#!/usr/bin/python

'''##################################################################
# Imports
##################################################################'''
import sys 
sys.path.append('../utils/')
from aspmine_imports import *

'''##################################################################
# ARGUMENTS and setup
##################################################################'''
startTime = datetime.now() # record runtime
parser = CustomArgumentParser(formatter_class=SmartFormatter, 
								description="Get best hits from specific gene name or ID", 
								 usage='%(prog)s -dbname [database name] -out [output filename] -gene [JGI gene name or ID] -strain [JGI organism keys] \n'
								 "Example: python %(prog)s -dbname aspminedb -out hitlist.csv -gene An15g00560 -strain Aspni_DSM_1")
parser.add_argument("-out", "-o", required=False, default = "hitlist.csv", help="Name of output file, default hitlist.csv")
parser.add_argument("-tab", "-t", required=False, action='store_true', help="Create table new or use existing tmp")
parser.add_argument("-dbname", "-d", required=False, default = "aspminedb", help="Database name")
parser.add_argument("-strain", "-s", nargs = '*',  required=True, default=[None], action='store', help="JGI organism keys")
parser.add_argument("-gene", "-g", nargs = '*',  required=True, default=[None], action='store', help="List of gene IDs or names")

args 	= parser.parse_args()
genes 	= args.gene
strains = args.strain
outfile = args.out
table 	= args.tab

'''------------------------------------------------------------------
# Print argument values to screen
------------------------------------------------------------------'''
print "#--------------------------------------------------------------\n\
# ARGUMENTS\n\
# Database\t: %s\n\
# Outfile\t: %s\n\
# Genes\t\t: %s\n\
# Strains\t: %s\n\
# Table\t\t: %s\n\
#--------------------------------------------------------------" % (args.dbname, outfile, genes, strains, table)

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

##################################################################
# FUNCTIONS
##################################################################


##################################################################
def make_database(strain_string, gene_string):
##################################################################
	print "# INFO: creating subset table"

	strain_string	= strain_string.replace("blast", "a.blast")
	gene_string		= gene_string.replace("blast", "a.blast")
	
	'''------------------------------------------------------------------
	# This table includes two way hits, as it gets both sseq and qseq data
	# These hits could be the same but in rare cases they might differ
	# Common case: A -> B AND B -> A # not nessesaraly with the same score
	# Rare case: A -> B BUT B ! -> ! A
	------------------------------------------------------------------'''

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
		sys.exit("# ERROR: executing following query failed:\n%s\n # ERROR %d: %s" % (cursor._last_executed, e.args[0],e.args[1]))

	cursor.execute("SELECT count(*) FROM tmp")
	result = cursor.fetchall()

	if not type(result) is tuple:
		sys.exit("# ERROR : result of query is not a tuple but a %s" % type(result))

	print "# INFO: created table tmp using gene string:\n# %s\n# and strain string:\n# %s" % (gene_string, strain_string)
		
##################################################################
def get_data ():
##################################################################
	result = ""

	try: 
		cursor.execute("select blast_sseq_id,blast_qseq_id,blast_pident from tmp r \
						where blast_pident=(select max(cast( blast_pident as DECIMAL(5,2))) from tmp s \
						where s.blast_qseq_jg1=r.blast_qseq_jg1 and s.blast_sseq_jg1=r.blast_sseq_jg1  \
						and s.blast_qseq_jg1 != s.blast_sseq_jg1\
						group by s.blast_sseq_jg1,s.blast_qseq_jg1);")
	except mdb.Error, e:
		sys.exit("# ERROR: executing following query failed:\n%s\n # ERROR %d: %s" % (cursor._last_executed, e.args[0],e.args[1]))

	result = cursor.fetchall()

	if not type(result) is tuple:
		sys.exit("# ERROR : result of query is not a tuple but a %s" % type(result))

	if len(result) < 1:
		sys.exit("# ERROR: query returned no rows, verify strains and gene names")	
	#print "# INFO: get_data returns %s rows from tmp" % len(result)	

	f = open(outfile,'wb')
	writer = csv.writer(f, dialect = 'excel', delimiter=";")
	writer.writerows(result)
	f.close()

	return result

##################################################################
def test_table(strain_string, gene_string):
##################################################################
	'''------------------------------------------------------------------
	# Check if table tmp already contains the needed genesa and strains 
	# The user might think that the table is already there.
	# For the case where that is wrong, error and tell to run new table
	------------------------------------------------------------------'''


	try: 
		cursor.execute("SELECT blast_sseq_id,blast_qseq_id FROM tmp WHERE ( " + gene_string + ") AND (" + strain_string + ") ;")
	except mdb.Error, e:
			sys.exit("# ERROR: executing following query failed:\n# %s\n# ERROR %d: %s" % (cursor._last_executed, e.args[0],e.args[1]))

	result = cursor.fetchall()

	if len(result) < 1:
		if not table:
			sys.exit("# ERROR: table tmp does not contain the needed genes/strains\n# Run with -t option to redo the table")
		else:
			print "# INFO: table tmp does not contain the needed genes/strains"	
	else:
		print "# INFO: table tmp was found to contain entries for strain and genes"		
			
##################################################################
# MAIN
##################################################################

'''------------------------------------------------------------------
# Create strain and gene query strings
------------------------------------------------------------------'''
(strain_string, gene_string) = ("","")

for s in range(0,len(strains)):
	strain_string += "blast_qseq_jg1 = \'" + strains[s] + "\' or blast_sseq_jg1 = \'" + strains[s] + "\'"
	if s != len(strains)-1:
		strain_string += " or "


for i in range(0,len(genes)):
	gene_string += "blast_qseq_id like \'%" + genes[i] + "%\' or blast_sseq_id like \'%" + genes[i] + "%\'"
	if i != len(genes)-1:
		gene_string += " or "

'''------------------------------------------------------------------
# TEST content of existing tmp table
# create table if tmp does not contain the right genes and strains
------------------------------------------------------------------'''
if table:
	print "# INFO: Creating table tmp"
	make_database(strain_string, gene_string)

if cursor.execute("SHOW TABLES LIKE 'tmp' ;"):
	test_table(strain_string, gene_string)
else:
	print "# INFO: Table tmp does not exist, will create table"
	make_database(strain_string, gene_string)

'''------------------------------------------------------------------
# Get data from tmp table
------------------------------------------------------------------'''
result = get_data()
nr_rows = len(result)

percent = [float(row[2]) for row in result ]
mean_id = float(sum(percent))/len(percent) if len(percent) > 0 else float('nan')

min_id	= min(percent)
max_id	= max(percent)

min_string = [row[0:2] for row in result if float(row[2]) == min_id]
max_string = [row[0:2] for row in result if float(row[2]) == max_id]

min_string = ["%s -- VS -- %s" % tup for tup in min_string]
max_string = ["%s -- VS -- %s" % tup for tup in max_string]
'''------------------------------------------------------------------
# Final information to user
------------------------------------------------------------------'''
print "#--------------------------------------------------------------\n\
# ARGUMENTS\n\
# Runtime\t\t: %s\n\
# Nr. rows in table\t: %s\n\
# Average percent id\t: %s\n\
# Maximum percent id\t: %s \n\
# Comparisons with max\t: %s\n\
# Minimum percent id\t: %s \n\
# Comparisons with min\t: %s\n\
#--------------------------------------------------------------" % \
( (datetime.now()-startTime), nr_rows, mean_id, max_id, '\n#\t\t\t: '.join(map(str, max_string)), min_id, '\n#\t\t\t: '.join(map(str, min_string)))

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
''' Get command line arguments '''
#------------------------------------------------------------------
parser = CustomArgumentParser(formatter_class=SmartFormatter, 
								description="Create protein length tale and plots", 
								 usage='%(prog)s -out filename -dbname [database name] -sec [taxonimic section] -plot [draw plots n/y] -tab [create table n/y] -R [name of R script]\n'
								 "Example: python %(prog)s -dbname aspminedb -plot y -tab n -out length.csv [uses existing datafile]")
#parser.add_argument("-limit", "-l", required=False, default = 100, help="limit length selects, dafault 100, to select all set to 'all' ")
parser.add_argument("-out", "-o", required=False, default = "lengths.csv", help="Name of output file")
parser.add_argument("-dbname", "-d", required=False, default = "aspminedb", help="Database name")
parser.add_argument("-sec", required=False, help="Section name")
parser.add_argument("-plot", "-p", action='store_true', required=False, help="Create plot")
parser.add_argument("-tab", "-t", action='store_true', required=False, help="Create table")
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
print "# ARGUMENTS :\n# Database\t:%s\n# Outfile\t:%s \n# Create Plot\t:%s \n# Create Table\t:%s \n# Rscript\t:%s" % (args.dbname, outfile, plot, table, rscript)
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







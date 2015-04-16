#!/usr/bin/python

'''------------------------------------------------------------------
# Imports
------------------------------------------------------------------'''
import sys,os
sys.path.append(os.path.join(os.path.dirname(__file__), "../utils/"))
from aspmine_imports import *

'''------------------------------------------------------------------
# Code uses local functions: 
------------------------------------------------------------------

CustomArgumentParser
DBconnect
executeQuery

------------------------------------------------------------------'''

'''------------------------------------------------------------------
# Description:
------------------------------------------------------------------

Get all antismash annotations for PKS, NRPS and so on
Get all backbone homologs in other organisms
Create all against all matrix of genes
Cluster genes 

------------------------------------------------------------------'''


'''------------------------------------------------------------------
# ARGUMENTS and setup
------------------------------------------------------------------'''
startTime = datetime.now() # record runtime
parser = CustomArgumentParser(formatter_class=SmartFormatter, usage='%(prog)s -dbname [database name] -back [backbone type]')
parser.add_argument("-dbname", "-d", required=False, default = "aspminedb", help="Database name")
parser.add_argument("-back", "-b", required=True, default = "aspminedb", choices=["PKS","NRPS","all"], help="Backbone type")

""" TABLES """
parser.add_argument("-clean1", "-c1", required=False, action='store_true', help="Compile new t_antismashBackbones")

""" ANALYSIS """
parser.add_argument("-ana1", "-a1", required=False, action='store_true', help="Cluster VS organism - not best candidate, looping")
#parser.add_argument("-csize", "-cs", required=False, default = "5", help="Cluster size cutoff - default 5")

""" PARSE ARGUMENTS """
args 	= parser.parse_args()
clean1 = args.clean1
ana1 = args.ana1
backbone = args.back
startTime = datetime.now() # record runtime

'''------------------------------------------------------------------
# Print argument values to screen
------------------------------------------------------------------'''
print "#--------------------------------------------------------------\n\
# ARGUMENTS\n\
# Database\t: %s\n\
# Compile new t_antismashBackbones -clean1\t: %s\n\
# Analysis, cluster VS organism - with best candidate, looping\t\t: %s\n\
#--------------------------------------------------------------" % (args.dbname, clean1, ana1)

'''------------------------------------------------------------------
# Connection to database
------------------------------------------------------------------'''
cursor = DBconnect("localhost", "asp", "1234", args.dbname) # costum function

"""----------------------------------------------------------------------------
	TABLES
----------------------------------------------------------------------------"""

'''------------------------------------------------------------------
# Create connection table: t_antismashBackbones
------------------------------------------------------------------'''
print "# INFO: Processing - t_antismashBackbones"

# CREATE table if specified
if clean1:
	startTimet_antismashBackbones = datetime.now()

	print "# INFO: Connection antismash with t_antismash2blast - t_antismashBackbones"
	print "# INFO: Dropping and re-creating t_antismashBackbones"

	cursor.execute("DROP TABLE IF EXISTS t_antismashBackbones ")

	if backbone=="all":
		query = "CREATE TABLE t_antismashBackbones as \
		select t.* from antismash as a \
		join t_antismash2blast as t on (a.org_id=t.org_id and a.sm_protein_id=t.sm_protein_id);"
	else:
		query = "CREATE TABLE t_antismashBackbones as \
		select t.* from antismash as a \
		join t_antismash2blast as t on (a.org_id=t.org_id and a.sm_protein_id=t.sm_protein_id) where sm_short=\'" + backbone + "\';"

	executeQuery(cursor, query)

	executeQuery(cursor, "CREATE INDEX `i_sm_protein_id_org_id` ON `t_antismashBackbones` (`sm_protein_id`,`org_id`)")
	executeQuery(cursor, "CREATE INDEX `i_q_seqkey` ON `t_antismashBackbones` (q_seqkey)")
	executeQuery(cursor, "CREATE INDEX `i_h_seqkey` ON `t_antismashBackbones` (h_seqkey)")
	print "# INFO Runtime t_antismashBackbones: ", (datetime.now()-startTimet_antismashBackbones)

# If not specified to create table test if the table actually does exist
if not	cursor.execute("Show tables LIKE 't_antismashBackbones'"):
	sys.exit("# ERROR: table does not exist, re-run with -clean1 option")

"""----------------------------------------------------------------------------
	ANALYSIS
----------------------------------------------------------------------------"""
"""--------------------------------------
ORGANISM VS CLUSTER
--------------------------------------"""
if ana1:
	# If not specified to create table test if the table actually does exist
	if not cursor.execute("Show tables LIKE 't_antismashBackbones'"):
		sys.exit("# ERROR: table does not exist, re-run with -clean1 option")

	print "# INFO: selecting from t_antismashBackbones"

	startTimet_ana1 = datetime.now()

	query = "SELECT a.q_org, a.h_org, a.q_seqkey, a.h_seqkey, a.pident, a.q_len, a.h_len, a.q_cov, a.h_cov \
	from t_antismashBackbones as a join t_antismashBackbones as b on ( a.q_seqkey=b.h_seqkey) \
	where a.org_id=b.org_id and a.pident > 40 and a.q_seqkey>b.q_seqkey;"
	(columns, result) = executeQuery(cursor, query)

	print "# INFO Runtime t_antismashBackbones ana1: ", (datetime.now()-startTimet_ana1)

	cursor2csv(columns, result, "t_antismashBackbones.csv")

	print "# INFO: wrote results to t_antismashBackbones.csv"

	#print "# INFO: running Rscript"
	#os.system("R CMD BATCH '--args t_antismashLoopAntismash.csv' blast_antismash_subset.R test.out ")

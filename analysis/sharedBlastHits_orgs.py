#!/usr/bin/python

'''------------------------------------------------------------------
# Imports
------------------------------------------------------------------'''
import sys,os
root = os.path.abspath(__file__).split("github")[0]
sys.path.append(root+"github/utils/")
from aspmine_imports import *
import pandas as pd
from pandas import Series, DataFrame
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
parser = CustomArgumentParser(formatter_class=SmartFormatter, 
	usage='%(prog)s -pcut [Cutoff for percent id] -lowcut [Cutoff for lowest coverage id] -highcut [Cutoff for highest coverage id]')
parser.add_argument("-dbname", "-d", required=False, default = "aspminedb", help="Database name")
parser.add_argument("-pcut", "-pc", required=False, default = "50", help="Cutoff for percent id")
parser.add_argument("-lowcut", "-lc", required=False, default = "70", help="Cutoff for lowest coverage id")
parser.add_argument("-highcut", "-hc", required=False, default = "80", help="Cutoff for highest coverage id")

""" TABLES """
parser.add_argument("-clean1", "-c1", required=False, action='store_true', help="Create protein count table from table proteins")
parser.add_argument("-clean2", "-c2", required=False, action='store_true', help="Create all counts table from table protein count and blast")
parser.add_argument("-clean3", "-c3", required=False, action='store_true', help="Create cut counts table from table protein count and blast")
parser.add_argument("-force", "-f", required=False, action='store_true', help="Force to create all counts table from table protein count and blast")

""" ANALYSIS """
parser.add_argument("-allCounts", "-ac", required=False, action='store_true', help="Get numbers for all hits")
parser.add_argument("-cutCounts", "-cc", required=False, action='store_true', help="Get numbers for hits that fullfill the cutoffs")

""" PARSE ARGUMENTS """
args 	= parser.parse_args()
dbname 	= args.dbname
pcut	= args.pcut
lowcut	= args.lowcut
highcut	= args.highcut

clean1 = args.clean1
clean2 = args.clean2
clean3 = args.clean3

allCounts	= args.allCounts
cutCounts	= args.cutCounts

startTime = datetime.now() # record runtime

'''------------------------------------------------------------------
# Print argument values to screen
------------------------------------------------------------------'''
print "#--------------------------------------------------------------\n\
# ARGUMENTS\n\
# Database\t: %s\n\
# Cutoff for percent id\t\t\t: %s\n\
# Cutoff for lowest coverage id\t\t: %s\n\
# Cutoff for highest coverage id\t: %s\n\
# Analysis, get numbers for all hits (-allCounts)\t\t\t: %s\n\
# Analysis, get numbers for hits that fullfill the cutoffs (-cutCounts)\t: %s\n\
# Recreate t_proteinCounts -clean1\t: %s\n\
# Recreate t_allCounts -clean2 only with -force\t: %s\n\
#--------------------------------------------------------------" % (dbname, pcut, lowcut, highcut, allCounts, cutCounts, clean1, clean2)

'''------------------------------------------------------------------
# Connection to database
------------------------------------------------------------------'''
cursor = DBconnect("localhost", "asp", "1234", dbname) # costum function

'''------------------------------------------------------------------
# Get organisms from BLAST table
------------------------------------------------------------------'''
(columns, orgpairs)	= executeQuery(cursor, "SELECT distinct q_org, h_org from t_testblast")
(columns, q_orgs)	= executeQuery(cursor, "SELECT distinct q_org from t_testblast")
(columns, h_orgs)	= executeQuery(cursor, "SELECT distinct h_org from t_testblast")

print "# INFO: BLAST table contains %s unique organism pairs, %s q_orgs and %s h_orgs" % ( len(orgpairs), len(q_orgs), len(h_orgs) )
#print "# INFO: Query organisms:\t", sorted(list(q_orgs))
#print "# INFO: Hit organisms:\t\t", sorted(list(h_orgs))
		


'''------------------------------------------------------------------
# Create allCounts table

python sharedBlastHits_orgs.py -c2
#--------------------------------------------------------------
# ARGUMENTS
# Database	: aspminedb
# Cutoff for percent id			: 50
# Cutoff for lowest coverage id		: 70
# Cutoff for highest coverage id	: 80
# Analysis, get numbers for all hits (-allCounts)			: False
# Analysis, get numbers for hits that fullfill the cutoffs (-cutCounts)	: False
# Recreate t_proteinCounts	: False
# Recreate t_allCounts	: True
#--------------------------------------------------------------
# INFO: BLAST table contains 0 unique organism pairs, 0 q_orgs and 0 h_orgs
# INFO: Dropping and re-creating t_allCounts
# INFO: Creating table t_allCounts from blast
# INFO Runtime allCounts:  2:19:12.269561

------------------------------------------------------------------'''
if clean1:
	startTime_clean = datetime.now()

	print "# INFO: Dropping and re-creating t_proteinCounts"

	cursor.execute("DROP TABLE IF EXISTS t_proteinCounts")

	print "# INFO: Creating table t_proteinCounts from blast"
	query = "CREATE TABLE t_proteinCounts as SELECT org_id, prot_orgkey, count(*) as pcount from proteins group by prot_orgkey;"
	executeQuery( cursor, query )
	print "# INFO Runtime allCounts: ", (datetime.now()-startTime_clean)	

'''------------------------------------------------------------------
# Create allCounts table
------------------------------------------------------------------'''
if clean2 and force:
	startTime_clean = datetime.now()

	print "# INFO: Dropping and re-creating t_allCounts"

	cursor.execute("DROP TABLE IF EXISTS t_allCounts")
	#executeQuery(cursor, "CREATE INDEX qh_i ON blast (q_org, h_org)")

	print "# INFO: Creating table t_allCounts from blast"
	query = "CREATE TABLE t_allCounts as SELECT q_org, h_org, count(distinct(q_seqkey)) as count FROM blast GROUP BY q_org, h_org;"
	executeQuery( cursor, query )
	print "# INFO Runtime t_allCounts: ", (datetime.now()-startTime_clean)	

'''------------------------------------------------------------------
# Create cutCounts table
CREATE TABLE t_cutCounts as SELECT q_org, h_org, count(distinct(q_seqkey)) as count FROM blast\
WHERE pident >= 50 and GREATEST(q_cov, h_cov) >= 80 and LEAST(q_cov, h_cov) >= 50\
GROUP BY q_org, h_org;
Query OK, 2549 rows affected (18 min 54.82 sec)

------------------------------------------------------------------'''
if clean3 and force:
	startTime_clean = datetime.now()

	print "# INFO: Dropping and re-creating t_cutCounts"

	cursor.execute("DROP TABLE IF EXISTS t_cutCounts")

	print "# INFO: Creating table t_cutCounts from blast"
	query = "CREATE TABLE t_cutCounts as SELECT q_org, h_org, count(distinct(q_seqkey)) as count FROM blast\
	WHERE pident >= %s and GREATEST(q_cov, h_cov) >= %s and LEAST(q_cov, h_cov) >= %s\
	GROUP BY q_org, h_org;" % (pcut, highcut, lowcut)
	executeQuery( cursor, query )
	print "# INFO Runtime t_cutCounts: ", (datetime.now()-startTime_clean)	


'''------------------------------------------------------------------
# SELECT from allCounts - all
------------------------------------------------------------------'''
if allCounts:
	print "# INFO: SELECT all hits from t_allCounts"
	st = datetime.now()
	if not cursor.execute("Show tables LIKE 't_allCounts'"):
		sys.exit("# ERROR: table does not exist, re-run with -clean option")

	query = "SELECT t_allCounts.*, pa.pcount as qtotal, pb.pcount as htotal from t_allCounts\
	join t_proteinCounts pa on q_org=pa.prot_orgkey\
	join t_proteinCounts pb on h_org=pb.prot_orgkey"
	(columns, allCounts_nr)	= executeQuery( cursor, query )
	#print allCounts_nr

	print "# INFO: constructing matrix"

	t_q_orgs	= list(x[0] for x in allCounts_nr)
	t_h_orgs	= list(x[1] for x in allCounts_nr)
	t_shared	= list(x[2] for x in allCounts_nr)
	t_shared	= ["0" if v is None else v for v in t_shared]

	t_q_total 	= list(x[3] for x in allCounts_nr)
	t_q_total	= ["0" if v is None else v for v in t_q_total]

	t_h_total	= list(x[4] for x in allCounts_nr)
	t_h_total	= ["0" if v is None else v for v in t_h_total]
	
	uniq_q_orgs = list(set(t_q_orgs))
	uniq_h_orgs = list(set(t_h_orgs))
	uniq_orgs 	= sorted(list(set (uniq_q_orgs + uniq_h_orgs)))
	
	hit_matrix = pd.DataFrame(index=uniq_orgs, columns=uniq_orgs)
	print uniq_orgs
	for org1 in uniq_orgs:
		#i  = t_h_orgs.index(org1)
		org2 = t_q_orgs[ t_h_orgs.index(org1) ]

		hit_matrix.iloc[uniq_orgs.index(org2), uniq_orgs.index(org1)] = 1
		hit_matrix.iloc[uniq_orgs.index(org1), uniq_orgs.index(org2)] = 1

		#hit_matrix.iloc[uniq_orgs.index(org2), uniq_orgs.index(org1)] = (float(t_shared[i])/float(t_q_total[i]))*100
		#hit_matrix.iloc[uniq_orgs.index(org1), uniq_orgs.index(org2)] = (float(t_shared[i])/float(t_h_total[i]))*100

		#print org, t_q_orgs.index(org)
		# if org in t_q_orgs:
		# 	i  = t_q_orgs.index(org)
		# 	o2 = t_h_orgs[i]
		# 	hit_matrix.iloc[uniq_orgs.index(o2), uniq_orgs.index(org)] = 0
		# 	hit_matrix.iloc[uniq_orgs.index(org), uniq_orgs.index(o2)] = 0


	#print hit_matrix
  	hit_matrix.to_csv('hit_matrix.csv', sep=';',  encoding='utf-8')
	"""
	hit_matrix = pd.DataFrame(index=uniq_q_orgs, columns=uniq_h_orgs)
	for org in t_q_orgs:
		i = t_q_orgs.index(org)
		o1 = t_q_orgs[i]
		o2 = t_h_orgs[i]
		#print uniq_h_orgs.index(o2),uniq_q_orgs.index(o1)
		hit_matrix.iloc[uniq_q_orgs.index(o1),uniq_h_orgs.index(o2)] = (float(t_shared[i])/float(t_h_total[i]))*100
		#hit_matrix.iloc[uniq_q_orgs.index(o2),uniq_h_orgs.index(o1)] = t_shared[i]/t_q_total[i]*100
  	print hit_matrix			
  	"""
	#print t_q_orgs, t_h_orgs, t_shared, t_h_total, t_q_total
	#orgnames = list( set( t_q_orgs + t_h_orgs )	)
	#hit_matrix = pd.DataFrame(index=orgnames, columns=orgnames)
	#print hit_matrix
	#print hit_matrix	
	# Add 1's to the right cells
	#for org1 in orgnames:
		#for org2 in orgnames:
			#if org1 in t_q_orgs and org2 in t_h_orgs:
				#hit_matrix.iloc[t_q_orgs.index(org1),t_h_orgs.index(org2)] = 1
			#else:
		#print row,orgnames.index(row[0])
		#print t_q_orgs.index(org),t_h_orgs.index(org)
  		#	hit_matrix.iloc[orgnames.index(org),orgnames.index(org)] = 0	
			
	# Check for the right dimension
	print "# INFO: Matrix dimensions: ",hit_matrix.shape

	cursor2csv(columns, allCounts_nr, "allCounts_nr.csv")
	print "# INFO: wrote results to allCounts_nr.csv"
	print "# INFO: Runtime analysis of t_testblast: ", (datetime.now()-st)

	print "# INFO: running Rscript"
	#os.system("R CMD BATCH '--args allCounts_nr.csv' blast_antismash_subset.R test.out ")


'''------------------------------------------------------------------
# SELECT from allCounts - cutoffs
!!!THIS DOES OT WORK!!!
It will require to re-run the all counts table every time. 
------------------------------------------------------------------'''
"""
if cutCounts:
	print "# INFO: SELECT cutoff restricted hits from allCounts"

	st = datetime.now()
	if not cursor.execute("Show tables LIKE 't_allCounts'"):
		sys.exit("# ERROR: table does not exist, re-run with -clean option")

	print "# INFO: selecting cutoff determined hits from allCounts"
	query = "SELECT * from allCounts WHERE pident >= %s and GREATEST(q_cov, h_cov) >= %s and LEAST(q_cov, h_cov) >= %s;" % (pcut, highcut, lowcut)
	(columns, cutCounts_nr)	= executeQuery(cursor, query )
	#cutCounts_nr = cutCounts_nr[:-1] # Remove sumup row as result of ROLLUP

	cursor2csv(columns, cutCounts_nr, "cutCounts_nr.csv")
	print "# INFO: wrote results to cutCounts_nr.csv"
	print "# INFO: Runtime analysis of t_testblast: ", (datetime.now()-st)

	print "# INFO: running Rscript"
	#os.system("R CMD BATCH '--args cutCounts_nr.csv' blast_antismash_subset.R test.out ")
"""
"""


SELECT ta.*, \
( SELECT count(*) from proteins where prot_orgkey=ta.q_org group by prot_orgkey) as q_count, \
( SELECT count(*) from proteins where prot_orgkey=ta.h_org group by prot_orgkey) as h_count from \
( SELECT q_org, h_org, count(*) as count FROM t_testblast GROUP BY q_org, h_org) as ta;

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
"""
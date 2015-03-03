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
# ARGUMENTS and setup
------------------------------------------------------------------'''
startTime = datetime.now() # record runtime
parser = CustomArgumentParser(formatter_class=SmartFormatter, usage='%(prog)s -dbname [database name]')
parser.add_argument("-dbname", "-d", required=False, default = "aspminedb", help="Database name")

""" TABLES """
parser.add_argument("-clean1", "-c1", required=False, action='store_true', help="Compile new antismash2organism")
parser.add_argument("-clean2", "-c2", required=False, action='store_true', help="Compile new antismash2blast")
parser.add_argument("-clean3", "-c3", required=False, action='store_true', help="Compile new antismash2blast_reduced")
parser.add_argument("-loop", required=False, action='store_true', help="Compile antismashLoopAntismash")

""" ANALYSIS """
parser.add_argument("-ana1", "-a1", required=False, action='store_true', help="Cluster VS organism - not best candidate, looping")
parser.add_argument("-ana2", "-a2", required=False, action='store_true', help="Cluster VS Cluster - not best candidate, looping")
parser.add_argument("-anaLoop", "-aloop", required=False, action='store_true', help="Cluster VS organism - with best candidate, looping")

parser.add_argument("-csize", "-cs", required=False, default = "5", help="Cluster size cutoff - default 5")

""" PARSE ARGUMENTS """
args 	= parser.parse_args()
clean1 = args.clean1
clean2 = args.clean2
clean3 = args.clean3
csize = args.csize
ana1 = args.ana1
ana2 = args.ana2
anaLoop = args.anaLoop
loop = args.loop
startTime = datetime.now() # record runtime

if clean1:
	(clean2, clean3, loop) = (True,True,True)
if clean2:
	(clean3, loop) = (True,True)
if clean3:
	loop=True

'''------------------------------------------------------------------
# Print argument values to screen
------------------------------------------------------------------'''
print "#--------------------------------------------------------------\n\
# ARGUMENTS\n\
# Database\t: %s\n\
# Compile new antismash2organism -clean1\t: %s\n\
# Compile new antismash2blast -clean2\t\t: %s\n\
# Compile new antismash2blast_reduced -clean3\t: %s\n\
# Cutoff on cluster size\t\t: %s\n\
# Analysis, cluster VS organism - not best candidate, looping\t\t: %s\n\
# Analysis, cluster VS cluster - not best candidate, looping\t\t: %s\n\
# Analysis, cluster VS organism - with best candidate, looping\t\t: %s\n\
#--------------------------------------------------------------" % (args.dbname, clean1, clean2, clean3, csize, ana1, ana2, anaLoop)

'''------------------------------------------------------------------
# Connection to database
------------------------------------------------------------------'''
cursor = DBconnect("localhost", "asp", "1234", args.dbname) # costum function


"""----------------------------------------------------------------------------
	TABLES
----------------------------------------------------------------------------"""

'''------------------------------------------------------------------
# Create connection table: antismash -> organism
------------------------------------------------------------------'''
print "# INFO: Processing - antismash2organism"

# CREATE table if specified
if clean1:
	startTimeantismash2organism = datetime.now()

	print "# INFO: Connection antismash with organism - antismash2organism"
	print "# INFO: Dropping and re-creating antismash2organism"

	cursor.execute("DROP TABLE IF EXISTS antismash2organism ")
	query = "CREATE TABLE antismash2organism as \
	select ta.org_id, sm_protein_id, clust_id, name \
	from ( select org_id, sm_protein_id, concat(org_id, '_', clust_backbone, '_', clust_size) as clust_id  from antismash) as ta \
	join organism using (org_id);"

	executeQuery(cursor, query)

	executeQuery(cursor, "CREATE INDEX `sm_protein_id` ON `antismash2organism` (`sm_protein_id`,`name`)")
	executeQuery(cursor, "CREATE INDEX `org_id`  ON `antismash2organism`  (`org_id`)")
	executeQuery(cursor, "CREATE INDEX `name`  ON `antismash2organism` (`name`)")
	executeQuery(cursor, "CREATE INDEX `sm_protein_id_2`  ON `antismash2organism` (`sm_protein_id`)")
	print "# INFO Runtime antismash2organism: ", (datetime.now()-startTimeantismash2organism)

# If not specified to create table test if the table actually does exist
if not	cursor.execute("Show tables LIKE 'antismash2organism'"):
	sys.exit("# ERROR: table does not exist, re-run with -clean1 option")

# SELECT data from table
query = "SELECT distinct antismash2organism.name, organism.name from antismash2organism right join organism using (name) where antismash2organism.name is null;"
(columns, lackingAnti) = executeQuery(cursor, query)

query = "SELECT distinct antismash2organism.name, organism.name from antismash2organism right join organism using (name) where organism.name is null;"
(columns, lackingOrg) = executeQuery(cursor, query)

if lackingAnti:
	orgs = [o[1] for o in lackingAnti]
	print ("# WARNING: There are orgs in organism not found in antismash %s" % orgs )

if lackingOrg:
	orgs = [o[1] for o in lackingOrg]
	sys.exit("# ERROR: There are orgs in antismash not found in organism, %s" % orgs)

'''------------------------------------------------------------------
# Create connection table: antismash-organism -> blast
------------------------------------------------------------------'''

print "# INFO: Processing - antismash2blast and antismash2blast_reduced"

# CREATE table if specified
if clean2:
	startTimeantismash2blast = datetime.now()
	print "# INFO: Connection antismash2organism with BLAST - antismash2blast, this will take a while"
	print "# INFO: Dropping and re-creating antismash2blast"

	cursor.execute("DROP TABLE IF EXISTS antismash2blast")
	query = "CREATE TABLE antismash2blast as select * from antismash2organism join blast on (sm_protein_id = q_seqkey and q_org = name) where h_cov+q_cov >= 180;"

	executeQuery(cursor, query)

	print "# INFO: Createing index"
	executeQuery(cursor, "CREATE INDEX `sm_protein_id` ON `antismash2blast` (`sm_protein_id`,`name`)")
	executeQuery(cursor, "CREATE INDEX `name` ON `antismash2blast` (`name`)")
	executeQuery(cursor, "CREATE INDEX `sm_protein_id_2`  ON `antismash2blast` (`sm_protein_id`)")
	executeQuery(cursor, "CREATE INDEX `seqids`  ON `antismash2blast` (`h_seqid`, `q_seqid`)")
	executeQuery(cursor, "CREATE INDEX `qseq_horg_pident`  ON `antismash2blast` (`q_seqid`, `h_org`, `pident`)")
	print "# INFO Runtime antismash2blast: ", (datetime.now()-startTimeantismash2blast)

if clean3:
	startTimeantismash2blast_reduced = datetime.now()
	print "# INFO: Reducing antismash2blast to best match for each organism-cluster pair - antismash2blast_reduced"
	print "# INFO: Dropping and re-creating antismash2blast_reduced"

	cursor.execute("DROP TABLE IF EXISTS antismash2blast_reduced")

	query="CREATE TABLE antismash2blast_reduced as \
		select * from (\
		select ta.* from antismash2blast as ta\
		join ( select q_seqid, h_seqid, h_org, max(pident) as max from antismash2blast group by q_seqid, h_org) as tb \
		on ta.q_seqid=tb.q_seqid and ta.h_org=tb.h_org and  pident=max) as tc group by q_seqid, h_org"
	executeQuery(cursor, query)

	print "# INFO Runtime antismash2blast_reduced: ", (datetime.now()-startTimeantismash2blast_reduced)
	

query = "SELECT count(distinct(sm_protein_id)) from antismash2blast_reduced;"
(columns, smGenesRed) = executeQuery(cursor, query)

query = "SELECT count(distinct(sm_protein_id)) from antismash2blast;"
(columns, smGenes) = executeQuery(cursor, query)

if smGenesRed and smGenes:
	if smGenesRed[0][0] != smGenes[0][0] :
		print "# WARNING: number of genes in antismash2blast_reduced and antismash2blast are not equal, reduced= %s, %s" % (smGenesRed[0][0], smGenes[0][0]) 
		print "# INFO: Executing diagnose query - this will take some time" 
		diagnose_query="SELECT abp from \
					(SELECT distinct(sm_protein_id) as abp from antismash2blast) as ta\
					left join \
					(SELECT distinct(sm_protein_id) as abrp from antismash2blast_reduced) as tb\
					on abp=abrp where (abrp is NULL);"
		(columns, abGenes) = executeQuery(cursor, diagnose_query)				
		missingGenes = [o for o in abGenes]
		#print missingGenes

	else:
		print "# INFO: number of antismash genes in antismash2blast_reduced and antismash2blast %s" % smGenesRed[0] 

"""--------------------------------------
antismashLoopAntismash - best candicate cluster table
--------------------------------------"""
if loop:
	
	print "# INFO: Looping antismash2blast reduced back to antismash - antismashLoopAntismash"
	print "# INFO: Dropping and re-creating antismashLoopAntismash"

	cursor.execute("DROP TABLE IF EXISTS antismashLoopAntismash")

	startTimeantismashLoopAntismash = datetime.now()
	query = "CREATE table antismashLoopAntismash as \
		SELECT antismash2blast_reduced.org_id as q_orgid, clust_id as q_clustid, SUBSTRING_INDEX(clust_id, '_', -1) as q_clust_size,\
		antismash2blast_reduced.name as q_orgname,  organism.org_id as h_orgid, h_org as h_orgname, h_seqkey, \
		clust_backbone as h_clust_backbone, O1.section as q_sec, O1.real_name as q_realname, O2.section as h_sec, O2.real_name as h_realname, \
		antismash.sm_protein_id as h_clust_protein_id, q_cov, h_cov, pident\
		from antismash2blast_reduced   \
		join organism on h_org=organism.name\
		left join antismash on organism.org_id=antismash.org_id and h_seqkey=antismash.sm_protein_id  \
		join organism O1 on (O1.name=antismash2blast_reduced.name) \
		join organism O2 on (O2.name=h_org)"

	executeQuery(cursor, query)

	print "# INFO: Temporary table, members per candidate cluster from antismashLoopAntismash - antismashLoopAntismashCandidates_tmp"
	print "# INFO: Dropping and re-creating antismashLoopAntismashCandidates_tmp"
	cursor.execute("DROP TABLE IF EXISTS antismashLoopAntismashCandidates_tmp")

	query = "CREATE table antismashLoopAntismashCandidates_tmp as \
			SELECT *, count(*) as candidateMembers  \
			from antismashLoopAntismash  \
			where h_clust_backbone is not null \
			group by q_clustid, h_clust_backbone"
	executeQuery(cursor, query)		
	executeQuery(cursor, "CREATE INDEX `index3` ON `antismashLoopAntismashCandidates_tmp` (`q_clustid`,`h_orgname`,candidateMembers)")

	print "# INFO: Getting best candidate cluster from antismashLoopAntismashCandidates_tmp - antismashLoopAntismashCandidates"
	print "# INFO: Dropping and re-creating antismashLoopAntismashCandidates"
	cursor.execute("DROP TABLE IF EXISTS antismashLoopAntismashCandidates")

	query = "CREATE TABLE antismashLoopAntismashCandidates as \
			SELECT tb.q_orgname, tb.q_clust_size, tb.q_clustid, tb.h_orgname, tb.candidateMembers, tb.h_realname, tb.q_realname, \
			tb.candidateMembers/tb.q_clust_size as clustCov, tc.h_clust_backbone, tc.q_sec, tc.h_sec, tb.q_orgid, tb.h_orgid \
			from ( select ta.*, max(ta.candidateMembers) as max_candidateMembers \
			from  antismashLoopAntismashCandidates_tmp as ta \
			group by ta.q_clustid, h_orgname ) as tb \
			join antismashLoopAntismashCandidates_tmp as tc \
			on tc.q_clustid=tb.q_clustid and tc.h_orgname=tb.h_orgname and tc.candidateMembers=max_candidateMembers ; "

	executeQuery(cursor, query)			

	print "# INFO Runtime antismashLoopAntismash: ", (datetime.now()-startTimeantismashLoopAntismash)


"""----------------------------------------------------------------------------
	ANALYSIS
----------------------------------------------------------------------------"""

"""--------------------------------------
antismashLoopAntismash - best candicate cluster analysis
--------------------------------------"""
if anaLoop:

	startTimeanaLoop = datetime.now()

	if not cursor.execute("Show tables LIKE 'antismashLoopAntismash'"):
		sys.exit("# ERROR: table does not exist, re-run with -loop option")

	print "# INFO: selecting from antismashLoopAntismash"
	query = "SELECT * from antismashLoopAntismashCandidates"

	(columns, result) = executeQuery(cursor, query)
	cursor2csv(columns, result, "antismashLoopAntismash.csv")

	print "# INFO: wrote results to antismashLoopAntismash.csv"
	print "# INFO Runtime analysis of antismashLoopAntismash: ", (datetime.now()-startTimeanaLoop)

	print "# INFO: running Rscript"
	os.system("R CMD BATCH '--args antismashLoopAntismash.csv' blast_antismash_subset.R test.out ")


"""--------------------------------------
ORGANISM VS CLUSTER
--------------------------------------"""
if ana1:
	# If not specified to create table test if the table actually does exist
	if not cursor.execute("Show tables LIKE 'antismash2blast'") or not cursor.execute("Show tables LIKE 'antismash2blast_reduced'"):
		sys.exit("# ERROR: table does not exist, re-run with -clean2 option")

	print "# INFO: selecting from antismash2blast_reduced"

	query = "SELECT ta.*, ROUND(nrHits/clust_size, 2) as clustCov, O1.section as sec1, O2.section as sec2  \
			from (select org_id, clust_id, name, h_org, count(*) as nrHits,\
			SUBSTRING_INDEX(clust_id, '_', -1) as clust_size from antismash2blast_reduced group by clust_id, h_org ) as ta \
			join organism O1 using (name)\
			join organism O2 on (O2.name=h_org)\
			where clust_size >= " + str(csize)+";"
	(columns, result) = executeQuery(cursor, query)
	cursor2csv(columns, result, "hitsPerOrgPerClusterGene.csv")

	print "# INFO: wrote results to hitsPerOrgPerClusterGene.csv"


"""--------------------------------------
CLUSTER VS CLUSTER
--------------------------------------"""

if ana2:
	print "# INFO: Dropping and re-creating clusterVScluster"
	cursor.execute("DROP TABLE IF EXISTS clusterVScluster")
	print "# INFO: selecting from clusterVScluster"

	query="CREATE TABLE clusterVScluster as\
	SELECT ta.clust_id ta_clust_id , ta.q_org ta_q_org , ta.q_seqid ta_q_seqid, ta.h_org ta_h_org, ta.h_seqid ta_h_seqid, \
	ta.pident ta_pident, ta.q_cov ta_q_cov, ta.h_cov ta_h_cov, SUBSTRING_INDEX(ta.clust_id, '_', -1) as ta_clust_size, \
	tb.clust_id tb_clust_id , tb.q_org tb_q_org , tb.q_seqid tb_q_seqid, tb.h_org tb_h_org, tb.h_seqid tb_h_seqid, \
	tb.pident tb_pident, tb.q_cov tb_q_cov, tb.h_cov tb_h_cov, SUBSTRING_INDEX(tb.clust_id, '_', -1) as tb_clust_size \
	from antismash2blast as ta \
	LEFT join antismash2blast as tb \
	on ta.q_seqid=tb.h_seqid and ta.h_seqid=tb.q_seqid;"

	query="SELECT *, ROUND(over/ta_clust_size, 2) as ta_clustCov, ROUND(over/tb_clust_size, 2) as tb_clustCov \
	from ( SELECT *, count(*) as over from clusterVScluster group by ta_clust_id, tb_clust_id) as ta;"

	(columns, result) = executeQuery(cursor, query)
	cursor2csv(columns, result, "hitsPerClusterPerClusterGene.csv")






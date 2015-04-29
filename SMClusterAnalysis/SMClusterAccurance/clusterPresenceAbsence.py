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
parser = CustomArgumentParser(formatter_class=SmartFormatter, usage='%(prog)s -dbname [database name]')
parser.add_argument("-dbname", "-d", required=False, default = "aspminedb", help="Database name")
parser.add_argument("-Rscript", "-R", required=False, default = "/home/tcve/github/tcve/blast_antismash_subset.R", help="Rscript name/path")
parser.add_argument("-specificOrgs", "-sp", nargs = '*',  required=False, default=[None], action='store', help="List of specificOrgs for plots")
parser.add_argument("-section", "-sec", required=False, default = "*", help="Aspergillus section")
parser.add_argument("-afolder", "-af",  required=False, default="afolder_clusterPresenceAbsence", help="Name of folder for analysis results")

""" TABLES """
parser.add_argument("-clean1", "-c1", required=False, action='store_true', help="Compile new t_antismash2organism")
parser.add_argument("-clean2", "-c2", required=False, action='store_true', help="Compile new t_antismash2blast")
parser.add_argument("-clean3", "-c3", required=False, action='store_true', help="Compile new t_antismash2blast_reduced")
parser.add_argument("-loop", required=False, action='store_true', help="Compile t_antismashLoopAntismash")

""" ANALYSIS """
parser.add_argument("-a1", required=False, action='store_true', help="Cluster VS organism - not best candidate, looping")
parser.add_argument("-a2", required=False, action='store_true', help="Cluster VS Cluster - not best candidate, looping")
parser.add_argument("-a3", required=False, action='store_true', help="Cluster VS organism - with best candidate, looping")

parser.add_argument("-csize", "-cs", required=False, default = "5", help="Cluster size cutoff - default 5")

""" PARSE ARGUMENTS """
args 	= parser.parse_args()
Rscript = args.Rscript
clean1 = args.clean1
clean2 = args.clean2
clean3 = args.clean3
csize = args.csize
a1 = args.a1
a2 = args.a2
a3 = args.a3
loop = args.loop
specificOrgs = args.specificOrgs 
section = args.section

#sys.exit(afolder)

startTime = datetime.now() # record runtime

if clean1:
	(clean2, clean3, loop) = (True,True,True)
if clean2:
	(clean3, loop) = (True,True)
if clean3:
	loop=True

root = os.path.abspath(__file__).split("/")[1:-1]
Rpath = "/" + "/".join(root) + "/"
afolder = os.getcwd() + "/" + args.afolder
os.system("mkdir -p " + afolder + "" )
afolderString = afolder.replace("/", "\\/")
	
'''------------------------------------------------------------------
# Print argument values to screen
------------------------------------------------------------------'''
print "#--------------------------------------------------------------\n\
# ARGUMENTS\n\
# Database\t: %s\n\
# Name of Rscript -Rscript\t: %s\n\
# List of specificOrgs -sp\t: %s\n\
# Compile new t_antismash2organism -clean1\t: %s\n\
# Compile new t_antismash2blast -clean2\t\t: %s\n\
# Compile new t_antismash2blast_reduced -clean3\t: %s\n\
# Cutoff on cluster size\t\t: %s\n\
# Analysis, cluster VS organism - not best candidate, looping\t\t: %s\n\
# Analysis, cluster VS cluster - not best candidate, looping\t\t: %s\n\
# Analysis, cluster VS organism - with best candidate, looping\t\t: %s\n\
#--------------------------------------------------------------" % (args.dbname, Rscript, specificOrgs, clean1, clean2, clean3, csize, a1, a2, a3)

if specificOrgs != [None] and (not a3 or not a1):
	print "# WARNING: specificOrgs is only applicable to analysis 3 and 1, -a3/-a1"


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
print "# INFO: Processing - t_antismash2organism"

# CREATE table if specified
if clean1:
	startTimet_antismash2organism = datetime.now()

	print "# INFO: Connection antismash with organism - t_antismash2organism"
	print "# INFO: Dropping and re-creating t_antismash2organism"

	cursor.execute("DROP TABLE IF EXISTS t_antismash2organism ")
	query = "CREATE TABLE t_antismash2organism as \
	select ta.org_id, sm_protein_id, clust_id, name \
	from ( select org_id, sm_protein_id, concat(org_id, '_', clust_backbone, '_', clust_size) as clust_id  from antismash) as ta \
	join organism using (org_id);"

	executeQuery(cursor, query)

	executeQuery(cursor, "CREATE INDEX `sm_protein_id` ON `t_antismash2organism` (`sm_protein_id`,`name`)")
	executeQuery(cursor, "CREATE INDEX `org_id`  ON `t_antismash2organism`  (`org_id`)")
	executeQuery(cursor, "CREATE INDEX `name`  ON `t_antismash2organism` (`name`)")
	executeQuery(cursor, "CREATE INDEX `sm_protein_id_2`  ON `t_antismash2organism` (`sm_protein_id`)")
	print "# INFO Runtime t_antismash2organism: ", (datetime.now()-startTimet_antismash2organism)

# If not specified to create table test if the table actually does exist
if not	cursor.execute("Show tables LIKE 't_antismash2organism'"):
	sys.exit("# ERROR: table does not exist, re-run with -clean1 option")

# SELECT data from table
query = "SELECT distinct t_antismash2organism.name, organism.name from t_antismash2organism right join organism using (name) where t_antismash2organism.name is null;"
(columns, lackingAnti) = executeQuery(cursor, query)

query = "SELECT distinct t_antismash2organism.name, organism.name from t_antismash2organism right join organism using (name) where organism.name is null;"
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

print "# INFO: Processing - t_antismash2blast and t_antismash2blast_reduced"

# CREATE table if specified
if clean2:
	startTimet_antismash2blast = datetime.now()
	print "# INFO: Connection t_antismash2organism with BLAST - t_antismash2blast, this will take a while"
	print "# INFO: Dropping and re-creating t_antismash2blast"

	cursor.execute("DROP TABLE IF EXISTS t_antismash2blast")
	query = "CREATE TABLE t_antismash2blast as select * from t_antismash2organism join blast on \
	(sm_protein_id = q_seqkey and q_org = name) where h_cov+q_cov >= 160 and LEAST(h_cov, q_cov) >= 70;"

	executeQuery(cursor, query)

	print "# INFO: Createing index"
	executeQuery(cursor, "CREATE INDEX `sm_protein_id` ON `t_antismash2blast` (`sm_protein_id`,`name`)")
	executeQuery(cursor, "CREATE INDEX `name` ON `t_antismash2blast` (`name`)")
	executeQuery(cursor, "CREATE INDEX `sm_protein_id_2`  ON `t_antismash2blast` (`sm_protein_id`)")
	executeQuery(cursor, "CREATE INDEX `seqids`  ON `t_antismash2blast` (`h_seqid`, `q_seqid`)")
	executeQuery(cursor, "CREATE INDEX `qseq_horg_pident`  ON `t_antismash2blast` (`q_seqid`, `h_org`, `pident`)")
	print "# INFO Runtime t_antismash2blast: ", (datetime.now()-startTimet_antismash2blast)

if clean3:
	startTimet_antismash2blast_reduced = datetime.now()
	print "# INFO: Reducing t_antismash2blast to best match for each organism-cluster pair - t_antismash2blast_reduced"
	print "# INFO: Dropping and re-creating t_antismash2blast_reduced"

	cursor.execute("DROP TABLE IF EXISTS t_antismash2blast_reduced")

	query="CREATE TABLE t_antismash2blast_reduced as \
		select * from (\
		select ta.* from t_antismash2blast as ta\
		join ( select q_seqid, h_seqid, h_org, max(pident) as max from t_antismash2blast group by q_seqid, h_org) as tb \
		on ta.q_seqid=tb.q_seqid and ta.h_org=tb.h_org and  pident=tb.max) as tc group by q_seqid, h_org;"
	executeQuery(cursor, query)

	print "# INFO Runtime t_antismash2blast_reduced: ", (datetime.now()-startTimet_antismash2blast_reduced)
	

query = "SELECT count(distinct(sm_protein_id)) from t_antismash2blast_reduced;"
(columns, smGenesRed) = executeQuery(cursor, query)

query = "SELECT count(distinct(sm_protein_id)) from t_antismash2blast;"
(columns, smGenes) = executeQuery(cursor, query)

if smGenesRed and smGenes:
	if smGenesRed[0][0] != smGenes[0][0] :
		print "# WARNING: number of genes in t_antismash2blast_reduced and t_antismash2blast are not equal, reduced= %s, %s" % (smGenesRed[0][0], smGenes[0][0]) 
		print "# INFO: Executing diagnose query - this will take some time" 
		diagnose_query="SELECT abp from \
					(SELECT distinct(sm_protein_id) as abp from t_antismash2blast) as ta\
					left join \
					(SELECT distinct(sm_protein_id) as abrp from t_antismash2blast_reduced) as tb\
					on abp=abrp where (abrp is NULL);"
		(columns, abGenes) = executeQuery(cursor, diagnose_query)				
		missingGenes = [o for o in abGenes]
		#print missingGenes

	else:
		print "# INFO: number of antismash genes in t_antismash2blast_reduced and t_antismash2blast %s" % smGenesRed[0] 

"""--------------------------------------
t_antismashLoopAntismash - best candicate cluster table
--------------------------------------"""
if loop:
	
	print "# INFO: Looping t_antismash2blast reduced back to antismash - t_antismashLoopAntismash"
	print "# INFO: Dropping and re-creating t_antismashLoopAntismash"

	cursor.execute("DROP TABLE IF EXISTS t_antismashLoopAntismash")

	startTimet_antismashLoopAntismash = datetime.now()
	query = "CREATE table t_antismashLoopAntismash as \
		SELECT t_antismash2blast_reduced.org_id as q_orgid, clust_id as q_clustid, SUBSTRING_INDEX(clust_id, '_', -1) as q_clust_size,\
		t_antismash2blast_reduced.name as q_orgname,  organism.org_id as h_orgid, h_org as h_orgname, h_seqkey, \
		clust_backbone as h_clust_backbone, O1.section as q_sec, O1.real_name as q_realname, O2.section as h_sec, O2.real_name as h_realname, \
		antismash.sm_protein_id as h_clust_protein_id, q_cov, h_cov, pident\
		from t_antismash2blast_reduced   \
		join organism on h_org=organism.name\
		left join antismash on organism.org_id=antismash.org_id and h_seqkey=antismash.sm_protein_id  \
		join organism O1 on (O1.name=t_antismash2blast_reduced.name) \
		join organism O2 on (O2.name=h_org)"

	executeQuery(cursor, query)

	print "# INFO: Temporary table, members per candidate cluster from t_antismashLoopAntismash - t_antismashLoopAntismashCandidates_tmp"
	print "# INFO: Dropping and re-creating t_antismashLoopAntismashCandidates_tmp"
	cursor.execute("DROP TABLE IF EXISTS t_antismashLoopAntismashCandidates_tmp")

	query = "CREATE table t_antismashLoopAntismashCandidates_tmp as \
			SELECT *, count(*) as candidateMembers  \
			from t_antismashLoopAntismash  \
			where h_clust_backbone is not null \
			group by q_clustid, h_clust_backbone"
	executeQuery(cursor, query)		
	executeQuery(cursor, "CREATE INDEX `index3` ON `t_antismashLoopAntismashCandidates_tmp` (`q_clustid`,`h_orgname`,candidateMembers)")

	print "# INFO: Getting best candidate cluster from t_antismashLoopAntismashCandidates_tmp - t_antismashLoopAntismashCandidates"
	print "# INFO: Dropping and re-creating t_antismashLoopAntismashCandidates"
	cursor.execute("DROP TABLE IF EXISTS t_antismashLoopAntismashCandidates")

	query = "CREATE TABLE t_antismashLoopAntismashCandidates as \
			SELECT tb.q_orgname, tb.q_clust_size, tb.q_clustid, tb.h_orgname, tb.candidateMembers, tb.h_realname, tb.q_realname, \
			tc.candidateMembers/tc.q_clust_size as clustCov, tc.h_clust_backbone, tc.q_sec, tc.h_sec, tb.q_orgid, tb.h_orgid \
			from ( select ta.*, max(ta.candidateMembers) as max_candidateMembers \
			from  t_antismashLoopAntismashCandidates_tmp as ta \
			group by ta.q_clustid, h_orgname ) as tb \
			join t_antismashLoopAntismashCandidates_tmp as tc \
			on tc.q_clustid=tb.q_clustid and tc.h_orgname=tb.h_orgname and tc.candidateMembers=max_candidateMembers ; "

	executeQuery(cursor, query)			

	print "# INFO Runtime t_antismashLoopAntismash: ", (datetime.now()-startTimet_antismashLoopAntismash)


"""----------------------------------------------------------------------------
	ANALYSIS
----------------------------------------------------------------------------"""
def main():

	print "# INFO: creating analysis folder:" + afolder

	os.system("sed \"s/replace_folderpath/\'" + afolderString + "\'/g\" " + Rpath + "/clusterPresenceAbsence_generic.R > " + afolder + "/tmp.R")

	if a1:
		( data1, data2, data3, file_string  ) = analysis1( afolder, afolderString )
		os.system("sed -i \"s/replace_analysisNumber/\'a1\'/g\" " + afolder + "/tmp.R")
	if a2:
		( data1, data2, data3, file_string  ) = analysis2( afolder, afolderString )
		os.system("sed -i \"s/replace_analysisNumber/\'a2\'/g\" " + afolder + "/tmp.R")
	if a3:
		( data1, data2, data3, file_string  ) = analysis3( afolder, afolderString )
		os.system("sed -i \"s/replace_analysisNumber/\'a3\'/g\" " + afolder + "/tmp.R")

	if a1 or a2 or a3:
		os.system("sed -i \"s/replace_data1/" + data1 + "/g\" " + afolder + "/tmp.R")
		os.system("sed -i \"s/replace_data2/" + data2 + "/g\" " + afolder + "/tmp.R")
		os.system("sed -i \"s/replace_data3/" + data3 + "/g\" " + afolder + "/tmp.R")
		os.system("sed -i \"s/replace_infile/\'" + file_string + "\'/g\" " + afolder + "/tmp.R")
		os.system("R CMD BATCH " + afolder + "/tmp.R " + afolder + "/test.out ")

"""--------------------------------------
t_antismashLoopAntismash - best candicate cluster analysis
--------------------------------------"""
def analysis3( afolder, afolderString ):

	startTime = datetime.now()
	file_string = afolderString + "\\/t_clusterPresenceAbsence_analysis3.csv"

	if not cursor.execute("Show tables LIKE 't_antismashLoopAntismash'"):
		sys.exit("# ERROR: table does not exist, re-run with -loop option")

	print "# INFO: selecting from t_antismashLoopAntismash"
	query = "SELECT * from t_antismashLoopAntismashCandidates"

	(columns, result) = executeQuery(cursor, query)
	cursor2csv(columns, result, afolder + "/t_clusterPresenceAbsence_analysis3.csv")

	print "# INFO: wrote results to t_clusterPresenceAbsence_analysis3.csv"
	print "# INFO: Runtime analysis of t_antismashLoopAntismash: ", (datetime.now()-startTime)

	if specificOrgs == [None]:
		print "# INFO: creating section argument"
		secString = "grepl(\'" + section + "\', tmpdat\$q_sec) \& grepl(\'" + section + "\', tmpdat\$h_sec)"
		#secString = "tmpdat\$q_sec==\'" + section + "\' \& tmpdat\$h_sec==\'" + section + "\'"
		data1 = "data1 <- subset(tmpdat, (" + secString + "), select=c(q_orgid,q_clustid,h_realname,clustCov))"
		data2 = "data2 <- subset(tmpdat, tmpdat\$clust_size>10 \& (" + secString + "), select=c(q_orgid,q_clustid,h_realname,clustCov)) "
		data3 = "data3 <- subset(tmpdat, (" + secString + "), select=c(q_orgid,q_clustid,h_realname,clustCov,q_orgname)) "

	if specificOrgs != [None]:
		print "# INFO: creating organism arguments"
		count = 0
		orgString = ''
		for org in specificOrgs:
			orgString = orgString + " tmpdat\$q_orgname==\'" + str(org) + "\' "
			count += 1
			if not count == len(specificOrgs):
				orgString = orgString + " | "
		data1 = "data1 <- subset(tmpdat, (" + orgString + "), select=c(q_orgid,q_clustid,h_realname,clustCov))"
		data2 = "data2 <- subset(tmpdat, tmpdat\$clust_size>10 \& (" + orgString + "), select=c(q_orgid,q_clustid,h_realname,clustCov)) "
		data3 = "data3 <- subset(tmpdat, (" + orgString + "), select=c(q_orgid,q_clustid,h_realname,clustCov,q_orgname)) "

	return ( data1, data2, data3, file_string )			

"""--------------------------------------
ORGANISM VS CLUSTER
--------------------------------------"""
def analysis1( afolder, afolderString ):
	startTime = datetime.now()
	file_string = afolderString + "\\/t_clusterPresenceAbsence_analysis1.csv"


	#sys.exit()
	# If not specified to create table test if the table actually does exist
	if not cursor.execute("Show tables LIKE 't_antismash2blast'") or not cursor.execute("Show tables LIKE 't_antismash2blast_reduced'"):
		sys.exit("# ERROR: table does not exist, re-run with -clean2 option")

	print "# INFO: selecting from t_antismash2blast_reduced"

	query = "SELECT ta.*, ROUND(nrHits/clust_size, 2) as clustCov, O1.section as sec1, O2.section as sec2  \
			from (select org_id, clust_id, name, h_org, count(*) as nrHits,\
			SUBSTRING_INDEX(clust_id, '_', -1) as clust_size from t_antismash2blast_reduced group by clust_id, h_org ) as ta \
			join organism O1 using (name)\
			join organism O2 on (O2.name=h_org)\
			where clust_size >= " + str(csize)+";"

	(columns, result) = executeQuery(cursor, query)
	cursor2csv(columns, result, afolder + "/t_clusterPresenceAbsence_analysis1.csv")

	print "# INFO: Runtime analysis of t_antismashLoopAntismash: ", (datetime.now()-startTime)
	print "# INFO: wrote results to t_clusterPresenceAbsence_analysis1.csv"

	if specificOrgs == [None]:
		print "# INFO: creating section argument"
		secString = "grepl(\'" + section + "\', tmpdat\$sec1) \& grepl(\'" + section + "\', tmpdat\$sec2)"
		#secString = "tmpdat\$sec1==\'" + section + "\' \& tmpdat\$sec2==\'" + section + "\'"
		data1 = "data1 <- subset(tmpdat, (" + secString + "), select=c(q_orgid,q_clustid,h_realname,clustCov))"
		data2 = "data2 <- subset(tmpdat, tmpdat\$clust_size>10 \& (" + secString + "), select=c(q_orgid,q_clustid,h_realname,clustCov)) "
		data3 = "data3 <- subset(tmpdat, (" + secString + "), select=c(q_orgid,q_clustid,h_realname,clustCov,q_orgname)) "

	if specificOrgs != [None]:
		print "# INFO: creating organism arguments"
		count = 0
		orgString = ''
		for org in specificOrgs:
			orgString = orgString + " tmpdat\$q_orgname==\'" + str(org) + "\' "
			count += 1
			if not count == len(specificOrgs):
				orgString = orgString + " | "
		data1 = "data1 <- subset(tmpdat, (" + orgString + "), select=c(q_orgid,q_clustid,h_realname,clustCov))"
		data2 = "data2 <- subset(tmpdat, tmpdat\$clust_size>10 \& (" + orgString + "), select=c(q_orgid,q_clustid,h_realname,clustCov)) "
		data3 = "data3 <- subset(tmpdat, (" + orgString + "), select=c(q_orgid,q_clustid,h_realname,clustCov,q_orgname)) "
		
	return ( data1, data2, data3, file_string)			

"""--------------------------------------
CLUSTER VS CLUSTER
--------------------------------------"""

def analysis2( afolder, afolderString ):
	startTime = datetime.now()
	file_string = afolderString + "\\/t_clusterPresenceAbsence_analysis2.csv"

	print "# INFO: Dropping and re-creating t_clusterVScluster"
	cursor.execute("DROP TABLE IF EXISTS t_clusterVScluster")

	print "# INFO: selecting from t_clusterVScluster"

	query="CREATE TABLE t_clusterVScluster as\
	SELECT ta.clust_id ta_clust_id, ta.q_org ta_q_org, ta.q_seqid ta_q_seqid, ta.h_org ta_h_org, ta.h_seqid ta_h_seqid, \
	ta.pident ta_pident, ta.q_cov ta_q_cov, ta.h_cov ta_h_cov, SUBSTRING_INDEX(ta.clust_id, '_', -1) as ta_clust_size, \
	tb.clust_id tb_clust_id, tb.q_org tb_q_org, tb.q_seqid tb_q_seqid, tb.h_org tb_h_org, tb.h_seqid tb_h_seqid, \
	tb.pident tb_pident, tb.q_cov tb_q_cov, tb.h_cov tb_h_cov, SUBSTRING_INDEX(tb.clust_id, '_', -1) as tb_clust_size \
	from t_antismash2blast as ta \
	LEFT join t_antismash2blast as tb \
	on ta.q_seqid=tb.h_seqid and ta.h_seqid=tb.q_seqid;"
	(columns, result) = executeQuery(cursor, query)

	query="SELECT *, ROUND(over/ta_clust_size, 2) as ta_clustCov, ROUND(over/tb_clust_size, 2) as tb_clustCov \
	from ( SELECT *, count(*) as over from t_clusterVScluster group by ta_clust_id, tb_clust_id) as ta;"

	(columns, result) = executeQuery(cursor, query)
	cursor2csv(columns, result, afolder + "/t_clusterPresenceAbsence_analysis2.csv")

	print "# INFO: Runtime analysis of t_clusterVScluster: ", (datetime.now()-startTime)
	print "# INFO: wrote results to t_clusterPresenceAbsence_analysis2.csv"

	#os.system("sed -i \"s/analysisNr/a2/g\" " + afolder + "/tmp.R")
	#os.system("sed \"s/afolder/" + afolderString + "/g\" " + Rpath + "/clusterPresenceAbsence_analysis2.R > " + afolder + "/tmp.R")
	#os.system("sed -i \"s/infile/" + file_string + "/g\" " + afolder + "/tmp.R")
	#os.system("R CMD BATCH " + afolder + "/tmp.R " + afolder + "/test.out ")


if __name__ == '__main__':
    main()
#!/usr/bin/python

'''##################################################################
# Imports
##################################################################'''
import sys,os
sys.path.append(os.path.join(os.path.dirname(__file__), "../utils/"))
from aspmine_imports import *

'''##################################################################
# Code uses local functions: 
##################################################################

CustomArgumentParser
DBconnect
executeQuery

##################################################################'''


'''##################################################################
# ARGUMENTS and setup
##################################################################'''
startTime = datetime.now() # record runtime
parser = CustomArgumentParser(formatter_class=SmartFormatter, usage='%(prog)s -dbname [database name]')
parser.add_argument("-dbname", "-d", required=False, default = "aspminedb", help="Database name")

parser.add_argument("-clean1", "-c1", required=False, action='store_true', help="Compile new antismash2organism")
parser.add_argument("-clean2", "-c2", required=False, action='store_true', help="Compile new antismash2blast")
parser.add_argument("-clean3", "-c3", required=False, action='store_true', help="Compile new antismash2blast_reduced")
parser.add_argument("-ana1", "-a1", required=False, action='store_true', help="Cluster VS organism")
parser.add_argument("-ana2", "-a2", required=False, action='store_true', help="Cluster VS Cluster")

parser.add_argument("-csize", "-cs", required=False, default = "5", help="Cluster size cutoff - default 5")

args 	= parser.parse_args()
clean1 = args.clean1
clean2 = args.clean2
clean3 = args.clean3
csize = args.csize
ana1 = args.ana1
ana2 = args.ana2
startTime = datetime.now() # record runtime

if clean2:
	clean3=True
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
# Analysis, cluster VS organism\t\t: %s\n\
# Analysis, cluster VS cluster\t\t: %s\n\
#--------------------------------------------------------------" % (args.dbname, clean1, clean2, clean3, csize, ana1, ana2)

'''------------------------------------------------------------------
# Connection to database
------------------------------------------------------------------'''
cursor = DBconnect("localhost", "asp", "1234", args.dbname) # costum function

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
------------------------------------------------------------------

antismash2blast

For each gene, how many hits are found in other organisms:
select org_id, clust_id, name, sm_protein_id h_org, count(*) from antismash2blast_reduced group by sm_protein_id, h_org

Each antismash gene and all its BLAST hits in otehr organisms where h_cov+q_cov >= 180

Get best hit for each antismash gene in each hit organism (maximum pident after the cutoff of h_cov+q_cov >= 180)
- saved in antismash2blast_reduced

For each antismash cluster - how many genes have a hit in another organism?
select org_id, clust_id, name, h_org, count(*) from antismash2blast_reduced group by clust_id, h_org
select *, ROUND(nrHits/clust_size, 2) as clustCov from (select org_id, clust_id, name, h_org, count(*) as nrHits, SUBSTRING_INDEX(clust_id, '_', -1) as clust_size 
from antismash2blast_reduced group by clust_id, h_org ) as ta ;

For each antismash gene - how many genes have a hit in another organism?
# VOID question! For the reduced table we have selected the best hit to a gene in each organism, by definition only one gene per organism!
# Only if two genes in orgnaism X are equally good matches will we have more than one per q_seq-h_org combo
!!! select org_id, clust_id, name, sm_protein_id h_org, count(*) from antismash2blast_reduced group by sm_protein_id, h_org

There is 8685 such cases, and theya re all cases where there is two equally good matches
select count(*), avg(pident), min(pident), max(pident) from (
	select org_id, clust_id, name, sm_protein_id, h_org, count(*) as c, pident from antismash2blast_reduced group by sm_protein_id, h_org ) as ta 
where c > 1;
+----------+-------------+-------------+-------------+
| count(*) | avg(pident) | min(pident) | max(pident) |
+----------+-------------+-------------+-------------+
|     8685 |   66.412028 |       25.71 |      100.00 |
+----------+-------------+-------------+-------------+

mysql> select count(*), avg(pident), min(pident), max(pident) from (
	select org_id, clust_id, name, sm_protein_id, h_org, count(*) as c, pident from antismash2blast_reduced group by sm_protein_id, h_org ) as ta 
where c = 1;
+----------+-------------+-------------+-------------+
| count(*) | avg(pident) | min(pident) | max(pident) |
+----------+-------------+-------------+-------------+
|   597806 |   67.346514 |       24.03 |      100.00 |
+----------+-------------+-------------+-------------+

For each antismash cluster - how many organisms have a hit?


For the hits, how many are found in a cluste rin the organism theya re found in?
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






"""

select antismash2blast_reduced.org_id as q_org_id, clust_id, antismash2blast_reduced.name, q_seqkey, organism.org_id as h_org_id, h_org, h_seqkey,
SUBSTRING_INDEX(clust_id, '_', -1) as clust_size 
from antismash2blast_reduced 
join organism on h_org=organism.name 
join antismash on
limit 10;


select tc.*, cluster from (
SELECT ta.*, O1.section as sec1, O2.section as sec2, O2.org_id as h_org_id  \
from (select org_id, clust_id, SUBSTRING_INDEX(clust_id, '_', -1) as clust_size, name, q_seqkey, h_org, h_seqkey from antismash2blast_reduced ) as ta \
join organism O1 using (name)\
join organism O2 on (O2.name=h_org)\
where clust_size >= 10 and h_org = "Aspni7" and h_org != "Afoetidus" and O1.section="Nigri" and O2.section="Nigri" ) as tc
left join antismash on antismash.org_id=h_org_id and h_seqkey=antismash.sm_protein_id
order by tc.clust_id
 limit 10;


select tc.*, cluster from (
SELECT ta.*, ROUND(nrHits/clust_size, 2) as clustCov, O1.section as sec1, O2.section as sec2, O2.org_id as h_org_id  \
from (select org_id, clust_id, sm_protein_id, name, h_org, h_seqkey, count(*) as nrHits,\
SUBSTRING_INDEX(clust_id, '_', -1) as clust_size from antismash2blast_reduced group by clust_id, h_seqkey ) as ta \
join organism O1 using (name)\
join organism O2 on (O2.name=h_org)\
where clust_size >= 10 and h_org = "Aspni7" and h_org != "Afoetidus" and O1.section="Nigri" and O2.section="Nigri" ) as tc
left join antismash on antismash.org_id=h_org_id and h_seqkey=antismash.sm_protein_id
order by tc.clust_id
 limit 10;


select * from (
select ta.* from antismash2blast as ta\
join ( select q_seqid, h_seqid, h_org, max(pident) as max from antismash2blast group by q_seqid, h_org ) as tb \
on ta.q_seqid=tb.q_seqid and ta.h_org=tb.h_org and  pident=max) as tc group by q_seqid, h_org

select * from (
SELECT ta.*, ROUND(nrHits/clust_size, 2) as clustCov, O1.section as sec1, O2.section as sec2  \
from (select org_id, clust_id, name, h_org, count(*) as nrHits,\
SUBSTRING_INDEX(clust_id, '_', -1) as clust_size from antismash2blast_reduced group by clust_id, h_org ) as ta \
join organism O1 using (name)
join organism O2 on (O2.name=h_org)) as b where clustCov > 1;

select tb.*, length(prot_seq) from (
select org_id, sm_protein_id, q_seqid, h_seqid, q_cov, h_cov, pident from (
select org_id, sm_protein_id from antismash where org_id=2 and clust_backbone=433825) as ta 
left join blast on q_org="Aspacu1" and h_org="Aspacu1" and sm_protein_id=q_seqkey) as tb
join proteins on prot_seqkey=sm_protein_id and tb.org_id=proteins.org_id;




mysql> select sm_protein_id, q_seqid, h_seqid, q_cov, h_cov, pident from (
	select org_id, sm_protein_id from antismash where clust_size = 11 and org_id=10 and clust_backbone=460533) as ta 
left join blast on q_org="Aspell1" and h_org="Aspell1" and sm_protein_id=q_seqkey;
+---------------+--------------------------------------------------------------------+--------------------------------------------------------------------+--------+--------+--------+
| sm_protein_id | q_seqid                                                            | h_seqid                                                            | q_cov  | h_cov  | pident |
+---------------+--------------------------------------------------------------------+--------------------------------------------------------------------+--------+--------+--------+
|         56456 | jgi|Aspell1|56456|CE56455_165                                      | jgi|Aspell1|56456|CE56455_165                                      | 100.00 | 100.00 | 100.00 |
|         56476 | jgi|Aspell1|56476|CE56475_553                                      | jgi|Aspell1|329731|e_gw1.118.63.1                                  |  57.98 |  96.76 |  48.49 |
|         56476 | jgi|Aspell1|56476|CE56475_553                                      | jgi|Aspell1|56476|CE56475_553                                      | 100.00 | 100.00 | 100.00 |
|        287513 | NULL                                                               | NULL                                                               |   NULL |   NULL |   NULL |
|        360028 | jgi|Aspell1|360028|estExt_Genewise1.C_1580064                      | jgi|Aspell1|360028|estExt_Genewise1.C_1580064                      | 100.00 | 100.00 | 100.00 |
|        401628 | jgi|Aspell1|401628|fgenesh1_kg.158_#_27_#_Locus9618v1rpkm5.42      | jgi|Aspell1|401628|fgenesh1_kg.158_#_27_#_Locus9618v1rpkm5.42      | 100.00 | 100.00 | 100.00 |
|        401629 | NULL                                                               | NULL                                                               |   NULL |   NULL |   NULL |
|        401631 | jgi|Aspell1|401631|fgenesh1_kg.158_#_30_#_Locus14544v2rpkm0.00_PRE | jgi|Aspell1|401631|fgenesh1_kg.158_#_30_#_Locus14544v2rpkm0.00_PRE | 100.00 | 100.00 | 100.00 |
|        421779 | jgi|Aspell1|421779|fgenesh1_pm.158_#_26                            | jgi|Aspell1|421779|fgenesh1_pm.158_#_26                            | 100.00 | 100.00 | 100.00 |
|        433159 | NULL                                                               | NULL                                                               |   NULL |   NULL |   NULL |
|        460533 | jgi|Aspell1|460533|MIX5791_115_51                                  | jgi|Aspell1|426414|gm1.1982_g                                      |  43.24 |  38.26 |  28.88 |
|        460533 | jgi|Aspell1|460533|MIX5791_115_51                                  | jgi|Aspell1|436157|gm1.11725_g                                     |  34.14 |  45.34 |  37.59 |
|        460533 | jgi|Aspell1|460533|MIX5791_115_51                                  | jgi|Aspell1|460533|MIX5791_115_51                                  | 100.00 | 100.00 | 100.00 |
|        486514 | jgi|Aspell1|486514|estExt_Genemark1.C_1580039                      | jgi|Aspell1|477196|MIX22454_10_50                                  |  90.43 |  90.42 |  43.73 |
|        486514 | jgi|Aspell1|486514|estExt_Genemark1.C_1580039                      | jgi|Aspell1|486514|estExt_Genemark1.C_1580039                      | 100.00 | 100.00 | 100.00 |
+---------------+--------------------------------------------------------------------+--------------------------------------------------------------------+--------+--------+--------+


select tb.*, length(prot_seq) from (
select org_id, sm_protein_id, q_seqid, h_seqid, q_cov, h_cov, pident from (
select org_id, sm_protein_id from antismash where clust_size = 11 and org_id=10 and clust_backbone=460533) as ta 
left join blast on q_org="Aspell1" and h_org="Aspell1" and sm_protein_id=q_seqkey) as tb
join proteins on prot_seqkey=sm_protein_id and tb.org_id=proteins.org_id;
+--------+---------------+--------------------------------------------------------------------+--------------------------------------------------------------------+--------+--------+--------+------------------+
| org_id | sm_protein_id | q_seqid                                                            | h_seqid                                                            | q_cov  | h_cov  | pident | length(prot_seq) |
+--------+---------------+--------------------------------------------------------------------+--------------------------------------------------------------------+--------+--------+--------+------------------+
|     10 |         56456 | jgi|Aspell1|56456|CE56455_165                                      | jgi|Aspell1|56456|CE56455_165                                      | 100.00 | 100.00 | 100.00 |              217 |
|     10 |         56476 | jgi|Aspell1|56476|CE56475_553                                      | jgi|Aspell1|329731|e_gw1.118.63.1                                  |  57.98 |  96.76 |  48.49 |              614 |
|     10 |         56476 | jgi|Aspell1|56476|CE56475_553                                      | jgi|Aspell1|56476|CE56475_553                                      | 100.00 | 100.00 | 100.00 |              614 |
|     10 |        287513 | NULL                                                               | NULL                                                               |   NULL |   NULL |   NULL |               67 |
|     10 |        360028 | jgi|Aspell1|360028|estExt_Genewise1.C_1580064                      | jgi|Aspell1|360028|estExt_Genewise1.C_1580064                      | 100.00 | 100.00 | 100.00 |              497 |
|     10 |        401628 | jgi|Aspell1|401628|fgenesh1_kg.158_#_27_#_Locus9618v1rpkm5.42      | jgi|Aspell1|401628|fgenesh1_kg.158_#_27_#_Locus9618v1rpkm5.42      | 100.00 | 100.00 | 100.00 |              424 |
|     10 |        401629 | NULL                                                               | NULL                                                               |   NULL |   NULL |   NULL |               79 |
|     10 |        401631 | jgi|Aspell1|401631|fgenesh1_kg.158_#_30_#_Locus14544v2rpkm0.00_PRE | jgi|Aspell1|401631|fgenesh1_kg.158_#_30_#_Locus14544v2rpkm0.00_PRE | 100.00 | 100.00 | 100.00 |              199 |
|     10 |        421779 | jgi|Aspell1|421779|fgenesh1_pm.158_#_26                            | jgi|Aspell1|421779|fgenesh1_pm.158_#_26                            | 100.00 | 100.00 | 100.00 |              229 |
|     10 |        433159 | NULL                                                               | NULL                                                               |   NULL |   NULL |   NULL |              110 |
|     10 |        460533 | jgi|Aspell1|460533|MIX5791_115_51                                  | jgi|Aspell1|426414|gm1.1982_g                                      |  43.24 |  38.26 |  28.88 |             2285 |
|     10 |        460533 | jgi|Aspell1|460533|MIX5791_115_51                                  | jgi|Aspell1|436157|gm1.11725_g                                     |  34.14 |  45.34 |  37.59 |             2285 |
|     10 |        460533 | jgi|Aspell1|460533|MIX5791_115_51                                  | jgi|Aspell1|460533|MIX5791_115_51                                  | 100.00 | 100.00 | 100.00 |             2285 |
|     10 |        486514 | jgi|Aspell1|486514|estExt_Genemark1.C_1580039                      | jgi|Aspell1|477196|MIX22454_10_50                                  |  90.43 |  90.42 |  43.73 |              564 |
|     10 |        486514 | jgi|Aspell1|486514|estExt_Genemark1.C_1580039                      | jgi|Aspell1|486514|estExt_Genemark1.C_1580039                      | 100.00 | 100.00 | 100.00 |              564 |
+--------+---------------+--------------------------------------------------------------------+--------------------------------------------------------------------+--------+--------+--------+------------------+




select org_id, clust_id, name, h_org, count(*) as nrHits,\
SUBSTRING_INDEX(clust_id, '_', -1) as clust_size from antismash2blast_reduced where clust_id="10_460533_11" group by clust_id, h_org
10_460533_11




SELECT *, ROUND(nrHits/clust_size, 2) as clustCov from (
	select org_id, clust_id, name, h_org, count(*) as nrHits,SUBSTRING_INDEX(clust_id, '_', -1) as clust_size from antismash2blast_reduced 
	group by clust_id, h_org ) as ta 
where org_id=10;


SELECT * from (
SELECT distinct ab.sm_protein_id as abp, abr.sm_protein_id as abrp from antismash2blast as ab \
left join ( SELECT distinct(sm_protein_id) from antismash2blast_reduced) as abr using (sm_protein_id) \
union \
SELECT distinct ab.sm_protein_id abp, abr.sm_protein_id as abrp from antismash2blast as ab \
right join ( SELECT distinct(sm_protein_id) from antismash2blast_reduced) as abr using (sm_protein_id)) as tmp \
where abrp is NULL or abp is NULL;


select abp, abrp from 
(SELECT distinct(sm_protein_id) as abp from antismash2blast) as ta
left join 
(SELECT distinct(sm_protein_id) as abrp from antismash2blast_reduced) as tb
on abp=abrp where (abrp is NULL)

 SELECT distinct abp, abrp from  (
 SELECT distinct(sm_protein_id) as abp from antismash2blast) as ta 
 left join (SELECT distinct(sm_protein_id) as abrp from antismash2blast_reduced) as tb 
 right join (SELECT distinct(sm_protein_id) as abrp from antismash2blast_reduced) as tc 
 on abp=abrp where abrp is NULL



 select ta.* from antismash2blast as ta\
join ( select q_seqid, h_seqid, h_org, max(pident) as max from antismash2blast where sm_protein_id=530716 group by q_seqid, h_org) as tb \
on ta.q_seqid=tb.q_seqid and ta.h_org=tb.h_org and  pident=max
"""





















# SELECT data from table
"""
# BEST match for each anti-gene in each organism
select q_seqid, h_seqid, h_org, max(pident) as max from antismash2blast group by q_seqid, h_org

select ta.pident tp, ta.q_seqid aqs, blast.q_seqid bqs, blast.pident bp, max

CREATE TABLE antismash2blast_reduced as 
select ta.q_org tqo, ta.q_seqkey tqs, ta.q_seqid tqid, ta.q_cov tqc, ta.h_org tho, ta.h_seqkey ths, ta.h_seqid thid, ta.h_cov thc, ta.pident tp,
blast.q_org bqo, blast.q_seqkey bqs, blast.q_seqid bqid, blast.q_cov bqc, blast.h_org bho, blast.h_seqkey bhs, blast.h_seqid bhid, blast.h_cov bhc, blast.pident bp,
tb.q_seqid mqid, tb.h_seqid mhid, tb.h_org mho, max
from blast 
right join antismash2blast as ta
on blast.h_seqid = ta.q_seqid and blast.q_seqid=ta.h_seqid  
left join (select q_seqid, h_seqid, h_org, max(pident) as max from antismash2blast group by q_seqid, h_org) as tb
on blast.pident=max and ta.q_seqid=tb.q_seqid and ta.h_org=tb.h_org;


CREATE TABLE antismash2blast_reducedbest as \
SELECT * FROM antismash2blast 

SELECT count(*) from antismash2blast as tb join
( SELECT q_seqid, h_seqid, max(pident) as max FROM antismash2blast where q_org = "Aspbr1" GROUP BY q_seqid, h_seqid) as ta
on tb.q_seqid=ta.q_seqid and tb.h_seqid=ta.h_seqid and tb.pident=ta.max limit 10;

SELECT * from antismash2blast as tb inner join
( SELECT q_seqid, h_org, max(pident) as max FROM antismash2blast GROUP BY q_seqid, h_org) as ta
on tb.q_seqid=ta.q_seqid and tb.h_org=ta.h_org and tb.pident=ta.max limit 10;

# reduce by taking best match for each q_seqkey and h_seqkey combo

# reduce hits by cutoff, coverage
# add coverage columns

(r.q_end-r.q_start)/r.q_len*100 as qcov, (r.h_end-r.h_start)/r.h_len*100 as scov

alter table blast add column q_cov decimal(10,2);
alter table blast add column h_cov decimal(10,2);
update blast set q_cov = ((q_end-q_start+1)/q_len*100);
update blast set h_cov = ((h_end-h_start+1)/h_len*100);

select count(*) from (SELECT q_seqid, h_org, max(pident) as max FROM antismash2blast GROUP BY q_seqid, h_org) as ta

"""


"""
select ta.pident, ta.q_seqid aqs, ta.h_seqid ahs, blast.q_seqid bqs, blast.h_seqid bhs, blast.pident  

select ta.pident, ta.q_seqid aqs, blast.q_seqid bqs, blast.pident  
from blast 
right join antismash2blast as ta
on blast.h_seqid = ta.q_seqid and blast.q_seqid=ta.h_seqid 
limit 10;


select 
ta.q_org aqo, ta.q_seqkey aqs, ta.h_org as aho, ta.h_seqkey ahs,  
blast.q_org bqo, blast.q_seqkey bqs, blast.h_org bho, blast.h_seqkey bhs,  
blast.pident, ta.pident 
from blast 
right join (
select q_seqkey, h_seqkey, q_org, h_org, q_seqid, h_seqid, pident, h_cov, q_cov from antismash2blast  
where h_cov+q_cov > 180) as ta 
on blast.h_seqid = ta.q_seqid and blast.q_seqid=ta.h_seqid 
limit 10;

select max(blast.pident), blast.q_seqid, blast.h_org,  ta.q_org, ta.h_seqid
from blast 
right join (
select q_seqkey, h_seqkey, q_org, h_org, q_seqid, h_seqid, pident, h_cov, q_cov from antismash2blast  
where h_cov+q_cov > 180) as ta 
on blast.h_seqid = ta.q_seqid and blast.q_seqid=ta.h_seqid 
group by blast.q_seqid, blast.h_org
limit 10;

select count(*) from antismash2blast; 3895404
select count(*), avg(h_cov), max(h_cov), min(h_cov), avg(q_cov), max(q_cov), min(q_cov), avg(pident), min(pident), max(pident) ,
count(distinct(q_seqid)), count(distinct(h_seqid))
from antismash2blast 
where h_cov > 90 and q_cov > 90; 
+----------+------------+------------+-------------+-------------+-------------+
| count(*) | avg(h_cov) | avg(q_cov) | avg(pident) | min(pident) | max(pident) |
+----------+------------+------------+-------------+-------------+-------------+
|  1216953 |  97.394853 |  97.497061 |   53.499076 |       24.02 |      100.00 |
+----------+------------+------------+-------------+-------------+-------------+

select count(*), avg(h_cov), max(h_cov), min(h_cov), avg(q_cov), max(q_cov), min(q_cov), avg(pident), min(pident), max(pident) ,
count(distinct(q_seqid)), count(distinct(h_seqid))
from antismash2blast 
where h_cov+q_cov > 180; 
+----------+------------+------------+-------------+-------------+-------------+
| count(*) | avg(h_cov) | avg(q_cov) | avg(pident) | min(pident) | max(pident) |
+----------+------------+------------+-------------+-------------+-------------+
|  1386997 |  96.677236 |  96.840297 |   52.461977 |       23.67 |      100.00 |
+----------+------------+------------+-------------+-------------+-------------+

 select count(*), ta.pident, ta.q_seqid aqs, blast.q_seqid bqs, blast.pident   
 from blast  right join ( 
 select q_seqkey, h_seqkey, q_org, h_org, q_seqid, h_seqid, pident, h_cov, q_cov from antismash2blast   where h_cov+q_cov > 180) as ta  
on blast.h_seqid = ta.q_seqid and blast.q_seqid=ta.h_seqid  limit 10;
+----------+--------+---------------------------------------------+----------------------------------------------------------+--------+
| count(*) | pident | aqs                                         | bqs                                                      | pident |
+----------+--------+---------------------------------------------+----------------------------------------------------------+--------+
|  1456900 |  73.19 | jgi|Aspfo1|223345|estExt_Genemark1.C_180044 | jgi|PenbrAgRF18_1|312792|estExt_Genewise1Plus.C_6_t10313 |  73.87 |
+----------+--------+---------------------------------------------+----------------------------------------------------------+--------+


select q_org, q_seqkey, organism.org_id 
from antismash2blast 
left join organism on (q_org = organism.name);

select * from (
select q_org, q_seqkey, h_seqkey, organism.org_id as orgid
from antismash2blast 
left join organism on (h_org = organism.name)
) as ta 
left join antismash on h_seqkey = antismash.sm_protein_id and 
ta.orgid = antismash.org_id;


Some hist are found several times for a query sequence
select count(*) from (select * from antismash2blast group by q_seqid, h_org) as ta;
select count(*) from (select * from antismash2blast group by q_seqid, h_seqid) as ta;	2924416
select count(distinct q_seqid, h_seqid) from antismash2blast ;	2924416
select count(*) from antismash2blast;	3895404

"""








"""
query = "SELECT distinct antismash2blast.name, antismash2organism.name from antismash2blast \
		right join antismash2organism using (name) where antismash2organism.name is null;"
(columns, lackingOrg) = executeQuery(cursor, query)

query = "SELECT distinct antismash2blast.name, antismash2organism.name from antismash2blast \
		right join antismash2organism using (name) where antismash2organism.name is null;"
(columns, lackingOrg) = executeQuery(cursor, query)

if lackingAnti:
	orgs = [o[1] for o in lackingAnti]
	print ("# INFO: There are orgs in antismash2organism not found in antismash2blast %s" % orgs )

if lackingOrg:
	orgs = [o[1] for o in lackingOrg]
	sys.exit("# INFO: There are orgs in antismash2blast not found in antismash2organism, %s" % orgs)
"""

"""


CREATE TABLE antismash2organism as select ta.org_id, sm_protein_id, clust_id, name 
from ( select org_id, sm_protein_id, concat(org_id, '_', clust_backbone, '_', clust_size) as clust_id  from antismash) as ta 
join organism using (org_id)

CREATE TABLE antismash2blast as select * from antismash2organism join blast on (sm_protein_id = q_seqkey and q_org = name);

CREATE TABLE antismash2blast as 
select * from antismash2organism 
join blast on (sm_protein_id = q_seqkey and q_org = name)

CREATE TABLE antiblast as 
select tb.*, max_pident from antitmp3 as tb
inner join(
select sm_protein_id, blast_sseq_jg1,
max(blast_pident) as max_pident 
from antitmp3 
group by sm_protein_id, blast_sseq_jg1
) as ta
on tb.sm_protein_id = ta.sm_protein_id and 
ta.max_pident = tb.blast_pident and 
tb.blast_sseq_jg1 = ta.blast_sseq_jg1;


n antismash not found in organism, %s" % (nrOrgs -nrOrgsAll))
"""
"""
+---------------------------+-----------------------+
| count(distinct(clust_id)) | count(distinct(name)) |
+---------------------------+-----------------------+
|                      3143 |                    46 |
+---------------------------+-----------------------+
+-----------------------+
| count(distinct(name)) |
+-----------------------+
|                    50 |
+-----------------------+
"""
"""
select * from antismash2organism 
join blast on (sm_protein_id = blast_qseq_jg2 and blast_qseq_jg1 = name)
limit 10;"""
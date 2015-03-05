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
parser.add_argument("-clean1", "-c1", required=False, action='store_true', help="Compile new t_antismash2organism")
parser.add_argument("-clean2", "-c2", required=False, action='store_true', help="Compile new t_antismash2blast")
#parser.add_argument("-stats", "-s", required=False, action='store_true', help="Compile stats table")

parser.add_argument("--out", "-o", required=False, default = "org_shared_gc.csv", help="Name of output file")
parser.add_argument("-R", default="org_shared_gc.R", required=False, help="Name of R script, org_shared_gc.R")

args 	= parser.parse_args()
clean1 = args.clean1
clean2 = args.clean2
outfile = args.out
rscript = args.R
startTime = datetime.now() # record runtime

'''------------------------------------------------------------------
# Print argument values to screen
------------------------------------------------------------------'''
print "#--------------------------------------------------------------\n\
# ARGUMENTS\n\
# Database\t: %s\n\
# Compile new t_antismash2organism\t: %s\n\
# Compile new t_antismash2blast\t: %s\n\
#--------------------------------------------------------------" % (args.dbname, clean1, clean2)

'''------------------------------------------------------------------
# Connection to database
------------------------------------------------------------------'''
cursor = DBconnect("localhost", "asp", "1234", args.dbname) # costum function

'''------------------------------------------------------------------
# Create connection table: antismash -> organism
------------------------------------------------------------------'''
print "# INFO: Connection antismash with organism - t_antismash2organism"

# CREATE table if specified
if clean1:
	print "# INFO: Dropping and recreating t_antismash2organism"

	cursor.execute("DROP TABLE IF EXISTS t_antismash2organism ")
	query = "CREATE TABLE t_antismash2organism as \
	select ta.org_id, sm_protein_id, clust_id, name \
	from ( select org_id, sm_protein_id, concat(org_id, '_', clust_backbone, '_', clust_size) as clust_id  from antismash) as ta \
	join organism using (org_id);"

	(columns, result)	 = executeQuery(cursor, query)

	executeQuery(cursor, "CREATE INDEX `sm_protein_id` ON `t_antismash2organism` (`sm_protein_id`,`name`)")
	executeQuery(cursor, "CREATE INDEX `org_id`  ON `t_antismash2organism`  (`org_id`)")
	executeQuery(cursor, "CREATE INDEX `name`  ON `t_antismash2organism` (`name`)")
	executeQuery(cursor, "CREATE INDEX `sm_protein_id_2`  ON `t_antismash2organism` (`sm_protein_id`)")

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
	print ("# INFO: There are orgs in organism not found in antismash %s" % orgs )

if lackingOrg:
	orgs = [o[1] for o in lackingOrg]
	sys.exit("# INFO: There are orgs in antismash not found in organism, %s" % orgs)

'''------------------------------------------------------------------
# Create connection table: antismash-organism -> blast
------------------------------------------------------------------'''
print "# INFO: Connection t_antismash2organism with BLAST - t_antismash2blast, this will take a while"

"""
select org_id, sm_protein_id,  concat(org_id, "_", clust_backbone, "_", clust_size) as clust_id  from antismash
CREATE TABLE t_antismash2blast as select * from t_antismash2organism join blast on (sm_protein_id = q_seqkey and q_org = name);
Query OK, 3895404 rows affected (2 min 47.47 sec)
"""

# CREATE table if specified
if clean2:
	print "# INFO: Dropping and recreating t_antismash2blast"

	cursor.execute("DROP TABLE IF EXISTS t_antismash2blast")
	query = "CREATE TABLE t_antismash2blast as select * from t_antismash2organism join blast on (sm_protein_id = q_seqkey and q_org = name);"

	(columns, result) = executeQuery(cursor, query)

	print "# INFO: Createing index"
	executeQuery(cursor, "CREATE INDEX `sm_protein_id` ON `t_antismash2blast` (`sm_protein_id`,`name`)")
	executeQuery(cursor, "CREATE INDEX `name`  ON `t_antismash2blast` (`name`)")
	executeQuery(cursor, "CREATE INDEX `sm_protein_id_2`  ON `t_antismash2blast` (`sm_protein_id`)")
	executeQuery(cursor, "CREATE INDEX `seqids`  ON `t_antismash2blast` (`h_seqid`, `q_seqid`)")

# If not specified to create table test if the table actually does exist
if not cursor.execute("Show tables LIKE 't_antismash2blast'"):
	sys.exit("# ERROR: table does not exist, re-run with -clean2 option")

# SELECT data from table

CREATE TABLE t_antismash2blast_reduced as \
SELECT * FROM t_antismash2blast 

SELECT count(*) from t_antismash2blast as tb join
( SELECT q_seqid, h_seqid, max(pident) as max FROM t_antismash2blast where q_org = "Aspbr1" GROUP BY q_seqid, h_seqid) as ta
on tb.q_seqid=ta.q_seqid and tb.h_seqid=ta.h_seqid and tb.pident=ta.max limit 10;

SELECT count(*) from t_antismash2blast as tb inner join
( SELECT q_seqid, h_seqid, max(pident) as max FROM t_antismash2blast where q_org = "Aspbr1" GROUP BY q_seqid, h_seqid) as ta
on tb.q_seqid=ta.q_seqid and tb.h_seqid=ta.h_seqid and tb.pident=ta.max limit 10;

SELECT q_seqid, h_seqid, max(pident[,h_org]) as max FROM t_antismash2blast where q_org = "Aspbr1" GROUP BY q_seqid, h_seqid
# reduce by taking best match for each q_seqkey and h_seqkey combo

# reduce hits by cutoff, coverage
# add coverage columns









"""
select q_org, q_seqkey, organism.org_id 
from t_antismash2blast 
left join organism on (q_org = organism.name);

select * from (
select q_org, q_seqkey, h_seqkey, organism.org_id as orgid
from t_antismash2blast 
left join organism on (h_org = organism.name)
) as ta 
left join antismash on h_seqkey = antismash.sm_protein_id and 
ta.orgid = antismash.org_id;


Some hist are found several times for a query sequence
select count(distinct q_seqid, h_seqid) from t_antismash2blast ; 2924416
select count(*) from t_antismash2blast; 3895404

"""








"""
query = "SELECT distinct t_antismash2blast.name, t_antismash2organism.name from t_antismash2blast \
		right join t_antismash2organism using (name) where t_antismash2organism.name is null;"
(columns, lackingOrg) = executeQuery(cursor, query)

query = "SELECT distinct t_antismash2blast.name, t_antismash2organism.name from t_antismash2blast \
		right join t_antismash2organism using (name) where t_antismash2organism.name is null;"
(columns, lackingOrg) = executeQuery(cursor, query)

if lackingAnti:
	orgs = [o[1] for o in lackingAnti]
	print ("# INFO: There are orgs in t_antismash2organism not found in t_antismash2blast %s" % orgs )

if lackingOrg:
	orgs = [o[1] for o in lackingOrg]
	sys.exit("# INFO: There are orgs in t_antismash2blast not found in t_antismash2organism, %s" % orgs)
"""

"""


CREATE TABLE t_antismash2organism as select ta.org_id, sm_protein_id, clust_id, name 
from ( select org_id, sm_protein_id, concat(org_id, '_', clust_backbone, '_', clust_size) as clust_id  from antismash) as ta 
join organism using (org_id)

CREATE TABLE t_antismash2blast as select * from t_antismash2organism join blast on (sm_protein_id = q_seqkey and q_org = name);

CREATE TABLE t_antismash2blast as 
select * from t_antismash2organism 
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
select * from t_antismash2organism 
join blast on (sm_protein_id = blast_qseq_jg2 and blast_qseq_jg1 = name)
limit 10;"""
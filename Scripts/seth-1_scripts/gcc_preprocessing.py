################################
# Imports ######################
################################

import sys,os
sys.path.append("/home/seth/300asp/utils")
from aspmine_imports import *
import pandas as pd
from seth_imports import *



parser = argparse.ArgumentParser(description="Redo tables and stuff", usage="%(prog)s")
#parser.add_argument("--plot", "-p", required=False, default = False, help="Only execute plotting")
parser.add_argument("--redo_bidir", "-rb", required=False, action = 'store_true', help="Redo t_genecluster_bidir_hits table")
parser.add_argument("--redo_tables", "-rt", required=False, action = 'store_true', help="Redo tables ")
parser.add_argument("--singletons", "-s", required=False, action = 'store_true', help="Get singletons")

args = parser.parse_args()
redo_tables = args.redo_tables
redo_bidir = args.redo_bidir
singletons = args.singletons

cursor = asp_con(path='192.38.13.196', user='setd', pw='1234')


####################
# Table management #
####################

if redo_tables:
	query = "DROP TABLE IF EXISTS t_genecluster_testing_full"

	cursor.execute(query)

	query = "CREATE table t_genecluster_testing_full as\
	select *,  (buf3/2)+(pident*members/q_size/200) as custom_score from\
	(select *, count(*) as members, max(buf2) as buf3 from\
		(select *, case when q_sm_short != 'none' or h_sm_short != 'none' then  buf+1 else buf+0 end as buf2\
			from t_genecluster_bidir_hits) tw\ #deleted: where pident >50
			group by q_clust,clust_id) tx where members >= 0.25*q_size and members >= 0.25*h_size and q_size >= 2 and h_size >=2;"
	
# Can be used to count: select count(distinct(q_clust)), count(distinct(clust_id)) from t_genecluster_testing_full;

# select * from t_genecluster_testing where q_clust = clust_id and org_id != 3 and q_size != members; This finds out if some members are not found in a self hit

	cursor.execute(query)

	# old version: erase buf3/2 in the custom score and the 200 as well

	query = "DROP TABLE IF EXISTS t_gcc_testing_full"
	cursor.execute(query)
# Here a count does not really make sense since I group them by q_clust and h_clust gets picked randomly
	query = "CREATE table t_gcc_testing_full as\
	select * from (select *, count(*) as family_hits\
		from (select *,  pident*members/q_size as custom_score\
			from (select *, count(*) as members from t_genecluster_bidir_hits\
				group by q_clust,clust_id) tx\ # deleted: where pident > 50
				where members >= 0.25*q_size and members >= 0.25*h_size and q_size >= 2 and h_size >=2) ty group by q_clust) tz where family_hits >10;"
	cursor.execute(query)


# This is generating the initial bidirectional hits table between geneclusters
if redo_bidir:
	query = "DROP TABLE IF EXISTS t_genecluster_bidir_hits_full"

	cursor.execute(query)

	print "Deleted table t_genecluster_bidir_hits_full"

	query = "CREATE table t_genecluster_bidir_hits_full as\
	select tx.org_id, tx.clust_id as q_clust, tx.clust_size as q_size, tx.sm_protein_id as q_gene, tx.sm_short as q_sm_short, tx.h_seqkey as h_gene, ty.sm_short as h_sm_short, ty.clust_id, ty.clust_size as h_size, tx.bitscore, tx.pident\
	from t_antismash2blast as tx join t_antismash2blast as ty on tx.h_seqkey = ty.sm_protein_id where ty.h_seqkey = tx.sm_protein_id AND pident >= 50;"

	cursor.execute(query)

	print "Created new t_genecluster_bidir_hits_full table"



#########################
# Extracting Singletons #
#########################

# Singletons are only hitting themselves and nothing else, so they should be extracted before other tables are generated
if singletons:
	query = "select q_clust from (select *, count(*) as family_hits from (select *,  pident*members/q_size as custom_score from (select *, count(*) as members from t_genecluster_bidir_hits where pident > 50 group by q_clust,clust_id) tx where members >= 0.5*q_size and members >= 0.5*h_size and q_size > 1) ty group by q_clust) tz where family_hits <= 1 and q_clust = clust_id;"

	cursor.execute(query)
	singletons = list(cursor.fetchall())

	print singletons

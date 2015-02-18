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
parser.add_argument("-cluster", "-c", required=False, action='store_true', help="Compile cluster table")
parser.add_argument("-stats", "-s", required=False, action='store_true', help="Compile stats table")

args 	= parser.parse_args()
cluster = args.cluster
stats 	= args.stats

startTime = datetime.now() # record runtime

'''------------------------------------------------------------------
# Print argument values to screen
------------------------------------------------------------------'''
print "#--------------------------------------------------------------\n\
# ARGUMENTS\n\
# Database\t: %s\n\
#--------------------------------------------------------------" % args.dbname

if not cluster and not stats:
	sys.exit("# INFO: no tables were selected, to see options, run with -h, exiting")

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


'''------------------------------------------------------------------
# antismash_cluster_type_counts
------------------------------------------------------------------'''
if cluster:
	print "# INFO: compiling antismash_cluster_type_counts"
	try:
		query = "select name as JGI_code, real_name as SpeciesName, \
		count(PKS) as nrPKS, count(PKSlike) as nrPKSlike, \
		count(NRPS) as nrNRPS, count(NRPSlike) as nrNRPDlike, \
		count(TC) as nrTC, count(HYBRID) as nrHYBRID, count(DMAT) as nrDMAT \
		from ( \
			select org_id,filename, clust_backbone, \
			case when clust_type='PKS' then clust_type end as PKS, \
			case when clust_type='NRPS' then clust_type end as NRPS, \
			case when clust_type='PKS-like' then clust_type end as PKSlike, \
			case when clust_type='TC' then clust_type end as TC, \
			case when clust_type='HYBRID' then clust_type end as HYBRID, \
			case when clust_type='DMAT' then clust_type end as DMAT, \
			case when clust_type='NRPS-like' then clust_type end as NRPSlike \
			from antismash where  sm_protein_id = clust_backbone \
		) as ta join organism using (org_id) group by org_id;"
		cursor.execute(query)
		result = cursor.fetchall()

	except mdb.Error, e:
			sys.exit("# ERROR: failed to execute antismash_cluster_type_counts, %d: %s" % (e.args[0],e.args[1]))

	print "# INFO: Writing to antismash_cluster_type_counts.csv"

	f = open("antismash_cluster_type_counts.csv",'wb')
	writer = csv.writer(f, dialect = 'excel', delimiter=";")
	columns = map(lambda x:x[0], cursor.description) 	
	writer.writerow(columns)
	writer.writerows(result)
	f.close()
	print "# INFO Runtime antismash_cluster_type_counts: ", (datetime.now()-startTime)


'''------------------------------------------------------------------
# genome_statistics
------------------------------------------------------------------'''
if stats:
	print "# INFO: compiling genome_statistics"
	try:
		query = "select JGIcode, name,\
		FORMAT(ROUND(nrProteins,0),0) as nrProteins,\
		ROUND(avgProteins,2) as avgProteinLength,\
		FORMAT(maxProteins,0) as maxProteinLength,\
		FORMAT(minProteins,0) as minProteinLength,\
		ROUND(bDNA/1000000,2) as GenomeMb \
		from ( select organism.org_id as JGIcode, name,\
		(select count(prot_seqkey) from proteins where org_id=organism.org_id group by org_id) as nrProteins,\
		(select avg(length(prot_seq)) from proteins where org_id=organism.org_id group by org_id) as avgProteins,\
		(select max(length(prot_seq)) from proteins where org_id=organism.org_id group by org_id) as maxProteins,\
		(select min(length(prot_seq)) from proteins where org_id=organism.org_id group by org_id) as minProteins,\
		(select sum(length(assembly_seq)) from assembly where org_id = organism.org_id) as bDNA\
		from organism group by organism.org_id) as ta;"
		cursor.execute(query)
		result = cursor.fetchall()

	except mdb.Error, e:
			sys.exit("# ERROR: failed to execute antismash_cluster_type_counts, %d: %s" % (e.args[0],e.args[1]))

	print "# INFO: writing to genome_statistics.csv"

	f = open("genome_statistics.csv",'wb')
	writer = csv.writer(f, dialect = 'excel', delimiter=";")
	columns = map(lambda x:x[0], cursor.description) 	
	writer.writerow(columns)
	writer.writerows(result)
	f.close()

	print "# INFO Runtime genome_statistics: ", (datetime.now()-startTime)

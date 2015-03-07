################################
# Imports ######################
################################

import sys,os
sys.path.append("../utils/")
from aspmine_imports import *


#################################
# Arguments #####################
#################################
# action sotre true means argument is boolean, cannot assign values
parser = argparse.ArgumentParser(description="Preprocessing databases for gene clustering", usage="%(prog)s --out filename")
parser.add_argument("--out", "-o", required=False, default = "testing.csv", help="Name of output file")
#parser.add_argument("-R", default="org_shared_gc.R", required=False, help="Name of R script, org_shared_gc.R")
#parser.add_argument("--sec", required=False, default = "Nigri", help="Section name")
parser.add_argument("--plot","-p", required=False, action='store_true', help="Executes R plotting functions on data")
parser.add_argument("--cutoff1", "-c1", type = float, required=False, help="Select a cutoff for gene cluster similarity, i.e. how many genes must be shared among two gene clusters in order to call them homolog/ the same")
parser.add_argument("--cutoff2", "-c2", type = int, required=False, default = 0, help="Select a cutoff for gene cluster sizes, defaults to 0")
parser.add_argument("--ana1", "-a1", required=False, action='store_true', help="Overview of Clusters")
parser.add_argument("--ana2", "-a2", required=False, action='store_true', help="Clusters shared among organisms")
parser.add_argument("--ana3", "-a3", required=False, action='store_true', help="Find which clusters are found in most organisms")
parser.add_argument("--bestclusters", "-bc", required=False, action='store_true', help="Find which clusters are found in most organisms, TOP 20, remember cutoffs!")
parser.add_argument("--uniqueclusters", "-uc", required=False, action='store_true', help="Find which clusters are unique among organisms, remember cutoffs!")
parser.add_argument("--max_gc_abundance", "-max", required=False, action='store_true', help="Calculates max values for gene clusters")
parser.add_argument("--normalize", "-d", required=False, action='store_true', help="Normalizes results from analysis 2 for heatmap")
#TODO: Fix type and choices

args = parser.parse_args()

outfile = args.out
cutoff1 = args.cutoff1
cutoff2 = args.cutoff2
plot = args.plot
ana1 = args.ana1
ana2 = args.ana2
ana3 = args.ana3
best = args.bestclusters
unique = args.uniqueclusters
max_gc_abundance = args.max_gc_abundance
done = args.done



section = 'Nigri' # TODO Needs implementation in mysql, here the table only has Nigri members....


print cutoff1
print cutoff2
print section
# IMPORTANT!  DELETE FROM antismashLoopAntismash where q_orgname = 'Aspac1' and h_orgname = 'Aspac1'; bc Aspac 1 entries are still weird...
def make_plot(outfile):
	print "# INFO: running Rscript"
	try:
		os.system("R CMD BATCH '--args %s %s %s %s' gc_shared_through_orgs.R test.out" % (outfile, section, cutoff1*100, cutoff2))
	except:
		print "Cannot connect to R"

# CONNECTING

def asp_con(path):
#------------------------------------------------------------------
# Connect to specific DB
#------------------------------------------------------------------
	try:
		if path == 'remote':
			db = mdb.connect(host="192.38.13.9", user="setd" ,passwd="1234",db= 'aspminedb') #host="localhost", user="root" ,passwd="Gru3n3r+T33",db= database)
	except mdb.Error, e:
		sys.exit("Cannot connect to database")#"# ERROR %d: %s" % (e.args[0],e.args[1]))
	try:
		cursor = db.cursor()
	except mdb.Error, e:
		sys.exit("# ERROR %d: %s" % (e.args[0],e.args[1]))
	return cursor

cursor = asp_con('remote')

if max_gc_abundance:
	query ="CREATE TABLE testing_max as SELECT q_orgname, max(shared) as maxsh from (SELECT q_orgname, h_orgname, count(*) as shared from(SELECT q_clustid, q_clust_size, q_orgname, h_orgname, h_clust_backbone, h_clust_protein_id, members from (SELECT *, count(*) as members from antismashLoopAntismash where h_clust_backbone is not null group by q_clustid, q_orgname, h_orgname, h_clust_backbone) ta where members > q_clust_size*%f and q_clust_size>%i) tb group by q_orgname ,h_orgname) tc group by q_orgname;" % (cutoff1, cutoff2)

if ana1:
	query = "SELECT q_clustid, q_clust_size, q_orgname, h_orgname, h_clust_backbone, h_clust_protein_id, members from (\
		SELECT *, count(*) as members from antismashLoopAntismash\
		where h_clust_backbone is not null group by q_clustid, q_orgname, h_orgname, h_clust_backbone) ta\
		where members > q_clust_size*%f and q_clust_size>%i;" % (cutoff1, cutoff2)

if ana2:
	query =" CREATE TABLE t_gc_shared_count as SELECT q_orgname, h_orgname, count(*) as shared from(\
		SELECT q_clustid, q_clust_size, q_orgname, h_orgname, h_clust_backbone, h_clust_protein_id, members from (\
		SELECT *, count(*) as members from antismashLoopAntismash\
		where h_clust_backbone is not null group by h_clust_backbone) ta\
		where members > q_clust_size*%f and q_clust_size>%i) tb group by q_orgname ,h_orgname;" % (cutoff1, cutoff2)

if best:
	query ="SELECT q_clustid, count(*) as found from (SELECT * from antismashLoopAntismashCandidates\
		where h_clust_backbone is not null and clustCov > %f and q_clust_size>%i group by h_clust_backbone) ta\
		group by q_clustid order by found DESC limit 20;" % (cutoff1, cutoff2)

if unique:
	query = " SELECT * from (SELECT q_clustid, count(*) as found from (\
		SELECT * from antismashLoopAntismashCandidates where h_clust_backbone is not null and clustCov > %f and q_clust_size>%i\
		 group by h_clust_backbone) ta group by q_clustid order by found) tb where found =1;" % (cutoff1, cutoff2)

# TODO finish here!
if done:
	query = "SELECT t_gc_shared_count.q_orgname, t_gc_shared_count.h_orgname, (t_gc_shared_count.shared/testing_max.maxsh)*100 as norm_shared  from t_gc_shared_count left join testing_max ON t_gc_shared_count.q_orgname = testing_max.q_orgname ;"

print query

try:
	cursor.execute(query)
	result = cursor.fetchall()
except:
	print "Error in Analysis"

f = open(outfile,'wb')
writer = csv.writer(f, dialect = 'excel')
writer.writerows(result)
f.close()

'''------------------------------------------------------------------
# Plotting
------------------------------------------------------------------'''

if plot:
	make_plot(outfile)



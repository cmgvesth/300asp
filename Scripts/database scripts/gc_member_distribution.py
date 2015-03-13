import sys,os
sys.path.append(os.path.join("../../utils/")) #sys.path.append(os.path.join(os.path.dirname(__file__), "../../utils/"))
from aspmine_imports import *
import pandas as pd
outfile = 'gc_asframe.csv'

def gc_preprocess():

#------------------------------------------------------------------
# Connect to specific DB
#------------------------------------------------------------------
	try:
		db = mdb.connect(host="192.38.13.9", user="setd" ,passwd="1234",db="aspminedb")
	except mdb.Error, e:
		sys.exit("Cannot connect to database")#"# ERROR %d: %s" % (e.args[0],e.args[1]))
	try:
		cursor = db.cursor()
	except mdb.Error, e:
		sys.exit("# ERROR %d: %s" % (e.args[0],e.args[1]))

	try:
		query = """SELECT h_org, found_in_gc from (select ta.h_org, ta.h_seqkey,\
			antismash.sm_protein_id, antismash.org_id, ta.clust_id as clust_id_x, concat(antismash.org_id, "_", antismash.clust_backbone, "_", antismash.clust_size ) as clust_id_y, count(*) as found_in_gc, blast_qseq_id  from (\
		 	select h_org, h_seqkey, blast_qseq_id, organism.org_id as orgid, clust_id from antiblast\
		 	left join organism on (h_org = organism.name)) as ta left join antismash on h_seqkey = antismash.sm_protein_id and ta.orgid = antismash.org_id where org_id!=0 group by clust_id_y) tc;
		"""
		cursor.execute(query)
		result = cursor.fetchall()

	except mdb.Error, e:
		sys.exit('Setting up initial table failed')

	genome = {}

	for line in result:
		a,b = line
		if a in genome:
			genome[a].append(int(b))
		else:
			genome[a] = list()
			genome[a].append(int(b))
		
# To do: If runs slow, improve here by filling up with NAs before to avoid transposing
	gframe = pd.DataFrame.from_dict(genome, orient='index', dtype=float)
	gframe = gframe.transpose()
	gframe.to_csv(outfile)



def gc_plot(outfile):
	print "# INFO: running Rscript"
	os.system("R CMD BATCH '--args %s' gc_plot.R " % outfile)

make_table()

make_plot(outfile)

#	f = open("gc_processed",'wb')
#	writer = csv.writer(f, dialect = 'excel')
#	writer.writerows(result)
#	f.close()
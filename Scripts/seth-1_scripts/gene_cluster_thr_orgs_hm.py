import sys,os
sys.path.append("../utils/")
from aspmine_imports import *
from seth_imports import *
import pandas as pd
from genecluster_viewer import *

parser = argparse.ArgumentParser(description="Analyze distribution of gene cluster members among species with heatmap.2", usage="%(prog)s")
#parser.add_argument("--plot", "-p", required=False, default = False, help="Only execute plotting")
parser.add_argument("--gc_id", "-g", required=False, type = str, default = '11_382281_6', help="Specify gene cluster for analysis. Add the s option!")
parser.add_argument("--batch_mode", "-batch", type = str, required=False, help="Specify csv file for batch gene cluster for analysis, needs genecluster ids as first element")
parser.add_argument("--cutoff", "-c", type = int, required=False, help="Minimum amount of members found in a hit-genecluster")
parser.add_argument("--geneclusterview", "-view", action = 'store_true', required=False, help="Uses Genoediagram to create a gene viewer of the cluster. Remember to specify a cutoff as well!")
parser.add_argument("--single_gc", "-s", action = 'store_true', required=False, help="Add if only a single genecluster should be analyzed")

args = parser.parse_args()
gc_id = args.gc_id
batch_mode = args.batch_mode
member_cutoff = args.cutoff
view = args.geneclusterview
single_gc = args.single_gc
#plot = args.plot

def gc_mapping(gc_id = gc_id):
	#if not plot:
	if gc_id == '11_382281_6':
		print 'Showing example for gene cluster 11_382281_6. \n To analyze your gene cluster of interest please provide a valid genecluster id with gc_id = ...'
	else:
		print 'Fetching data for gene cluster %s' % gc_id

	cursor = asp_con(path='192.38.13.196', user='setd', pw='1234')

	#################
	# Fetching data#
	################

	# Hit clusters are created and their blast percent identity extracted
# Attempting ordering so replaced this piece of code for real_name: concat(real_name, ' : ', tx.h_clust_id)
	query = "SELECT tx.q_id, real_name as hc_id, h_id, pident, tx.org_id FROM (\
		SELECT a2b.sm_protein_id as q_id,  concat(smash.org_id, '_', smash.clust_backbone, '_', smash.clust_size) AS h_clust_id, smash.sm_protein_id as h_id, a2b.pident as pident, smash.org_id\
		from t_antismash2blast_reduced AS a2b left JOIN antismash AS smash on h_seqkey = smash.sm_protein_id\
		WHERE a2b.clust_id = '%s' and smash.org_id is not NULL\
		order by org_id, h_clust_id, a2b.sm_protein_id) tx left join organism using (org_id) where h_clust_id in (SELECT h_clust_id FROM(SELECT h_clust_id, count(*) as new FROM (SELECT a2b.sm_protein_id as q_id,  concat(smash.org_id, '_', smash.clust_backbone, '_', smash.clust_size) AS h_clust_id, smash.sm_protein_id as h_id, a2b.pident as pident, smash.org_id from t_antismash2blast_reduced AS a2b left JOIN antismash AS smash on h_seqkey = smash.sm_protein_id WHERE a2b.clust_id = '%s' and smash.org_id is not NULL order by org_id, h_clust_id, a2b.sm_protein_id) tg group by h_clust_id) th where new > %s) ;" % (gc_id, gc_id, member_cutoff)
	cursor.execute(query)

	data = list(cursor.fetchall())

	df = pd.DataFrame(data, columns=['q_id', 'h_clust_id','h_id', 'pident','org_id'])

	storage_name = 'gc_m_distr_%s.csv' % gc_id

	df.to_csv(storage_name)

	print "Saved file as %s" % storage_name




	#####################################################
	# Creating syntenic order using seqrecord locations #
	#####################################################
	
	print "Building SeqRecord object for synteny of genecluster %s" % gc_id

	record = single_seq_builder(gc_id)

	l = []
	for i in record.features:
		l.append((i.location.start.real,i.qualifiers['locus_tag']))

	sor_members = sorted(l,key=lambda x:x[0])

	l_members = [x[1] for x in sor_members]

	l_members.insert(0, gc_id) # Appending gene cluster name to list for use in R

	l_members.insert(1, member_cutoff) # Appending cutoff to list for use in R


	nargs = '%s '*len(l_members)

	r_query = "R CMD BATCH '--args %s ' gc_map.R test.out" % nargs

	print r_query

	os.system(r_query % tuple(l_members))

####################
# Executing script #
####################

if batch_mode:
	with open('%s' % batch_mode) as f:
		reader = csv.reader(f)
		hits = [x[0] for x in reader]#list(reader)
	if single_gc:
		for i in hits:
			gc_mapping(str(i))
	if view:
		for i in hits:
			geneclusterview(str(i), member_cutoff)

if single_gc:
	gc_mapping(gc_id)


if view and not batch_mode:
	geneclusterview(gc_id, member_cutoff)
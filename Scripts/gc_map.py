import sys,os
sys.path.append("../utils/")
from aspmine_imports import *
from seth_imports import *
import pandas as pd


cursor = asp_con(path='192.38.13.196', user='setd', pw='1234')
#a2b.clust_id,
# old version query = "SELECT a2b.sm_protein_id as q_id,  concat(smash.org_id, '_', smash.clust_backbone, '_', smash.clust_size) AS h_clust_id, smash.sm_protein_id as h_id, a2b.pident as pident, smash.org_id from t_antismash2blast_reduced AS a2b left JOIN antismash AS smash on h_seqkey = smash.sm_protein_id WHERE a2b.clust_id = '11_382281_6' and smash.org_id is not NULL order by org_id, h_clust_id, a2b.sm_protein_id;"
query = "SELECT tx.q_id, concat(real_name, ' : ', tx.h_clust_id) as hc_id, h_id, pident, tx.org_id FROM (SELECT a2b.sm_protein_id as q_id,  concat(smash.org_id, '_', smash.clust_backbone, '_', smash.clust_size) AS h_clust_id, smash.sm_protein_id as h_id, a2b.pident as pident, smash.org_id from t_antismash2blast_reduced AS a2b left JOIN antismash AS smash on h_seqkey = smash.sm_protein_id WHERE a2b.clust_id = '11_382281_6' and smash.org_id is not NULL order by org_id, h_clust_id, a2b.sm_protein_id) tx left join organism using (org_id);"
cursor.execute(query)

data = list(cursor.fetchall())

df = pd.DataFrame(data, columns=['q_id', 'h_clust_id','h_id', 'pident','org_id'])

df.to_csv('thmap.csv')

record = single_seq_builder()

l = []
for i in record.features:
	l.append((i.location.start.real,i.qualifiers['locus_tag']))

sor_members = sorted(l,key=lambda x:x[0])

l_members = [x[1] for x in sor_members]

nargs = '%s '*len(l_members)

r_query = "R CMD BATCH '--args %s ' gc_map.R test.out" % nargs

os.system(r_query % tuple(l_members))

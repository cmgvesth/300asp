select org_id, sm_protein_id,  concat(org_id, "_", clust_backbone, "_", clust_size) as clust_id  from antismash




CREATE TABLE antismash2organism as
select ta.org_id, sm_protein_id, clust_id, name 
from ( select org_id, sm_protein_id,  
concat(org_id, "_", clust_backbone, "_", clust_size) as clust_id  from antismash) as ta 
join organism using (org_id);



select * from antismash2organism 
join blast on (sm_protein_id = blast_qseq_jg2 and blast_qseq_jg1 = name)
limit 10;
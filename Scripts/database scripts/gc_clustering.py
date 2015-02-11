db = "aspminedb"
"""
select distinct antiorg, org, name, blast_qseq_jg1 from (select distinct antismash.org_id as antiorg ,organism.org_id as org, name  from organism left join antismash using (org_id)) as ta
left join (select distinct blast_qseq_jg1 from blast) as tb on (blast_qseq_jg1=name) ; 
CREATE TABLE antitmp as select org_id, sm_protein_id, concat(org_id, "_", clust_backbone, "_", clust_size) as clust_id from antismash where clust_backbone='1079950' and org_id=22;

CREATE TABLE antitmp2 as select antitmp.org_id, sm_protein_id, clust_id, name from antitmp join organism using (org_id);


proper way to create table:
create table antitmp3 as select * from antitmp2 join blast on (sm_protein_id = blast_qseq_jg2 and blast_qseq_jg1 = name);



CREATE TABLE antiblast as select tb.*, max_pident from antitmp3 as tb
inner join(
	select sm_protein_id, blast_sseq_jg1, max(blast_pident) as max_pident from antitmp3 group by sm_protein_id, blast_sseq_jg1
) as ta
on tb.sm_protein_id = ta.sm_protein_id and ta.max_pident = tb.blast_pident and tb.blast_sseq_jg1 = ta.blast_sseq_jg1;



 select blast_sseq_jg1, blast_sseq_jg2, organism.org_id from antiblast left join organism on (blast_sseq_jg1 = organism.name);

 select * from (select blast_sseq_jg1, blast_sseq_jg2, organism.org_id as orgid from antiblast left join organism on (blast_sseq_jg1 = organism.name)) as ta left join antismash on blast_sseq_jg2 = antismash.sm_protein_id and ta.orgid = antismash.org_id;

 select ta.blast_sseq_jg1, ta.blast_sseq_jg2, ta.blast_qseq_id, antismash.sm_protein_id, antismash.org_id from (select blast_sseq_jg1, blast_sseq_jg2, blast_qseq_id, organism.org_id as orgid from antiblast left join organism on (blast_sseq_jg1 = organism.name)) as ta left join antismash on blast_sseq_jg2 = antismash.sm_protein_id and ta.orgid = antismash.org_id;


"""

select count(*) from antismash where org_id = 41;

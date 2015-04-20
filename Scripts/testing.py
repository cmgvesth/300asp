SELECT h_clust_id, real_name FROM (SELECT h_clust_id, h_org, In_cluster FROM (\
		SELECT org_id, In_cluster, h_clust_id, h_org, new FROM (\
		SELECT h_org, GROUP_CONCAT(h_seqkey) as new,a2b.clust_id AS q_clust_id,\
		concat(smash.org_id, '_', smash.clust_backbone, '_', smash.clust_size) AS h_clust_id, smash.org_id, count(*) AS In_cluster from t_antismash2blast_reduced AS a2b\
		JOIN antismash AS smash on h_seqkey = smash.sm_protein_id WHERE a2b.clust_id = '%s' and smash.org_id GROUP BY h_clust_id) ta ORDER BY org_id, In_cluster DESC) tx\
		JOIN organism on tx.org_id = organism.org_id) ty WHERE In_cluster > 2) ty LEFT JOIN organism on h_org = name;
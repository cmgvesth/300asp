notes.py



asubdat$total <- NULL
amdat <- melt(asubdat)
amdat$n <- ordered(asubdat$SpeciesName, levels = rev(levels(asubdat$SpeciesName)))
myColors <- brewer.pal(8,"GnBu")

names(myColors) <-levels(amdat$variable)

pdf('antismash_cluster_type_counts.pdf', width=8.5, height=8)
a2 <- ggplot(amdat, aes(x = n, y=value, fill=variable)) +
geom_bar(stat="identity") +
coord_flip() +
scale_fill_manual(name = "Backbone type",values = myColors, labels=c("PKS", "PKS like", "NRPS", "NRPD", "TC", "Hybrid", "DMAT"))+
theme(axis.text = element_text(size = 12, colour="grey20")) +
ggtitle("Nr. and type of predicted SM clusters") +
xlab("Organism") +
ylab("Nr. sec. metabolism clusters")
a2
dev.off()


notes.txt

Selecting gene annotations for gene clusters

select ta.*, prot_orgkey, prot_tailkey, prot_seqname from (select  org_id, sm_protein_id, sm_short, sm_desc from antismash where clust_backbone = 5256 and org_id = 21) ta join proteins on ta.org_id = proteins.org_id and  sm_protein_id = prot_seqkey;

Finding accession numbers if you have them

select sm_protein_id,sm_short, sm_desc, concat(proteins.org_id,'_', prot_seqkey,'_', clust_size) as clust_id from proteins left join antismash as smash on smash.org_id = proteins.org_id and prot_seqkey = clust_backbone where prot_tailkey = 'AN7909';

startTimet_antismashLoopAntismash = datetime.now()
print "# INFO Runtime t_antismashLoopAntismash: ", (datetime.now()-startTimet_antismashLoopAntismash)


sed -i 's/_/ /g' Custom.K6.FullTaxonomy.nwk 

 gsub("_"," ",taxo$tip.label)

 SELECT q_clustid, count(*) as found from (SELECT * from t_antismashLoopAntismashCandidates where h_clust_backbone is not null and clustCov > 0.7 and q_clust_size> 5 group by h_clust_backbone) ta group by q_clustid order by found DESC limit 20;


select sm_protein_id, protein_has_ipr.* from (select * from antismash where clust_backbone = 433535 and antismash.org_id = 11) ta 

select ipr.ipr_id, ipr_desc, sm_protein_id from (
	select sm_protein_id, protein_has_ipr.* from (
		select * from antismash where clust_backbone = 433535 and antismash.org_id = 11) ta
	join protein_has_ipr on sm_protein_id = protein_id and ta.org_id = protein_has_ipr.org_id) tb join ipr on tb.ipr_id = ipr.ipr_id group by ipr.ipr_id , sm_protein_id;



SELECT organism.real_name, tx.In_cluster FROM (SELECT org_id, In_cluster FROM (SELECT a2b.clust_id AS q_clust_id, concat(smash.org_id, '_', smash.clust_backbone, '_', smash.clust_size) AS h_clust_id, smash.org_id, count(*) AS In_cluster from t_antismash2blast_reduced AS a2b JOIN antismash AS smash on h_seqkey = smash.sm_protein_id WHERE a2b.clust_id = 11_433535_8 and smash.org_id GROUP BY h_clust_id) ta ORDER BY org_id, In_cluster DESC) tx join organism on tx.org_id = organism.org_id;
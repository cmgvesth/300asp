import sys, os
def asp_con(path, user, pw):
	import MySQLdb as mdb
#------------------------------------------------------------------
# Connect to specific DB
#------------------------------------------------------------------
	try:
		db = mdb.connect(host=path, user=user ,passwd=pw,db= 'aspminedb') #host="localhost", user="root" ,passwd="Gru3n3r+T33",db= database)
	except mdb.Error, e:
		sys.exit("# ERROR %d: %s" % (e.args[0],e.args[1]))
	try:
		cursor = db.cursor()
	except mdb.Error, e:
		sys.exit("# ERROR %d: %s" % (e.args[0],e.args[1]))
	return cursor

def ap_to_dict(dictionary, app):
	# TODO: Think about this function, if it makes sense
	for line in app:
		a,b = line
	if a in dictionary:
		dictionary[a].append(b)
	else:
		dictionary[a] = list()
		dictionary[a].append(b)

def single_seq_builder(gene_cluster, save = False):
	from Bio import SeqIO
	from Bio.Seq import Seq
	from Bio.SeqRecord import SeqRecord
	from Bio.Alphabet import generic_dna
	from Bio.SeqFeature import FeatureLocation, SeqFeature

	gc_id = gene_cluster

	if os.path.isfile("GC_"+gc_id+".gb") == True:
		print "Seqrecord file for %s already exists, loading file" % gc_id
		gc_record = SeqIO.read("GC_"+gc_id+".gb", "gb")
		return gc_record

	else:
		print"Creating SeqRecord object for %s" % gc_id
		print"Files are only created if you set single_seq_builder(gene_cluster, save = True)"
		cursor = asp_con('')

		
			# Fetching length of whole cluster
		query = "SELECT min(gff_start), max(gff_end) from (select org_id, sm_protein_id, concat(org_id, '_', clust_backbone, '_', clust_size) AS q_clust_id from antismash) tx join gff on tx.org_id = gff.org_id and sm_protein_id = gff_protein_id where q_clust_id = '%s' and gff_type ='CDS';" % gc_id

		cursor.execute(query)
		space = cursor.fetchall()[0]

		start, end = space
		dummy = end-start
		#query = "Get sequence FROM assembly"

		#gc_id,sequence = cursor.fetchall()
		# For now I just fill it with random things
		gc_seq = Seq("ATCG"*(dummy/4+40))

		if gc_id[1]=='_':
			organism_id = gc_id[0]
		else:
			organism_id = gc_id[:1]

		gc_record = SeqRecord(gc_seq, id=gc_id, name = organism_id) # name used for organism_id!
		gc_record.seq.alphabet = generic_dna # Assign an alphabet to the sequence (in this case DNA)

		query = "SELECT org_id, sm_protein_id from t_antismash2organism where clust_id = '%s';" % gc_id # Assigning features, leaving out exons for now

		cursor.execute(query)
		complete_cluster = cursor.fetchall()

		for i in complete_cluster:
			org, protein = i
			query = "SELECT min(gff_start), max(gff_end), gff_strand  from antismash as ta join gff on ta.org_id=gff.org_id and ta.sm_protein_id = gff_protein_id where ta.org_id = %s and ta.sm_protein_id = %s and gff.gff_type='CDS';" % (org, protein)
			cursor.execute(query)
			f_start,f_end,strand = cursor.fetchall()[0]

	# TODO for information fetch query = SELECT clust_id, sm_protein_id, GROUP_CONCAT(DISTINCT go_name) as go_desc FROM (SELECT clust_id, sm_protein_id, go_term_id from t_antismash2organism as a2o left join protein_has_go on sm_protein_id = protein_id  where clust_id = '28_10278_8') tx left join go using (go_term_id) group by sm_protein_id;


			if f_start is None:
				print "Could not get information about protein %s, organism %s from gff table" % (protein, org)
				continue
			rel_start = f_start-start
			rel_end = f_end -start
			#my_start_pos = SeqFeature.ExactPosition(rel_start) # Exactposition doesnt work for some reason.... maybe implement later
			#my_end_pos = SeqFeature.ExactPosition(rel_end) # 2. Use the locations do define a FeatureLocation
			# my_feature_location = FeatureLocation(my_start_pos,my_end_pos)
			my_feature_type = "CDS" # Could also fetch this from db 

			my_feature = SeqFeature(FeatureLocation(rel_start, rel_end),type = my_feature_type, strand = int(strand+'1'), qualifiers={'codon_start':'1', 'locus_tag':str(protein)}) # Create a SeqFeature

			gc_record.features.append(my_feature) # Append your newly created SeqFeature to your SeqRecord

			#optional: print the SeqRecord to STDOUT in genbank format, with your new feature added.
		#print "\nThis bit is the SeqRecord, printed out in genbank format, with a feature added.\n"
		#print(gc_record.format("gb"))
			#seqr_collection[gc_id] = gc_record
			#seqr_list.append(gc_record)

		if save == True:
			SeqIO.write(gc_record, "GC_"+gc_id+".gb", "gb")

		return gc_record

'''
def batch_seq_builder(q_cluster):
	from Bio import SeqIO
	from Bio.Seq import Seq
	from Bio.SeqRecord import SeqRecord
	from Bio.Alphabet import generic_dna
	from Bio.SeqFeature import FeatureLocation, SeqFeature

	cursor = asp_con('192.38.13.196', 'setd', '1234')

	query = "SELECT h_clust_id, real_name FROM (\
		SELECT * FROM (SELECT h_clust_id, iname, In_cluster FROM (\
		SELECT org_id, In_cluster, h_clust_id, iname, new FROM (\
		SELECT smash.org_id as iname, GROUP_CONCAT(h_seqkey) as new,a2b.clust_id AS q_clust_id,\
		concat(smash.org_id, '_', smash.clust_backbone, '_', smash.clust_size) AS h_clust_id, smash.org_id, count(*) AS In_cluster from t_antismash2blast_reduced AS a2b\
		JOIN antismash AS smash on h_seqkey = smash.sm_protein_id WHERE a2b.clust_id = '%s' and smash.org_id GROUP BY h_clust_id) ta ORDER BY org_id, In_cluster DESC) tx\
		JOIN organism on tx.org_id = organism.org_id) ty WHERE In_cluster > 2) tz LEFT JOIN organism on iname = organism.org_id;" % q_cluster
	# all_clusters = ['28_10013_15']
	cursor.execute(query)
# TODO CREATE TUPLE WITH ID AND NAME to sue in the for loop
	#all_clusters = []
	#temp_clusters = list(cursor.fetchall())
	all_clusters = list(cursor.fetchall())
	#print all_clusters
	#for i in temp_clusters:
	#	all_clusters.append(i[0])	
	all_clusters.insert(0,(q_cluster, 'Aspeucalypticola'))

	for i in all_clusters:
		single_seq_builder(i)

'''

def seq_builder(q_cluster, cutoff, save = False):
	from Bio import SeqIO
	from Bio.Seq import Seq
	from Bio.SeqRecord import SeqRecord
	from Bio.Alphabet import generic_dna
	from Bio.SeqFeature import FeatureLocation, SeqFeature

	if q_cluster == '11_382281_6':
		print "Showing example case for genecluster '11_382281_6'"


	####################################################
	# Fetching all gene cluster hits for a given cluster
	####################################################

	cursor = asp_con('192.38.13.196', 'setd', '1234')

	query = "SELECT h_clust_id, real_name FROM (\
		SELECT * FROM (SELECT h_clust_id, iname, In_cluster FROM (\
		SELECT org_id, In_cluster, h_clust_id, iname, new FROM (\
		SELECT smash.org_id as iname, GROUP_CONCAT(h_seqkey) as new,a2b.clust_id AS q_clust_id,\
		concat(smash.org_id, '_', smash.clust_backbone, '_', smash.clust_size) AS h_clust_id, smash.org_id, count(*) AS In_cluster from t_antismash2blast_reduced AS a2b\
		JOIN antismash AS smash on h_seqkey = smash.sm_protein_id WHERE a2b.clust_id = '%s' and smash.org_id GROUP BY h_clust_id) ta ORDER BY org_id, In_cluster DESC) tx\
		JOIN organism on tx.org_id = organism.org_id) ty WHERE In_cluster > %s) tz LEFT JOIN organism on iname = organism.org_id;" % (q_cluster, cutoff)
	# all_clusters = ['28_10013_15']
	cursor.execute(query)
# TODO CREATE TUPLE WITH ID AND NAME to sue in the for loop
	#all_clusters = []
	#temp_clusters = list(cursor.fetchall())
	all_clusters = list(cursor.fetchall())
	#print all_clusters
	#for i in temp_clusters:
	#	all_clusters.append(i[0])	
	all_clusters.insert(0,(q_cluster, 'Aspeucalypticola'))
	#print all_clusters

	dummy = []

	for i in all_clusters:
		dummy.append(i[0])

	print "('"+"','".join(dummy)+"','.*')"
	############################################
	# Creating SeqIO objects out of geneclusters
	############################################
	seqr_collection = {}
	seqr_list = []

	for gc_id, organism_name in all_clusters:
		#gc_id = str(gc_id)
		# Fetching length of whole cluster
		query = "SELECT min(gff_start), max(gff_end) from (select org_id, sm_protein_id, concat(org_id, '_', clust_backbone, '_', clust_size) AS q_clust_id from antismash) tx join gff on tx.org_id = gff.org_id and sm_protein_id = gff_protein_id where q_clust_id = '%s' and gff_type ='CDS';" % gc_id

		cursor.execute(query)
		space = cursor.fetchall()[0]

		start, end = space
		dummy = end-start
		#query = "Get sequence FROM assembly"

		#gc_id,sequence = cursor.fetchall()
		# For now I just fill it with random things
		gc_seq = Seq("ATCG"*(dummy/4+40))
		gc_record = SeqRecord(gc_seq, id=gc_id, name= gc_id, description = organism_name)
		gc_record.seq.alphabet = generic_dna # Assign an alphabet to the sequence (in this case DNA)

		query = "SELECT org_id, sm_protein_id from t_antismash2organism where clust_id = '%s';" % gc_id # Assigning features, leaving out exons for now

		cursor.execute(query)
		complete_cluster = cursor.fetchall()

		for i in complete_cluster:
			org, protein = i
			query = "SELECT min(gff_start), max(gff_end), gff_strand  from antismash as ta join gff on ta.org_id=gff.org_id and ta.sm_protein_id = gff_protein_id where ta.org_id = %s and ta.sm_protein_id = %s and gff.gff_type='CDS';" % (org, protein)
			cursor.execute(query)
			f_start,f_end,strand = cursor.fetchall()[0]

# TODO for information fetch query = SELECT clust_id, sm_protein_id, GROUP_CONCAT(DISTINCT go_name) as go_desc FROM (SELECT clust_id, sm_protein_id, go_term_id from t_antismash2organism as a2o left join protein_has_go on sm_protein_id = protein_id  where clust_id = '28_10278_8') tx left join go using (go_term_id) group by sm_protein_id;


			if f_start is None:
				print "Could not get information about protein %s, organism %s from gff table" % (protein, org)
				continue
			rel_start = f_start-start
			rel_end = f_end -start
			#my_start_pos = SeqFeature.ExactPosition(rel_start) # Exactposition doesnt work for some reason.... maybe implement later
			#my_end_pos = SeqFeature.ExactPosition(rel_end) # 2. Use the locations do define a FeatureLocation
			# my_feature_location = FeatureLocation(my_start_pos,my_end_pos)
			my_feature_type = "CDS" # Could also fetch this from db 

			my_feature = SeqFeature(FeatureLocation(rel_start, rel_end),type = my_feature_type, strand = int(strand+'1'), qualifiers={'codon_start':'1', 'locus_tag':'sm'+str(protein)}) # Create a SeqFeature

			gc_record.features.append(my_feature) # Append your newly created SeqFeature to your SeqRecord

		#optional: print the SeqRecord to STDOUT in genbank format, with your new feature added.
	#print "\nThis bit is the SeqRecord, printed out in genbank format, with a feature added.\n"
	#print(gc_record.format("gb"))
		seqr_collection[gc_id] = gc_record
		seqr_list.append(gc_record)
		if save == True:
			SeqIO.write(gc_record, "GC_"+gc_id+".gb", "gb")
		gc_record = 0
	return seqr_list

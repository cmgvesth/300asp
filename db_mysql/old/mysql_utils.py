#!/usr/bin/python
#------------------------------------------------------------------
# Imports
#------------------------------------------------------------------
import MySQLdb as mdb
import sys, getopt, argparse, re, glob, os
import gzip
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from mysql_classes import Fasta	# custom class
import pprint
import numpy
##############################################################################################################################################
##############################################################################################################################################
####################################################################### UTILS ################################################################
##############################################################################################################################################
##############################################################################################################################################

#######################################################################
# Fetch or create organism _id, auto increment 
#######################################################################
def get_org_id(org_name, dbname):
	db = mdb.connect("localhost","asp","1234",dbname) # test of connection has to be performed prior ------------------------------------
	cursor = db.cursor()

	cursor.execute("SELECT org_id from organism where name LIKE \"%" + org_name + "%\";")
	org_id = cursor.fetchone() # first element of tuple ------------------------------------

	if not org_id:
		try:
			cursor.execute("REPLACE INTO organism(name) VALUES(\"" + org_name + "\");")
			db.commit()
			org_id = cursor.lastrowid
			db.close()
			
		except mdb.Error, e:
			sys.exit("# ERROR %d: %s" % (e.args[0],e.args[1]))
		
		return org_id
	else:
		return org_id[0]

#######################################################################
# get_org_seqname_from_prot_id
#######################################################################
def get_org_seqname_from_prot_id (prot_id, org_name, dbname):
	db = mdb.connect("localhost","asp","1234",dbname)
	cursor = db.cursor()

	try:
		cursor.execute("SELECT org_id,prot_seq_name, prot_protein_id FROM proteins JOIN organism USING (org_id) WHERE prot_protein_id =%s AND name LIKE \"%s\";", (prot_id, org_name))
		result = cursor.fetchall() 
	except mdb.Error, e:
		sys.exit("# ERROR utils seqname %d: %s" % (e.args[0],e.args[1]))

	if len(result) > 1:
		sys.exit("# ERROR utils seqname: more than one sequence/organism entry was retrieved for protein _id %s" % prot_id)
	elif len(result) < 1:
		sys.exit("# ERROR utils seqname: less than one sequence/organism entry was retrieved for protein _id %s" % prot_id)

	org_id = result[0][0]
	prot_seq_name = result[0][1]
	db.close()
	
	return org_id, prot_seq_name, result

#######################################################################
# Handle file format check
#######################################################################	
def check_file_format (wanted_nr_columns, records):
	if len(records) < 1:
		sys.exit("# ERROR utils: no records found in file")
		
	for r in records:	
		if r.startswith("#") or r == "" or r == "protein_id": continue	# Exclude comment and empty lines ------------------------------------

		nr_columns = len(r.split("\t"))
		if not nr_columns == wanted_nr_columns :
			sys.exit("# ERROR utils : wrong number of columns, not correct file format, %s" % r)
			
#######################################################################
# Decompress files
#######################################################################
def gzip_file (filepath):
	records = []
	try:
		decompressed_file = gzip.open(filepath, "rb")
		tmp_records = (decompressed_file.read()).split("\n")

	except:
		sys.exit("# ERROR utils gzip: not GZIP file" )

	for i in tmp_records:
		if i.startswith("#") or i == "" or i.startswith("protein"): continue	# Exclude comment and empty lines ------------------------------------
		else: records.append(i) 

	return records

#######################################################################
# Decompress fasta files
#######################################################################
def gzip_fasta_file (filepath):
	try:
		decompressed_file = gzip.open(filepath, "rb")
		tmp_records = list(SeqIO.parse(decompressed_file , "fasta"))
	except:
		sys.exit("# ERROR utils FASTA gzip: not GZIP file" )
	return tmp_records

##############################################################################################################################################
##############################################################################################################################################
####################################################################### DATA #################################################################
##############################################################################################################################################
##############################################################################################################################################

#######################################################################
# Handle FASTA data
#######################################################################
def fasta (org_name, filepath, action, filetype, dbname, prefix):
	counter = 0
	totalcounter = 0
	org_id = ''
	
	records = gzip_fasta_file(filepath) # function to unzip file and store in list ------------------------------------
		
	# If this is empty, there was nothing that matched a FASTA format in the file ------------------------------------
	if len(list(records)) < 1:	
		sys.exit("# ERROR: Not FASTA format or empty FASTA entry")

	# Find the file name and add prefix if specified  ------------------------------------
	filename = filepath.split("/")[-1]
	if prefix: filename = prefix+filename

	valueSortLine = ("jgi1","jgi2","jgi3","seq_name", "seq", "tp_id", "org_id")

	# Protein FASTA ------------------------------------
	if filetype == "pf":
		name = "protein-FASTA"
		valueSortLine = ("jgi1","jgi2","jgi3","seq_name", "seq", "tp_id", "org_id", "seq_key")
		insertQuery = "REPLACE INTO proteins(prot_jgi1, prot_jgi2, prot_jgi3, prot_seq_name, prot_seq, prot_protein_id, org_id, prot_seq_key) values(%s,%s,%s,%s,%s,%s,%s,%s);"
		
	# DNA FASTA ------------------------------------
	elif filetype == "tf":
		name = "transcript-FASTA"
		insertQuery = "REPLACE INTO transcripts(trans_jgi1, trans_jgi2, trans_jgi3, trans_seq_name, trans_seq, trans_transcript_id, org_id) values(%s,%s,%s,%s,%s,%s,%s);"

	elif filetype == "cf":
		name = "cds-FASTA"
		insertQuery = "REPLACE INTO cds(cds_jgi1, cds_jgi2, cds_jgi3, cds_seq_name, cds_seq, cds_transcript_id, org_id) values(%s,%s,%s,%s,%s,%s,%s);"

	elif filetype == "uma" or filetype == "ma":
		name = "assembly-FASTA"
		valueSortLine = ("jgi1","jgi2","jgi3","seq_name", "seq", "org_id")
		insertQuery = "REPLACE INTO assembly(assembly_jgi1, assembly_jgi2, assembly_jgi3, assembly_seq_name, assembly_seq, org_id) values(%s,%s,%s,%s,%s,%s);"
		
	else:
		sys.exit("# ERROR: filetype not recognized %s" % filepath.split("/")[-1])
	
	return load_fasta_files(name, valueSortLine, insertQuery, org_name, filepath, action, dbname)	


#######################################################################
# Load FASTA type files
#######################################################################
def load_fasta_files(name, valueSortLine, insertQuery, org_name,filepath, action, dbname):
	if action == "load":
		print "# LOADING %s file: %s" % (name, filepath.split("/")[-1])
	else:
		print "# TESTING %s file: %s" % (name, filepath.split("/")[-1])

	records = gzip_fasta_file(filepath) # function to unzip file and store in list ------------------------------------
	if len(list(records)) < 1:	# If this is empty, there was nothing that matched a FASTA format in the file	 ------------------------------------
		sys.exit("# ERROR: Not FASTA format or empty FASTA entry")

	values_to_insert = []
	counter = 0
	totalcounter = 0
	valuehash = {}
	org_id = get_org_id(org_name, dbname) # Only get the org_id once and not for all records ------------------------------------
	
	for r in records:
		counter += 1
		totalcounter += 1
		values = ()

		# Store variable names in dict and get values from record ------------------------------------

		if name == "protein-FASTA":
			non_ExtendedIUPACProtein = set(str(r.seq.lower())) - set("ACDEFGHIKLMNPQRSTVWYBXZJUO".lower()) # Make "sure" sequence is Protein ------------------------------------ 
			# ExtendedIUPACProtein: ACDEFGHIKLMNPQRSTVWYBXZJUO http://biopython.org/DIST/docs/api/Bio.Alphabet.IUPAC.ExtendedIUPACProtein-class.html
			
			if len(non_ExtendedIUPACProtein) != 0:
				if len(non_ExtendedIUPACProtein) == 1 and "*" in non_ExtendedIUPACProtein:
					warning = "Found ambiguouse symbol: *"
				else:
					sys.exit("# ERROR : Found amino ac_ids other than ACDEFGHIKLMNPQRSTVWYBXZJUO = %s" % ",".join(non_ExtendedIUPACProtein))

		if name == "assembly-FASTA" or name == "transcirpt-FASTA" or name == "cds-FASTA": 
			non_IUPACUnambiguousDNA = set(str(r.seq.lower())) - set("atcg") # Make "sure" sequence is DNA ------------------------------------ 
			# ExtendedIUPACDNA : GATCBDSW http://biopython.org/DIST/docs/api/Bio.Alphabet.IUPAC.ExtendedIUPACDNA-class.html
			# IUPACAmbiguousDNA : GATCRYWSMKHBVDN http://biopython.org/DIST/docs/api/Bio.Alphabet.IUPAC.IUPACAmbiguousDNA-class.html
			# IUPACUnambiguousDNA : GATC http://biopython.org/DIST/docs/api/Bio.Alphabet.IUPAC.IUPACUnambiguousDNA-class.html
			
			if len(non_IUPACUnambiguousDNA) > 0:
				if len(non_IUPACUnambiguousDNA) == 1 and "n" in non_IUPACUnambiguousDNA:
					warning = "Found ambiguouse base: N"
				else:
					sys.exit("# ERROR : Found bases other than A, T, G, C or N: " + ",".join(non_IUPACUnambiguousDNA))

		valuehash["seq_name"]	= (r.description,)
		#valuehash["id_line"]	= (r.description.split(" ")[0],)
		if len(r.description.split(" ")[0].split("|"))>3:
			valuehash["jgi1"]		= (r.description.split(" ")[0].split("|")[1],)
			valuehash["jgi2"]		= (r.description.split(" ")[0].split("|")[2],)
			valuehash["jgi3"]		= (r.description.split(" ")[0].split("|")[3],)
		else:
			valuehash["jgi1"]		= (r.description.split(" ")[0],)
			valuehash["jgi2"]		= ("",)
			valuehash["jgi3"]		= ("",)
		
		valuehash["seq"]		= (str(r.seq.lower()),)
		valuehash["name"]		= ("_".join(r.description.split("_")[:-1]),)
		valuehash["tp_id"]		= ('',)
		valuehash["seq_key"]	= ('',)
		
		if len(r.description.split("|")) > 1:
			valuehash["tp_id"]	=	(r.description.split("|")[2],)
			valuehash["seq_key"] = 	(r.description.split("|")[1]+"_"+r.description.split("|")[2],)	 # jgi|Aspac1|10177|gw1.1.1170.1

		for val in valueSortLine:	# create list of values to be inserted ------------------------------------
			if val == "org_id":	values += (org_id,)
			else: 				values += valuehash[val]
			
		values_to_insert.append(values) # Create sets of lists of values ------------------------------------

		counter_threshold = 1000
		# MySQL will only allow a certain data insert size so for assemblies 
		if name == "assembly-FASTA" : counter_threshold = 2 
		
		if action == "load" and ( counter == counter_threshold or totalcounter == len(records)):
			if totalcounter % 3000 == 0  or totalcounter == len(records) : print "# INFO: Inserting record number %s" % totalcounter

			try:	
				db = mdb.connect("localhost","asp","1234",dbname)
				cursor = db.cursor()
				cursor.executemany(insertQuery,values_to_insert)
				db.commit()
				db.close()
				
				counter = 0 
				values_to_insert = []

			except mdb.Error, e:
				sys.exit( "# ERROR utils %s %d: %s" % (name, e.args[0],e.args[1]))

	print "# FINISHED %s file: %s with total number of records %s" % (name,filepath.split("/")[-1], totalcounter)
	return org_id	

##############################################################################################################################################
##############################################################################################################################################
####################################################################### ANNOTATIONS ##########################################################
##############################################################################################################################################
##############################################################################################################################################


#######################################################################
# Load tab type files
#######################################################################
def load_tab_files(name, nrtab, valueLine, valueSortLine, insertQuery, org_name,filepath, action, dbname):
	if action == "load":
		print "# LOADING %s file: %s" % (name, filepath.split("/")[-1])
	else:
		print "# TESTING %s file: %s" % (name, filepath.split("/")[-1])

	records = gzip_file(filepath) # function to unzip file and store in list ------------------------------------
	check_file_format (nrtab, records) # function to verify that the tab file has the right number of fields ------------------------------------
	
	values_to_insert = []
	counter = 0
	totalcounter = 0
	valuehash = {}
	org_id = get_org_id(org_name, dbname) # Only get the org_id once and not for all records ------------------------------------

	
	for r in records:
		counter += 1
		totalcounter += 1
		values = ()
	
		for var in range(0,len(valueLine)):	# Store variable names in dict and get values from record ------------------------------------
			valuehash[valueLine[var]] = (r.split("\t")[var],)
			#print valuehash
			
		#print valuehash["gff_attributes"][0]
		for val in valueSortLine:	# create list of values to be inserted ------------------------------------
			if name == "GFF":
				string = str(valuehash["gff_attributes"][0])
			if val == "org_id":
				values += (org_id,)

			elif val == "gff_name":
				values += (valuehash["gff_attributes"][0].split("\"")[1],)

			elif val == "gff_protein_id":
				if (re.search(r"proteinId (\d+)", string)):
					values += (re.search(r"proteinId (\d+)", string).group(1),)
				else:
					values += ('',)
			elif val == "gff_transcript_id":
				if (re.search(r"transcriptId (\d+)", string)):
					values += (re.search(r"transcriptId (\d+)", string).group(1),)
				else:
					values += ('',)
			else:	
				values += valuehash[val]

		print values_to_insert
		values_to_insert.append(values) # Create sets of lists of values ------------------------------------
		
		if action == "load" and ( counter == 1000 or totalcounter == len(records)):
			# Print insert information 
			if totalcounter % 10000 == 0 : print "# INFO: Inserting record number %s" % totalcounter
			elif totalcounter % 2000 == 0: print "# INFO: Inserting record number %s" % totalcounter
			elif totalcounter == len(records): print "# INFO: Inserting record number %s" % totalcounter
						
			try:	
				db = mdb.connect("localhost","asp","1234",dbname)
				cursor = db.cursor()
				cursor.executemany(insertQuery,values_to_insert)

				db.commit()
				db.close()
				
				counter = 0 
				values_to_insert = []

			except mdb.Error, e:
				sys.exit( "# ERROR utils %s %d: %s" % (name, e.args[0],e.args[1]))

	print "# FINISHED %s file: %s with total number of records %s" % (name,filepath.split("/")[-1], totalcounter)
	#print "Org_id:", org_id
	return org_id	
	
#######################################################################
# GO
#######################################################################
def go (org_name,filepath, action, dbname):	
	name = "GO"
	nrtab = 5
	valueLine = ["protein_id", "goterm_id", "goName", "gotermType", "goAcc"]
	valueSortLine = ("protein_id", "goterm_id", "goName", "gotermType", "goAcc", "org_id", "protein_id", "org_id")
	insertQuery = "REPLACE INTO go (go_protein_id, go_term_id, go_name, go_termtype, go_acc, org_id, go_prot_seq_name) values(%s,%s,%s,%s,%s,%s,(SELECT prot_seq_name from proteins WHERE prot_protein_id =%s and org_id = %s));"
	return load_tab_files(name, nrtab, valueLine, valueSortLine, insertQuery, org_name,filepath, action, dbname)


#######################################################################
# SignalP
#######################################################################
def sigp (org_name,filepath, action, dbname):	
	name = "SignalP"
	nrtab = 5
	valueLine = ["protein_id", "nn_cutpos", "neuro_net_vote", "hmm_cutpos", "hmm_signalpep_probability"]
	valueSortLine = ("protein_id", "nn_cutpos", "neuro_net_vote", "hmm_cutpos", "hmm_signalpep_probability", "org_id", "protein_id", "org_id")
	insertQuery = "REPLACE INTO sigp (sigp_protein_id, sigp_nn_cutpos, sigp_neuro_net_vote, sigp_hmm_cutpos, sigp_hmm_signalpep_probability, org_id, sigp_prot_seq_name) VALUES(%s,%s,%s,%s,%s,%s,(SELECT prot_seq_name from proteins WHERE prot_protein_id =%s and org_id = %s));"
	return load_tab_files(name, nrtab, valueLine, valueSortLine, insertQuery, org_name,filepath, action, dbname)

#######################################################################
# InterPro
#######################################################################
def ipr (org_name,filepath, action, dbname):	
	name = "InterPro"
	nrtab = 10
	valueLine = ["protein_id", "ipr_id", "ipr_desc", "domaindb", "domain_id", "domaindesc", "numHits", "domainStarts", "domainEnds", "score"]
	valueSortLine = ("protein_id", "ipr_id", "ipr_desc", "domaindb", "domain_id", "domaindesc", "numHits", "domainStarts", "domainEnds", "score", "org_id", "protein_id", "org_id")
	insertQuery = "REPLACE INTO ipr (ipr_protein_id, ipr_id, ipr_desc, ipr_domaindb, ipr_domain_id, ipr_domaindesc, ipr_numHits, ipr_domain_starts, ipr_domain_ends, ipr_score, org_id, ipr_prot_seq_name) values(%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,(SELECT prot_seq_name from proteins WHERE prot_protein_id =%s and org_id = %s));"
	return load_tab_files(name, nrtab, valueLine, valueSortLine, insertQuery, org_name,filepath, action, dbname)

#######################################################################
# KEGG
#######################################################################
def kegg (org_name,filepath, action, dbname):	
	name = "KEGG"
	nrtab = 9
	valueLine = ["protein_id", "ecNum", "definition", "catalyticActivity", "cofactors", "associatedDiseases", "pathway", "pathway_class", "pathway_type"]
	valueSortLine = ("protein_id", "ecNum", "definition", "catalyticActivity", "cofactors", "associatedDiseases", "pathway", "pathway_class", "pathway_type", "org_id", "protein_id", "org_id")
	insertQuery = "REPLACE INTO kegg (kegg_protein_id, kegg_ecNum, kegg_definition, kegg_catalyticActivity, kegg_cofactors, kegg_associatedDiseases, kegg_pathway, kegg_pathway_class, kegg_pathway_type, org_id, kegg_prot_seq_name) VALUES(%s,%s,%s,%s,%s,%s,%s,%s,%s,%s , (SELECT prot_seq_name from proteins WHERE prot_protein_id =%s and org_id = %s));"
	return load_tab_files(name, nrtab, valueLine, valueSortLine, insertQuery, org_name,filepath, action, dbname)

#######################################################################
# KOG
#######################################################################
def kog (org_name,filepath, action, dbname):	
	name = "KOG"
	nrtab = 6
	valueLine = ["transcript_id", "protein_id", "kog_id", "kogdefline", "e", "b"]
	valueSortLine = ("transcript_id", "protein_id", "kog_id", "kogdefline", "org_id", "protein_id", "org_id")
	insertQuery = "REPLACE INTO kog (kog_transcript_id, kog_protein_id, kog_id, kog_defline, org_id, kog_prot_seq_name) values (%s,%s,%s,%s, %s , (SELECT prot_seq_name from proteins WHERE prot_protein_id =%s and org_id = %s));"
	return load_tab_files(name, nrtab, valueLine, valueSortLine, insertQuery, org_name,filepath, action, dbname)

#######################################################################
# GFF
#######################################################################
def gff (org_name,filepath, action, dbname):	
	name = "GFF"
	nrtab = 9
	valueLine = ["gff_seqorigin", "a","gff_type", "gff_start", "gff_end", "gff_score", "gff_strand", "gff_phase", "gff_attributes"]
	valueSortLine = ("gff_seqorigin", "gff_type", "gff_start", "gff_end", "gff_score", "gff_strand", "gff_phase", "gff_attributes","org_id", "gff_name","gff_protein_id","gff_transcript_id")
	insertQuery = "REPLACE INTO gff(gff_seqorigin, gff_type, gff_start, gff_end, gff_score, gff_strand, gff_phase, gff_attributes,org_id, gff_name,gff_protein_id,gff_transcript_id) values (%s,%s,%s,%s, %s,%s,%s,%s,%s,%s,%s,%s);"
	return load_tab_files(name, nrtab, valueLine, valueSortLine, insertQuery, org_name,filepath, action, dbname)

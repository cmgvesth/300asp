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
from utilsDataFormats import *
from utilsFileOperations import *

#######################################################################
# Fetch or create organism _id, auto increment 
#######################################################################
def get_org_id(org_name, dbname, test):
	# test of connection has to be performed prior ------------------------------------
	db = mdb.connect("localhost","asp","1234",dbname) 
	cursor = db.cursor()

	cursor.execute("SELECT org_id from organism where name LIKE \'" + org_name + "%\';")
	# first element of tuple ------------------------------------
	org_id = cursor.fetchone() 

	if not org_id:
		try:
			cursor.execute("REPLACE INTO organism(name) VALUES(\'" + org_name + "\');")
			db.commit()
			org_id = cursor.lastrowid
			
		except mdb.Error, e:
			sys.exit("# ERROR %d: %s" % (e.args[0],e.args[1]))

		# If function is being tested, delete the entry again
		if test == 1: 
			cursor.execute("DELETE FROM organism WHERE name = (\'" + org_name + "\');")
			db.commit()
		db.close()
		return org_id

	else:
		db.close()
		return org_id[0]

		
##############################################################################################################################################
##############################################################################################################################################
####################################################################### DATA #################################################################
##############################################################################################################################################
##############################################################################################################################################

#######################################################################
# Handle FASTA data
# INFO: This function calls the load_fasta_files function
#######################################################################

def fasta (org_name, filepath, action, filetype, dbname, prefix):
	(totalcounter, counter) = (0,0)
	org_id = ''
	filename = filepath.split("/")[-1]
	
	# Find the file name and add prefix if specified  ------------------------------------
	if prefix: filename = prefix+filename

	valueSortLine = ("orgkey","seqkey","tailkey","seqname", "seq", "org_id", "orgseq")

	# Protein FASTA ------------------------------------
	if filetype == "pf":
		insertQuery = "REPLACE INTO proteins(prot_orgkey, prot_seqkey, prot_tailkey, prot_seqname, prot_seq, org_id, prot_orgseq) values(%s,%s,%s,%s,%s,%s,%s);"
		
	# DNA FASTA ------------------------------------
	elif filetype == "tf":
		insertQuery = "REPLACE INTO transcripts(trans_orgkey, trans_seqkey, trans_tailkey, trans_seqname, trans_seq, org_id, trans_orgseq) values(%s,%s,%s,%s,%s,%s,%s);"

	elif filetype == "cf":
		insertQuery = "REPLACE INTO cds(cds_orgkey, cds_seqkey, cds_tailkey, cds_seqname, cds_seq, org_id, cds_orgseq) values(%s,%s,%s,%s,%s,%s,%s);"

	elif filetype == "uma" or filetype == "ma":
		insertQuery = "REPLACE INTO assembly(assembly_orgkey, assembly_seqkey, assembly_tailkey, assembly_seqname, assembly_seq, org_id, assembly_orgseq) values(%s,%s,%s,%s,%s,%s,%s);"
		
	else:
		sys.exit("# ERROR: filetype not recognized %s" % filename)
	
	return load_fasta_files(filetype, valueSortLine, insertQuery, org_name, filepath, action, dbname)	


#######################################################################
# Load FASTA type files
#######################################################################
def load_fasta_files(filetype, valueSortLine, insertQuery, org_name,filepath, action, dbname):
	filename = filepath.split("/")[-1]
	print "# INFO Processing %s file with action %s: %s " % (filetype, action, filename)

	# Init counters ------------------------------------
	(totalcounter, counter) = (0,0)

	# Init data structures ------------------------------------
	values_to_insert = []
	valuehash = {}

	# Unzip file and store in list ------------------------------------
	records = gzip_fasta_file(filepath) 
		
	# If this is empty, there was nothing that matched a FASTA format in the file ------------------------------------
	if len(list(records)) < 1:	
		sys.exit("# ERROR: Not FASTA format or empty FASTA entry")


	# Only get the org_id once and not for all records ------------------------------------
	org_id = get_org_id(org_name, dbname, 0) 
	
	# If this is empty, there was nothing that matched a FASTA format in the file	 ------------------------------------
	if len(list(records)) < 1:	sys.exit("# ERROR: Not FASTA format or empty FASTA entry")

	for r in records:
		counter += 1
		totalcounter += 1
		values = ()

		# Store variable names in dict and get values from record ------------------------------------
		if filetype == "pf":
			# Make "sure" sequence is amino acids ------------------------------------ 
			check_seq_protein(r.seq)

		if filetype == "ma" or filetype == "tp" or filetype == "cp" or filetype == "uma": 
			# Make "sure" sequence is DNA ------------------------------------ 
			check_seq_dna(r.seq)

		valuehash["seqname"]	= (r.description,)
		valuehash["seq"]		= (str(r.seq.lower()),)

		if len(r.description.split(" ")[0].split("|"))>3:
			valuehash["orgkey"]		= (r.description.split(" ")[0].split("|")[1],)
			valuehash["seqkey"]		= (r.description.split(" ")[0].split("|")[2],)
			valuehash["tailkey"]	= (r.description.split(" ")[0].split("|")[3],)
			# EXAMPLE: jgi|Aspac1|10177|gw1.1.1170.1 ------------------------------------
		else:
			valuehash["orgkey"]		= (org_name,)
			valuehash["seqkey"]		= (r.description.split(" ")[0],)
			valuehash["tailkey"]	= ('',)

		valuehash["orgseq"]		= (str(org_id)+"_"+valuehash["seqkey"][0],)
				
		# create list of values to be inserted ------------------------------------
		for val in valueSortLine:	
			# org_id is NOT FOUND IN FILE BUT OBTAINED FROM DB ------------------------------------
			if val == "org_id":	values += (org_id,)	
			else: 				values += valuehash[val]
			
		# Create sets of lists of values ------------------------------------
		values_to_insert.append(values) 
		
		# MySQL will only allow a certain data insert size ------------------------------------
		if filetype == "ma" : counter_threshold = 2 
		counter_threshold = 1000
		
		if action == "load" and ( counter == counter_threshold or totalcounter == len(records)):
			if totalcounter % 3000 == 0  or totalcounter == len(records) : print "# INFO: Inserting record number %s" % totalcounter

			try:	
				db = mdb.connect("localhost","asp","1234",dbname)
				cursor = db.cursor()
				cursor.executemany(insertQuery,values_to_insert)
				
				# add the changes to the db  ------------------------------------
				db.commit() 

				# dont leave the db connection open  ------------------------------------
				db.close()	
				
				# reset counter and values to insert ------------------------------------
				counter = 0 
				values_to_insert = []

			except mdb.Error, e:
				sys.exit( "# ERROR utils %s %d: %s" % (filetype, e.args[0],e.args[1]))

	# If user is only testing or has chosen verbose mode - print example variable values ------------------------------------			
	if action == "test":
		print "# Query: %s\n# Values: %s" % (insertQuery,values_to_insert[0])		

	# Make sure output is the right format ------------------------------------	
	if not (isinstance( org_id, long ) or  isinstance( org_id, int )):
		sys.exit("# ERROR: function fasta does not return an integer number")

	print "# FINISHED %s file: %s with total number of records %s\n" % (filetype, filename, totalcounter)

	return org_id

##############################################################################################################################################
##############################################################################################################################################
####################################################################### ANNOTATIONS ##########################################################
##############################################################################################################################################
##############################################################################################################################################

#######################################################################
# Parse lines
#######################################################################
def parse_annotation_line(line):
	segments = line.replace("\t", " ").split(" ")
	
	# Pull the ID value from the first segment
	protein_id = str(segments[0])	
	
	# Init data structures  ------------------------------------
	id_list = []
	term_list = []
	type_list = []
	acc_list = []
  	protein_list = []

	# Init indexes  ------------------------------------
	i = 1
	j = len(segments) - 1

	# Read IDs Note: i put i<j here to force this to terminate at some point   ------------------------------------
	# Move 'i' forward because of break
	while(i<j): 
		if(segments[i] != "|"):
			id_list.append(int(segments[i]))
			if(segments[i+1] != "|"):
				break
		i += 1
	i += 1	

	# Read GO accession numbers  ------------------------------------
	# Move 'j' backward because of break
	while(i<j):
		if(segments[j] != "|"):
			acc_list.append(segments[j])
			if(segments[j-1] != "|"):
				break
		j-= 1
	j -= 1	

	# Read GO types  ------------------------------------
	# No need to move j at this point
	while(i<j):	
		if(segments[j] != "|"):
			type_list.append(segments[j])
			protein_list.append(protein_id)
			if(segments[j-1] != "|"):
				break
		j -= 1	
	

	# Read go terms  ------------------------------------
	remaining = " ".join(segments[i:j])
	term_list = remaining.split(" | ")

	# Check Results  ------------------------------------
	l1 = len(id_list)
	if(l1 == 0):
		sys.exit("# ERROR: annotation parser ran into an error parsing the line (is the line is empty?)")
	if(l1 != len(term_list) or l1 != len(type_list) or l1 != len(acc_list)):
		sys.exit("# ERROR: annotation parser ran into an error parsing the line, there must be the same number of terms in each line")

	# Build and return result set
	# valueLine1 = ["goterm_id", "goName", "gotermType", "goAcc"]
	# valueLine2 = ["protein_id", "goterm_id", "org_id"]

	#return [(protein_id, ann_id, ann_term, ann_type, ann_acc) for (ann_id, ann_term, ann_type, ann_acc) in zip(id_list, term_list, reversed(type_list), reversed(acc_list))]
	return zip(protein_list, id_list, term_list, reversed(type_list), reversed(acc_list))
#######################################################################
# Load tab type files
#######################################################################
def load_tab_files(name, nrtab, valueLine, insertQuery, org_name,filepath, action, dbname):
	# Init counters ------------------------------------
	counter = 0
	totalcounter = 0

	# Init data structures ------------------------------------
	values_to_insert = []
	valuehash = {}
	filename = filepath.split("/")[-1]
	
	if action == "load":	print "# LOADING %s file: %s" % (name, filename)
	else:					print "# TESTING %s file: %s" % (name, filename)

	# function to unzip file and store in list ------------------------------------
	records = gzip_tab_file(filepath) 

	# function to verify that the tab file has the right number of fields ------------------------------------
	check_tab_format (nrtab, records) 
	
	# Only get the org_id once and not for all records ------------------------------------
	org_id = get_org_id(org_name, dbname, 0) 

	for r in records:
		counter += 1
		totalcounter += 1
		nr_values = len(valueLine)


		if name == "GO":
			annotations = parse_annotation_line(r)
			
			for valueset in annotations:
				valuehash[valueLine[0]] = (valueset[0],)
				valuehash[valueLine[1]] = (org_id,)

				for var in range(1,len(valueLine)):	
					valuehash[valueLine[var]] = (valueset[var-1],)
		else:
			# Store variable names in dict and get values from record ------------------------------------
			for var in range(0,len(valueLine)):	
				valuehash[valueLine[var]] = (r.split("\t")[var],)
		


		if action == "load" and ( counter == 1000 or totalcounter == len(records)):
			# Print insert information  ------------------------------------
			if totalcounter % 10000 == 0 : print "# INFO: Inserting record number %s" % totalcounter
			elif totalcounter % 2000 == 0: print "# INFO: Inserting record number %s" % totalcounter
			elif totalcounter == len(records): print "# INFO: Inserting record number %s" % totalcounter
						
			try:	
				db = mdb.connect("localhost","asp","1234",dbname)
				cursor = db.cursor()
				cursor.executemany(insertQuery,values_to_insert)

				# add the changes to the db  ------------------------------------
				db.commit() 

				# dont leave the db connection open  ------------------------------------
				db.close()	
				
				# reset counter and values to insert ------------------------------------
				counter = 0 
				values_to_insert = []

			except mdb.Error, e:
				sys.exit( "# ERROR utils %s %d: %s" % (name, e.args[0],e.args[1]))

	if action == "test" and ( counter == 1000 or totalcounter == len(records)):
		print "# INFO: query %s value %s" % ( insertQuery, values_to_insert )
						
	print "# FINISHED %s file: %s with total number of records %s" % (name,filename, totalcounter)

	if not (isinstance( org_id, long ) or  isinstance( org_id, int )):
		sys.exit("# ERROR: load_tab_files does not return an integer number")

	return org_id	
	
#######################################################################
# GO
#######################################################################
def go (org_name,filepath, action, dbname):	
	name = "GO"
	nrtab = 5
	valueLine = ["org_id", "protein_id", "goterm_id", "goName", "gotermType", "goAcc"]

	insertQuery1 = "REPLACE INTO go (go_term_id, go_name, go_termtype, go_acc) values(%s,%s,%s,%s);"
	insertQuery2 = "REPLACE INTO protein_has_go (go_protein_id, go_term_id, org_id) values(%s,%s,%s);"

	load1 = load_tab_files(name, nrtab, valueLine, insertQuery1, org_name,filepath, action, dbname)
	load2 = load_tab_files(name, nrtab, valueLine, insertQuery2, org_name,filepath, action, dbname)

	return load1, load2

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

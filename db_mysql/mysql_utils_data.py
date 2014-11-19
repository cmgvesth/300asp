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
# Returns integer org_id 
#######################################################################
def get_org_id(org_name, dbname, test):
	# test of connection has to be performed prior ------------------------------------
	db = mdb.connect("localhost","asp","1234",dbname) 
	cursor = db.cursor()

	cursor.execute("SELECT org_id from organism where name LIKE \'" + org_name + "%\';")
	# first element of tuple ------------------------------------
	result = cursor.fetchone() 
	org_id = result[0]

	if not org_id:
		try:
			cursor.execute("REPLACE INTO organism(name) VALUES(\'" + org_name + "\');")
			db.commit()
			result = cursor.lastrowid
			org_id = result[0]

		except mdb.Error, e:
			sys.exit("# ERROR %d: %s" % (e.args[0],e.args[1]))

		# If function is being tested, delete the entry again
		if test == 1: 
			cursor.execute("DELETE FROM organism WHERE name = (\'" + org_name + "\');")
			db.commit()

	if not (isinstance( org_id, long ) or  isinstance( org_id, int )):
		sys.exit("# ERROR: get_org_id does not return an integer number")

	db.close()
	return org_id

#######################################################################
# TEST function - get_org_id
# As the org_id is auto incremented in the db it is not possible to test for an exact value
# The best we can do is to test that it is the right format
#######################################################################
def test__get_org_id():
	(org_name, dbname, prot_id) = ("Aspni1", "aspminedb", "19949")
	# Test get_org_id
	org_id = get_org_id(org_name, dbname, 1)
		
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
# TEST function - fasta
# This test makes sure the fucntions does not error on the test input
# It does not test the actual outcome of the function
#######################################################################
def test__fasta ():
	(org_name, dbname, prot_id, prefix) = ("Aspni1", "aspminedb", "19949", "")

	# Create temporary dna FASTA file
	fasta_dna_entry = ">scaffold_1\nTTTTATAAAAA\nAAACCAACtagttta\nttaatACtta\n>scaffold_7\ntaatgtaaaccta\natgtactaaTGT\nAATACCTCC"
	#fasta_dna_entry = ">jgi|Aspac1|48858|fgenesh1_pm.1_#_1\nMSPAFRIQAISIPNNYLPPRKSFTPLLSPLSPSLRLFGRYSLHRHLLVELGEWEAAYEKASAFVADLTTA\nRGFERVEIDSGGHSYVIFDLRRRDLSYWDVMAQEWLVAPGSYNVFIGASSRDLRLNGTFSLTVTA*\n>jgi|Aspac1|37826|Genemark1.2_g\nMDCDETKPECSNCVNHSVRCIYDSPAPKTKKPVTLPAQTHKTENKVQVQFVEFNFLHVQPAPAPVPATPA\nRLEQPFWEALQRKQPVALMLVAHFAALMRWMDGTWVSVGWPEHILRGAWMFLEEEERGLVRWLIDDVGVL\nGLSG*\n>jgi|Aspac1|37827|Genemark1.3_g\nMVRGIIFQLVFVGLVVDFVLRLSKRGRPQALWSNRPLLLLGGATALSLLLIYIRSVYRTIELLHGWTSST\nMHNEMLLIGLDGAIMVPAFSVYNLLHPGYLLPKVQREVGYLDARGLQMEMEAVKDGRYQIVEVEGGKDDS\nVTLLDRS*"
	tmp_dna_filename = 'tmp_' + org_name + '_dna.fsa.gz'
	tmp_dna_file = gzip.open(tmp_dna_filename, 'wb') 
	tmp_dna_file.write(fasta_dna_entry)
	tmp_dna_file.close()

	# Create temporary gene dna FASTA file
	fasta_gene_entry = ">jgi|Aspac1|48858|fgenesh1_pm.1_#_1\nATGTCCCTTCA\nATCTCTGG\nCGCGACAGCCTAA\n>jgi|Aspac1|37826|Genemark1.2_g\nATGGATT\nTGCGGTTG\nGGGTTGA\n>jgi|Aspac1|37827|Genemark1.3_g\nATGGTCCCA\nGTTAAGA\n"
	#fasta_dna_entry = ">jgi|Aspac1|48858|fgenesh1_pm.1_#_1\nMSPAFRIQAISIPNNYLPPRKSFTPLLSPLSPSLRLFGRYSLHRHLLVELGEWEAAYEKASAFVADLTTA\nRGFERVEIDSGGHSYVIFDLRRRDLSYWDVMAQEWLVAPGSYNVFIGASSRDLRLNGTFSLTVTA*\n>jgi|Aspac1|37826|Genemark1.2_g\nMDCDETKPECSNCVNHSVRCIYDSPAPKTKKPVTLPAQTHKTENKVQVQFVEFNFLHVQPAPAPVPATPA\nRLEQPFWEALQRKQPVALMLVAHFAALMRWMDGTWVSVGWPEHILRGAWMFLEEEERGLVRWLIDDVGVL\nGLSG*\n>jgi|Aspac1|37827|Genemark1.3_g\nMVRGIIFQLVFVGLVVDFVLRLSKRGRPQALWSNRPLLLLGGATALSLLLIYIRSVYRTIELLHGWTSST\nMHNEMLLIGLDGAIMVPAFSVYNLLHPGYLLPKVQREVGYLDARGLQMEMEAVKDGRYQIVEVEGGKDDS\nVTLLDRS*"
	tmp_gene_filename = 'tmp_' + org_name + '_gene.fsa.gz'
	tmp_gene_file = gzip.open(tmp_gene_filename, 'wb') 
	tmp_gene_file.write(fasta_gene_entry)
	tmp_gene_file.close()

	# Create temporary protein FASTA file
	fasta_protein_entry = ">jgi|Aspac1|43358|fgenesh1_pm.1_#_1\nMSPAFRIQDLTTA\nRGFGGLTVTA*\n>jgi|Aspac1|37826|Genemark1.2_g\nMDCDETKP\nRLEDVGVL\nGLSG*\n>jgi|Aspac1|37827|Genemark1.3_g\nMVRGITSST\nMHNEMLS\nVTLLDRS*"
	#fasta_protein_entry = ">scaffold_1\nTTTTTAATTTTCGAACTAACTAAAAAAACTTTTATTTAATCTGATAAAAA\nAAACCAACtattattatattaaattctaattatatattttatttagttta\nttaattttttaataatattaCAGAGAAGAGCAGTAATTACTGTAGACtta\n>scaffold_7\ntaatgtaaaccctaatgtaaaccctaatgtaaaccctaatgtaaacccta\natgtaaaccctaatgtaaccctaacccccctaatctgacaaccctaaTGT\nAATACCTTCGGCGGAAGGTTTATCCGGTGCTGACGGGCCGAGGTTTCTCC"
	tmp_protein_filename = 'tmp_' + org_name + '_protein.fsa.gz'
	tmp_protein_file = gzip.open(tmp_protein_filename, 'wb') 
	tmp_protein_file.write(fasta_protein_entry)
	tmp_protein_file.close()

	# Test different versions of FASTA function
	action = "test"
	(org_id, insert_values) = fasta(org_name, tmp_protein_filename, action, "pf", dbname, prefix)
	eval = cmp( insert_values, ['Aspac1', '43358', 'fgenesh1_pm.1_#_1', 'jgi|Aspac1|43358|fgenesh1_pm.1_#_1', 'mspafriqdlttargfggltvta*', "13762L", '13762_43358'])

	org_id = fasta(org_name, tmp_gene_filename, action, "tf", dbname, prefix)
	eval = cmp( insert_values, ['Aspac1', '48858', 'fgenesh1_pm.1_#_1', 'jgi|Aspac1|48858|fgenesh1_pm.1_#_1', 'atgtcccttcaatctctggcgcgacagcctaa', "13763L", '13763_48858'])

	org_id = fasta(org_name, tmp_gene_filename, action, "cf", dbname, prefix)
	eval = cmp( insert_values, ['Aspac1', '48858', 'fgenesh1_pm.1_#_1', 'jgi|Aspac1|48858|fgenesh1_pm.1_#_1', 'atgtcccttcaatctctggcgcgacagcctaa', "13764L", '13764_48858'])

	org_id = fasta(org_name, tmp_dna_filename, action, "ma", dbname, prefix)
	eval = cmp( insert_values, ['Aspni1', 'scaffold_1', '', 'scaffold_1', 'ttttataaaaaaaaccaactagtttattaatactta', "13765L", '13765_scaffold_1'])


#######################################################################
# Load FASTA type files
#######################################################################
def load_fasta_files(filetype, valueSortLine, insertQuery, org_name,filepath, action, dbname):
	filename = filepath.split("/")[-1]
	print "# INFO Processing %s file with action %s: %s " % (filetype, action, filename)

	# Init counters ------------------------------------
	(totalcounter, counter) = (0,0)

	# Init data structures ------------------------------------
	insert_values = []
	valuedict = {}

	# Unzip file and store in list ------------------------------------
	records = gzip_fasta_file(filepath) 
		
	# Only get the org_id once and not for all records ------------------------------------
	if action == "test": org_id = get_org_id(org_name, dbname, 1) 
	else: org_id = get_org_id(org_name, dbname, 0) 
	
	for r in records:
		counter += 1
		totalcounter += 1
		values = ()

		# Store variable names in dict and get values from record ------------------------------------
		# Make "sure" sequence is amino acids or DNA depemnding in filetype ------------------------------------ 
		if filetype == "pf":	check_seq_protein(r.seq)
		else:					check_seq_dna(r.seq)

		valuedict["seqname"]	= (r.description,)
		valuedict["seq"]		= (str(r.seq.lower()),)

		if len(r.description.split(" ")[0].split("|"))>3:
			valuedict["orgkey"]		= (r.description.split(" ")[0].split("|")[1],)
			valuedict["seqkey"]		= (r.description.split(" ")[0].split("|")[2],)
			valuedict["tailkey"]	= (r.description.split(" ")[0].split("|")[3],)
			# EXAMPLE: jgi|Aspac1|10177|gw1.1.1170.1 ------------------------------------
		else:
			valuedict["orgkey"]		= (org_name,)
			valuedict["seqkey"]		= (r.description.split(" ")[0],)
			valuedict["tailkey"]	= ('',)

		valuedict["orgseq"]		= (str(org_id)+"_"+valuedict["seqkey"][0],)
				
		# create list of values to be inserted ------------------------------------
		for val in valueSortLine:	
			# org_id is NOT FOUND IN FILE BUT OBTAINED FROM DB ------------------------------------
			if val == "org_id":	values += (org_id,)	
			else: 				values += valuedict[val]
			
		# Create sets of lists of values ------------------------------------
		insert_values.append(values) 
		
		# MySQL will only allow a certain data insert size ------------------------------------
		if filetype == "ma" : counter_threshold = 2 
		counter_threshold = 1000
		
		if action == "load" and ( counter == counter_threshold or totalcounter == len(records)):
			if totalcounter % 3000 == 0  or totalcounter == len(records) : print "# INFO: Inserting record number %s" % totalcounter

			try:	
				db = mdb.connect("localhost","asp","1234",dbname)
				cursor = db.cursor()
				cursor.executemany(insertQuery,insert_values)
				
				# add the changes to the db, close the db connection 
				db.commit() 
				db.close()	
				
				# reset counter and values to insert ------------------------------------
				counter = 0 
				insert_values = []

			except mdb.Error, e:
				sys.exit( "# ERROR utils %s %d: %s" % (filetype, e.args[0],e.args[1]))


	# Make sure output is the right format ------------------------------------	
	if not (isinstance( org_id, long ) or  isinstance( org_id, int )):
		sys.exit("# ERROR: function fasta does not return an integer number")

	print "# FINISHED %s file: %s with total number of records %s\n" % (filetype, filename, totalcounter)

	# If user is only testing or has chosen verbose mode - print example variable values ------------------------------------			
	if action == "test":
		print "# Query: %s\n# Values: %s" % (insertQuery,insert_values[0])		

	return (org_id, str(insert_values[0]))

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
def load_tab_files(name, nrtab, valueLine, main_insertQuery, connect_insertQuery, org_name,filepath, action, dbname):
	# Init counters ------------------------------------
	(totalcounter, counter) = (0,0)

	# Init data structures ------------------------------------
	main_insert_values = []
	connect_insert_values = []
	valuedict = {}
	filename = filepath.split("/")[-1]
	
	if action == "load":	print "# LOADING %s file: %s" % (name, filename)
	else:					print "# TESTING %s file: %s" % (name, filename)

	# function to unzip file and store in list ------------------------------------
	records = gzip_tab_file(filepath) 

	# function to verify that the tab file has the right number of fields ------------------------------------
	check_tab_format (nrtab, records) 
	
	# Only get the org_id once and not for all records ------------------------------------
	if action == "test": org_id = get_org_id(org_name, dbname, 1) 
	else: org_id = get_org_id(org_name, dbname, 0) 

	for r in records:
		counter += 1
		totalcounter += 1
		

		# Function to deal with multi value tab fields (protein_id 	valueA | valueB)
		annotations = parse_annotation_line(r)
		
		for annset in annotations:
			# First element of valueline is the org id, the secund is the protein id
			valuedict[valueLine[0]] = (org_id,)
			valuedict[valueLine[1]] = (annset[0],)

			# Store the remaing values in the valuedict
			for var in range(1,len(valueLine)):	
				valuedict[valueLine[var]] = (annset[var-1],)

			# Main annotation table - not organism specific, no org or protein information	
			# Clear the values list and sort valuedict by the given value line
			# Add each value to a list and add each list to the insert list
			values = ()
			for value in [valuedict[k] for k in valueLine[2:]]:
				values += value	
			main_insert_values.append(values)

			# Connection annotation table - only org, protein and annotation ids, no annotation information	
			# Clear the values list and sort valuedict by the given value line
			# Add each value to a list and add each list to the insert list
			values = ()
			for value in [valuedict[k] for k in valueLine[0:3]]:
				values += value	
			connect_insert_values.append(values)

		if action == "load" and ( counter == 1000 or totalcounter == len(records)):
			# Print insert information  ------------------------------------
			if totalcounter % 10000 == 0 : print "# INFO: Inserting record number %s" % totalcounter
			elif totalcounter % 2000 == 0: print "# INFO: Inserting record number %s" % totalcounter
			elif totalcounter == len(records): print "# INFO: Inserting record number %s" % totalcounter
						
			try:	
				db = mdb.connect("localhost","asp","1234",dbname)
				cursor = db.cursor()
				cursor.executemany(main_insertQuery,	main_insert_values)
				cursor.executemany(connect_insertQuery,	connect_insert_values)

				# add the changes to the db, close the db connection 
				db.commit() 
				db.close()	
				
				# reset counter and values to insert ------------------------------------
				counter = 0 
				insert_values = []

			except mdb.Error, e:
				sys.exit( "# ERROR utils %s %d: %s" % (name, e.args[0],e.args[1]))

	# If user is only testing or has chosen verbose mode - print example variable values ------------------------------------			
	if action == "test":
		print "# Query: %s\n# Values: %s" % (insertQuery,insert_values)		

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
	valueLine = ["org_id", "protein_id", "go_term_id", "go_name", "go_termtype", "go_acc"]

	main_insertQuery = "REPLACE INTO go (go_term_id, go_name, go_termtype, go_acc) values(%s,%s,%s,%s);"
	connect_insertQuery = "REPLACE INTO protein_has_go (org_id, protein_id, go_term_id ) values(%s,%s,%s);"

	org_id = load_tab_files(name, nrtab, valueLine, main_insertQuery, connect_insertQuery, org_name,filepath, action, dbname)
	#load_tab_files(name, nrtab, valueLine, insertQuery2, org_name,filepath, action, dbname)

	return org_id

def test__go():
	(org_name, dbname, prot_id, prefix) = ("Aspni1", "aspminedb", "19949", "")

	# Create temporary dna FASTA file
	tab_entry = ""
	#fasta_dna_entry = ">jgi|Aspac1|48858|fgenesh1_pm.1_#_1\nMSPAFRIQAISIPNNYLPPRKSFTPLLSPLSPSLRLFGRYSLHRHLLVELGEWEAAYEKASAFVADLTTA\nRGFERVEIDSGGHSYVIFDLRRRDLSYWDVMAQEWLVAPGSYNVFIGASSRDLRLNGTFSLTVTA*\n>jgi|Aspac1|37826|Genemark1.2_g\nMDCDETKPECSNCVNHSVRCIYDSPAPKTKKPVTLPAQTHKTENKVQVQFVEFNFLHVQPAPAPVPATPA\nRLEQPFWEALQRKQPVALMLVAHFAALMRWMDGTWVSVGWPEHILRGAWMFLEEEERGLVRWLIDDVGVL\nGLSG*\n>jgi|Aspac1|37827|Genemark1.3_g\nMVRGIIFQLVFVGLVVDFVLRLSKRGRPQALWSNRPLLLLGGATALSLLLIYIRSVYRTIELLHGWTSST\nMHNEMLLIGLDGAIMVPAFSVYNLLHPGYLLPKVQREVGYLDARGLQMEMEAVKDGRYQIVEVEGGKDDS\nVTLLDRS*"
	tmp_tab_filename = 'tmp_' + org_name + '_tab.gz'
	tmp_tab_file = gzip.open(tmp_tab_filename, 'wb') 
	tmp_tab_file.write(fasta_tab_entry)
	tmp_tab_file.close()

	# Test different versions of FASTA function
	action = "test"
	org_id = go (org_name, path, action, dbname)

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

def test__sigp():
	"""
	org_id = sigp 	(org_name, path, action, dbname)
	"""

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
def test__ipr():
	"""
	org_id = ipr 	(org_name, path, action, dbname)
	"""

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
def test__kegg():
	"""
	org_id = kegg 	(org_name, path, action, dbname)
	"""

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
def test__kog():
	"""
	org_id = kog 	(org_name, path, action, dbname)
	"""

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
def test__gff():
	"""
	org_id = gff 	(org_name, path, action, dbname)
	"""

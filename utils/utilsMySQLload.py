#!/usr/bin/python
'''##################################################################
# Imports
##################################################################'''
import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), "../utils/"))
from aspmine_imports import *

#######################################################################
# Fetch or create organism _id, auto increment
# Returns integer org_id 
#######################################################################
def get_org_id(org_name, dbname, test):
	# test of connection has to be performed prior ------------------------------------
	db = mdb.connect("localhost","asp","1234",dbname) 
	cursor = db.cursor()
	result = []
	org_id = ""

	# the limit at the event ensures that only one row is selected
	cursor.execute("SELECT org_id from organism where name LIKE \'" + org_name + "%\' limit 1;")
	result = cursor.fetchone() 
	if result:
		org_id = result[0]

	if not org_id:
		try:
			cursor.execute("REPLACE INTO organism(name) VALUES(\'" + org_name + "\');")
			db.commit()
			# this function always returns one single value
			result = cursor.lastrowid
			if result:
				org_id = result

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
	example_values = []

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
		example_values = [ str(k)[:10] for k in values ]
		
		# MySQL will only allow a certain data insert size ------------------------------------
		counter_threshold = 1000
		if filetype == "ma" : counter_threshold = 2 
		
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
		print "# Query: %s\n# Values (substrings): %s" % (insertQuery, str(example_values))

	return (org_id, str(example_values))

##############################################################################################################################################
##############################################################################################################################################
####################################################################### ANNOTATIONS ##########################################################
##############################################################################################################################################
##############################################################################################################################################
#######################################################################
# Parse generic annotation lines
#######################################################################
def parse_generic_annotation_line(line):
	annotations = []
	# pattern for KEGG files
	if re.search(r"\\N", line):
		annotations = line.replace("\N", "NULL").split("\t")
	else:	
		annotations = line.split("\t")
	if len(annotations) < 1:
		sys.exit("# ERROR: annotation line is empty, %s" % line)	
	return annotations

#######################################################################
# Parse GFF annotation lines
#######################################################################
def parse_gff_annotation_line(line):
	annotations = []
	(seq_id, exon_id) = ("NULL","NULL")
	# pattern for KEGG files
	line = line.replace("\"", "")
	annotations = line.split("\t")

	if len(annotations) < 1:
		sys.exit("# ERROR: annotation line is empty, %s" % line)	

	k = annotations[-1]
	if re.search("proteinId", k):
		seq_id = (k.split(";")[1]).split(" ")[-1]
		exon_id = (k.split(";")[2]).split(" ")[-1]

	if re.search("transcriptId", k):
		seq_id = (k.split(";")[1]).split(" ")[-1]

	annotations.append(seq_id)
	annotations.append(exon_id)	
	annotations.append((k.split(";")[0]).split(" ")[-1])
	#print annotations
	return annotations



#######################################################################
# Parse InterPro annotation lines
#######################################################################
def parse_IPR_annotation_line(line):
	annotations = []
	tab_fields = line.split("\t")
	tab_fields = [t.rstrip(",") for t in tab_fields]
	tab_fields = [t.replace(";", "") for t in tab_fields]

	# If the domain count fieled is more than 1, process line accordingly
	if int(tab_fields[6]) > 1:
		starts_list	= tab_fields[7].split(",")
		ends_list	= tab_fields[8].split(",")
		score_list	= tab_fields[9].split(",")

		# If the number of values is not the same something went wrong
		if len(starts_list) != len(ends_list) or len(starts_list) != len(score_list):
			sys.exit("# ERROR: number of elements in starts, ends or scores is not the same %s" % str(tab_fields))

		for e in range(0,len(starts_list)):
			annotations.append( tuple(tab_fields[0:7] + [ starts_list[e] , ends_list[e] , score_list[e] ]) )	
	# Add value line as annotation		
	else:
		annotations.append( tuple(tab_fields[0:]) )
	#annotations = [ k.replace(";", "") for k in annotations]	
	return annotations
			
#######################################################################
# Parse GO annotation lines
#######################################################################
def parse_GO_annotation_line(line):
	annotations = []
	tab_fields = line.split("\t")
	tab_fields = [t.replace(" | ", "|") for t in tab_fields]

	# If there is more than one value in the first field, process line accordingly
	if len(tab_fields[1].split("|")) > 1:
		term_list	= tab_fields[1].split("|")
		name_list	= tab_fields[2].split("|")
		type_list	= tab_fields[3].split("|")
		acc_list	= tab_fields[4].split("|")

		# If the number of values is not the same something went wrong
		if len(term_list) != len(name_list) or len(term_list) != len(type_list):
			sys.exit("# ERROR: number of elements in terms, name or type is not the same %s" % str(tab_fields))

		for e in range(0,len(term_list)):
			annotations.append( tuple([tab_fields[0]] + [ term_list[e] , name_list[e] , type_list[e] , acc_list[e]]) )	
	# Add value line as annotation		
	else:
		annotations.append( tuple(tab_fields[0:]) )
	return annotations

#######################################################################
# Load tab type files
#######################################################################
def load_tab_files(name, nrtab, main_values, connect_values, line_values, main_insertQuery, connect_insertQuery, org_name,filepath, action, dbname):
	# Init counters ------------------------------------
	(totalcounter, counter, print_flag) = (0,0,0)

	# Init data structures ------------------------------------
	main_insert_values = []
	connect_insert_values = []
	annotations = []
	valuedict = {}
	filename = filepath.split("/")[-1]
	
	if action == "load":	print "# LOADING %s file: %s" % (name, filename)
	else:					print "# TESTING %s file: %s" % (name, filename)

	# function to unzip file and store in list ------------------------------------
	records = gzip_tab_file(filepath) 

	# function to verify that the tab file has the right number of fields ------------------------------------
	check_tab_format (nrtab, records) 
	
	# Only get the org_id once and not for all records ------------------------------------
	if action == "test": 
		org_id = get_org_id(org_name, dbname, 1) 
		print "# INFO: number of annotation records %s for organism id %s" % (len(records), org_id)
		print "# INFO: example line: %s" % str(records[3])
	else: org_id = get_org_id(org_name, dbname, 0) 


	for r in records:
		counter += 1
		totalcounter += 1
		valuedict["org_id"] = (org_id,)

		# Function to deal with multi value tab fields (protein_id 	valueA | valueB)
		if 		name == "GO": 		annotations = parse_GO_annotation_line(r)
		elif 	name == "InterPro": annotations = parse_IPR_annotation_line(r)
		elif 	name == "GFF":		annotations = parse_gff_annotation_line(r)
		else:						annotations = parse_generic_annotation_line(r)

		if name == "GO" or name == "InterPro":
			for annotation in annotations:
				for key in line_values:
					valuedict[key] = (annotation[line_values.index(key)],) 
		elif name == "GFF":

			for key in line_values:
				#print key, line_values.index(key)

				if re.search("proteinId", annotations[line_values.index(key)]):
					valuedict["gff_protein_id"] = (annotations[line_values.index(key)+1],)
					valuedict["gff_trans_id"] = ("NULL",)

				elif re.search("transcriptId", annotations[line_values.index(key)]):
					valuedict["gff_trans_id"] = (annotations[line_values.index(key)+1],)
					valuedict["gff_protein_id"] = ("NULL",)

				elif key=="gff_seq_id": 
					continue	

				valuedict[key] = (annotations[line_values.index(key)],) 
				#if re.search("transcriptId", line_values["gff_attributes"]): 
				#	valuedict["gff_trans_id"] = (annotations[line_values["gff_seq_id"]],)

		else:			
			for key in line_values:
				valuedict[key] = (annotations[line_values.index(key)],) 
		#print valuedict		
		#sys.exit()		
		# Main annotation table - not organism specific, no org or protein information	
		# Clear the values list and sort valuedict by the given value line
		# Add each value to a list and add each list to the insert list
		values = ()
		for value in [valuedict[k] for k in main_values]:
			values += value	
		main_insert_values.append(values)

		# Connection annotation table - only org, protein and annotation ids, no annotation information	
		# Clear the values list and sort valuedict by the given value line
		# Add each value to a list and add each list to the insert list
		values = ()
		for value in [valuedict[k] for k in connect_values]:
			values += value	
		connect_insert_values.append(values)

		if action == "load" and ( counter == 1000 or totalcounter == len(records)):
			# Print insert information  ------------------------------------
			if totalcounter % 10000 == 0 and print_flag == 1: 
				print "# INFO: Inserting record number %s" % totalcounter
			elif totalcounter % 2000 == 0 and print_flag == 0: 
				print_flag = 1
				print "# INFO: Inserting record number %s" % totalcounter
			elif totalcounter == len(records): print "# INFO: Inserting record number %s" % totalcounter
						
			try:	
				db = mdb.connect("localhost","asp","1234",dbname)
				cursor = db.cursor()
				cursor.executemany(main_insertQuery,	main_insert_values)
				if name != "GFF":
					cursor.executemany(connect_insertQuery,	connect_insert_values)

				# add the changes to the db, close the db connection 
				db.commit() 
				db.close()	
				
				# reset counter and values to insert ------------------------------------
				counter = 0 
				main_insert_values = []
				connect_insert_values = []

			except mdb.Error, e:
				print "# ERROR: %s" % str(main_insert_values[0])
				sys.exit( "# ERROR utils %s %d: %s" % (name, e.args[0],e.args[1]))

	# If user is only testing or has chosen verbose mode - print example variable values ------------------------------------			
	if action == "test":
		print "# Query: %s\n# Values: %s" % (main_insertQuery,main_insert_values[0])		
		print "# Query: %s\n# Values: %s" % (connect_insertQuery,connect_insert_values[0])		
		#print annotations
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
	
	line_values 	= ["protein_id", "go_term_id", "go_name", "go_termtype", "go_acc"]
	main_values 	= ["go_term_id", "go_name", "go_termtype", "go_acc"]
	connect_values 	= ["org_id", "protein_id", "go_term_id"]

	main_insertQuery 	= "REPLACE INTO go (go_term_id, go_name, go_termtype, go_acc) values(%s);" % ("%s," * len(main_values)).rstrip(",")
	connect_insertQuery = "REPLACE INTO protein_has_go (org_id, protein_id, go_term_id ) values(%s);" % ("%s," * len(connect_values)).rstrip(",")

	org_id = load_tab_files(name, nrtab, main_values, connect_values, line_values, \
		main_insertQuery, connect_insertQuery, org_name,filepath, action, dbname)

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
# InterPro
#######################################################################
def ipr (org_name,filepath, action, dbname):	
	name = "InterPro"
	nrtab = 10

	line_values 	= ["protein_id", "ipr_id", "desc", "domaindb", "domain_id", "domaindesc", "numHits", "domain_starts", "domain_ends", "score"]
	main_values 	= ["ipr_id", "desc", "domaindb", "domain_id", "domaindesc"]
	connect_values 	= ["org_id", "protein_id", "ipr_id", "domain_starts", "domain_ends", "score"]
	
	main_insertQuery 	= "REPLACE INTO ipr (ipr_id, ipr_desc, ipr_domaindb, ipr_domain_id, ipr_domaindesc) VALUES(%s);" % ("%s," * len(main_values)).rstrip(",")
	connect_insertQuery = "REPLACE INTO protein_has_ipr (org_id, protein_id, ipr_id, ipr_domain_start, ipr_domain_end, ipr_score ) values(%s);" % ("%s," * len(connect_values)).rstrip(",")

	org_id = load_tab_files(name, nrtab, main_values, connect_values, line_values, main_insertQuery, connect_insertQuery, org_name,filepath, action, dbname)
	return org_id

def test__ipr():
	"""
	org_id = ipr 	(org_name, path, action, dbname)
	"""

#######################################################################
# KEGG
#######################################################################
def kegg (org_name,filepath, action, dbname):	
	#proteinId	ecNum	definition	catalyticActivity	cofactors	associatedDiseases	pathway	pathway_class	pathway_type
	# the combination of ecNum and pathway is unique
	name = "KEGG"
	nrtab = 9

	line_values 	= ["protein_id", "ecNum", "definition", "catalyticActivity", "cofactors", "associatedDiseases", "pathway", "pathway_class", "pathway_type"]
	main_values 	= ["ecNum", "definition", "catalyticActivity", "cofactors", "associatedDiseases", "pathway", "pathway_class", "pathway_type"]
	connect_values 	= ["org_id", "protein_id", "ecNum", "pathway","ecNum", "pathway"]


	main_insertQuery = "REPLACE INTO kegg (kegg_ecNum, kegg_definition, kegg_catalyticActivity, kegg_cofactors, kegg_associatedDiseases, kegg_pathway, kegg_pathway_class, kegg_pathway_type) VALUES(%s);" % ("%s," * len(main_values)).rstrip(",")
	connect_insertQuery = "REPLACE INTO protein_has_kegg (org_id, protein_id, kegg_ecNum, kegg_pathway, kegg_id ) values(%s,%s,%s,%s, (SELECT kegg_id from kegg where kegg_ecNum=%s and kegg_pathway=%s));"

	org_id = load_tab_files(name, nrtab, main_values, connect_values, line_values, main_insertQuery, connect_insertQuery, org_name,filepath, action, dbname)
	return org_id

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

	line_values 	= ["transcript_id", "protein_id", "kog_id", "kogdefline", "class", "group"]
	main_values 	= ["kog_id", "kogdefline", "class", "group"]
	connect_values 	= ["org_id", "protein_id", "kog_id"]

	main_insertQuery = "REPLACE INTO kog (kog_id, kog_defline, kog_Class, kog_Group) VALUES(%s);" % ("%s," * len(main_values)).rstrip(",")
	connect_insertQuery = "REPLACE INTO protein_has_kog (org_id, protein_id, kog_id) VALUES(%s);" % ("%s," * len(connect_values)).rstrip(",")
	org_id = load_tab_files(name, nrtab, main_values, connect_values, line_values, main_insertQuery, connect_insertQuery, org_name,filepath, action, dbname)
	return org_id

def test__kog():
	"""
	org_id = kog 	(org_name, path, action, dbname)
	"""

#######################################################################
# SignalP
#######################################################################
def sigp (org_name,filepath, action, dbname):	
	name = "SignalP"
	nrtab = 5

	line_values 	= ["protein_id", "nn_cutpos", "neuro_net_vote", "hmm_cutpos", "hmm_signalpep_probability"]
	main_values 	= ["org_id", "protein_id", "nn_cutpos", "neuro_net_vote", "hmm_cutpos", "hmm_signalpep_probability"]
	connect_values 	= []
	
	main_insertQuery 	= "REPLACE INTO sigp (org_id, protein_id, sigp_nn_cutpos, sigp_neuro_net_vote, sigp_hmm_cutpos, sigp_hmm_signalpep_probability) VALUES(%s);" % ("%s," * len(main_values)).rstrip(",")
	connect_insertQuery = "SELECT sigp_nn_cutpos from sigp limit 1"

	org_id = load_tab_files(name, nrtab, main_values, connect_values, line_values, main_insertQuery, connect_insertQuery, org_name,filepath, action, dbname)
	return org_id

def test__sigp():
	"""
	org_id = sigp 	(org_name, path, action, dbname)
	"""

#######################################################################
# GFF
#######################################################################
def gff (org_name,filepath, action, dbname):	
	name = "GFF"
	nrtab = 9
	#valueLine = []
	#valueSortLine = ("gff_seqorigin", "gff_type", "gff_start", "gff_end", "gff_score", "gff_strand", "gff_phase", "gff_attributes","org_id", "gff_name","gff_protein_id","gff_transcript_id")
	#insertQuery = "REPLACE INTO gff(gff_seqorigin, gff_type, gff_start, gff_end, gff_score, gff_strand, gff_phase, gff_attributes,org_id, gff_name,gff_protein_id,gff_transcript_id) values (%s,%s,%s,%s, %s,%s,%s,%s,%s,%s,%s,%s);"
	#return load_tab_files(name, nrtab, valueLine, valueSortLine, insertQuery, org_name,filepath, action, dbname)

	line_values 	= ["gff_seqorigin", "a","gff_type", "gff_start", "gff_end", "gff_score", "gff_strand", "gff_phase", "gff_attributes", "gff_seq_id", "gff_exonnr", "gff_name"]
	main_values 	= ["gff_seqorigin", "gff_type", "gff_start", "gff_end", "gff_score", "gff_strand", "gff_phase", "gff_attributes", "org_id", "gff_protein_id", "gff_trans_id" ,"gff_exonnr","gff_name"]
	connect_values 	= []

	main_insertQuery = "REPLACE INTO gff(gff_seqorigin, gff_type, gff_start, gff_end, gff_score, gff_strand, gff_phase, gff_attributes,org_id, gff_protein_id, gff_trans_id, gff_exonnr, gff_name) VALUES(%s);" % ("%s," * len(main_values)).rstrip(",")
	connect_insertQuery = "SELECT gff_seqorigin from gff limit 1"
	org_id = load_tab_files(name, nrtab, main_values, connect_values, line_values, main_insertQuery, connect_insertQuery, org_name,filepath, action, dbname)
	return org_id


def test__gff():
	"""
	org_id = gff 	(org_name, path, action, dbname)
	"""

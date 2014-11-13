#!/usr/bin/python

import sys, os, re
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio import SeqIO
import MySQLdb as mdb


#------------------------------------------------------------------
# One class per data type
#------------------------------------------------------------------

#------------------------------------------------------------------
#------------------------------------------------------------------
class Fasta:
#------------------------------------------------------------------
#------------------------------------------------------------------

	#------------------------------------------------------------------
	# Initialize
	#------------------------------------------------------------------
	def __init__(self, org_name, entry, filename, filetype, dbname):
		self.entry = entry
		self.desc= "This is a FASTA object from JGI data"
		self.filetype = filetype
		self.filename = filename
		
		self.seq_name = self.entry.description
		id_line = self.seq_name.split(" ")[0]
		self.jgi = id_line.split("_")[-1]
		self.seq = self.entry.seq.lower()

		self.name = "_".join(self.seq_name.split("_")[:-1])

		self.warning = ''
		self.error = ''

		self.org_name = org_name

		if len(self.seq_name.split("|")) > 1:
			self.tpid = self.seq_name.split("|")[2]
		else:
			self.tpid = ''
		#------------------------------------------------------------------
		# Make "sure" sequence is DNA 
		#------------------------------------------------------------------
		if filetype == "uma" or filetype == "ma" or filetype == "tf" or filetype == "cf": 
			non_IUPACUnambiguousDNA = set(self.seq) - set("atcg")
			# ExtendedIUPACDNA : GATCBDSW http://biopython.org/DIST/docs/api/Bio.Alphabet.IUPAC.ExtendedIUPACDNA-class.html
			# IUPACAmbiguousDNA : GATCRYWSMKHBVDN http://biopython.org/DIST/docs/api/Bio.Alphabet.IUPAC.IUPACAmbiguousDNA-class.html
			# IUPACUnambiguousDNA : GATC http://biopython.org/DIST/docs/api/Bio.Alphabet.IUPAC.IUPACUnambiguousDNA-class.html
			
			if len(non_IUPACUnambiguousDNA) > 0:
				if len(non_IUPACUnambiguousDNA) == 1 and "n" in non_IUPACUnambiguousDNA:
					self.warning = "Found ambiguouse base: N"
				else:
					self.error = ("Found bases other than A, T, G, C or N: " + ",".join(non_IUPACUnambiguousDNA))
			
		elif filetype == "pf":
			#------------------------------------------------------------------
			# Make "sure" sequence is Protein 
			#------------------------------------------------------------------
			non_ExtendedIUPACProtein = set(self.seq) - set("ACDEFGHIKLMNPQRSTVWYBXZJUO".lower())
			# ExtendedIUPACProtein: ACDEFGHIKLMNPQRSTVWYBXZJUO http://biopython.org/DIST/docs/api/Bio.Alphabet.IUPAC.ExtendedIUPACProtein-class.html
			 
			if len(non_ExtendedIUPACProtein) != 0:
				if len(non_ExtendedIUPACProtein) == 1 and "*" in non_ExtendedIUPACProtein:
					self.warning = "Found ambiguouse symbol: *"
				else:
					self.error = ("Found amino acids other than ACDEFGHIKLMNPQRSTVWYBXZJUO = " + ",".join(non_ExtendedIUPACProtein))

	#------------------------------------------------------------------
	# Return errors 
	#------------------------------------------------------------------
	def errorHandling(self):
		return self.error, self.warning


	#------------------------------------------------------------------
	# Fetch or create oganism id, auto increment 
	#------------------------------------------------------------------
	def get_org_id(self, org_name, dbname):
		db = mdb.connect("localhost","asp","1234",dbname) # test of connection has to be performed prior
		cursor = db.cursor()
		cursor.execute("SELECT org_id from organism where name LIKE \"%" + org_name + "%\";")
		org_id = cursor.fetchone() # first element of tuple

		if not org_id:
			try:
				cursor.execute("INSERT IGNORE INTO organism(name) VALUES(\"" + org_name + "\");")
				org_id = cursor.lastrowid
				db.commit()
				db.close()
				
			except mdb.Error, e:
				sys.exit("# ERROR %d: %s" % (e.args[0],e.args[1]))
		
			return org_id
		else:
			return org_id[0]

	#------------------------------------------------------------------
	# Load data into table 
	#------------------------------------------------------------------
	def load(self, org_id, jgi, seq_name, seq, filetype, dbname, tpid):
		db = mdb.connect("localhost","asp","1234",dbname) # test of connection has to be performed prior
		cursor = db.cursor()
		values = "\",\"".join(map(str,(org_id ,jgi ,seq_name , seq,tpid)))
		try:
			
			if filetype == "pf":
				cursor.execute("INSERT IGNORE INTO proteins(org_id, prot_jgi, prot_seq_name, prot_seq, prot_proteinId) VALUES(\"" + values + "\");")
			elif filetype == "tf":
				cursor.execute("INSERT IGNORE INTO transcripts(org_id, trans_jgi, trans_seq_name, trans_seq,trans_transcriptId) VALUES(\"" + values + "\");")
			elif filetype == "cf":
				cursor.execute("INSERT IGNORE INTO cds(org_id, cds_jgi, cds_seq_name, cds_seq, cds_transcriptId) VALUES(\"" + values + "\");")
					
			elif filetype == "uma" or filetype == "ma":
				values = "\",\"".join(map(str,(org_id ,jgi ,seq_name , seq)))
				cursor.execute("INSERT IGNORE INTO assembly(org_id, assembly_jgi, assembly_seq_name, assembly_seq) VALUES(\"" + values + "\");")
			db.commit()
			
				
		except mdb.Error, e:
			sys.exit("# ERROR %d: %s" % (e.args[0],e.args[1]))
			
		db.close()
"""
#------------------------------------------------------------------
#------------------------------------------------------------------
class Gff:
#------------------------------------------------------------------
#------------------------------------------------------------------
	#------------------------------------------------------------------
	# Initialize
	#------------------------------------------------------------------
	def __init__(self, org_name, entry, filename, filetype, dbname):
		self.entry = entry
		self.desc= "This is a GFF object from JGI data"
		self.filetype = filetype
		
		fields = entry.split("\t")
		self.seqorigin = fields[0]
		self.gtype = fields[2]
		self.start = fields[3]
		self.end = fields[4]
		self.score = fields[5]
		self.strand = fields[6]
		self.phase = fields[7]
		self.attributes = fields[-1]
		self.name = self.attributes.split("\"")[1]
		self.dbname = dbname
		self.filename = filename
		self.length = abs(int(self.start) - int(self.end))
		self.org_name = org_name

		self.proteinId = ''
		self.transcriptId = ''

		if re.search(r"proteinId (\d+)", self.attributes):
			self.proteinId =  re.search(r"proteinId (\d+)", self.attributes).group(1)
		elif re.search(r"transcriptId (\d+)", self.attributes):
			self.transcriptId =  re.search(r"transcriptId (\d+)", self.attributes).group(1)
			
		self.error = ''
		
	#------------------------------------------------------------------
	# Return errors 
	#------------------------------------------------------------------
	def errorHandling(self):
		return self.error

	#------------------------------------------------------------------
	# Fetch or create oganism id, auto increment 
	#------------------------------------------------------------------
	def get_org_id(self, org_name, dbname):
		db = mdb.connect("localhost","asp","1234",dbname) # test of connection has to be performed prior
		cursor = db.cursor()
		cursor.execute("SELECT org_id from organism where name LIKE \"%" + org_name + "%\";")
		org_id = cursor.fetchone() # first element of tuple

		if not org_id:
			try:
				cursor.execute("INSERT IGNORE INTO organism(name) VALUES(\"" + org_name + "\");")
				org_id = cursor.lastrowid
				db.commit()
				db.close()
				
			except mdb.Error, e:
				sys.exit("# ERROR %d: %s" % (e.args[0],e.args[1]))

			return org_id
                
		else:
			return org_id[0]

	#------------------------------------------------------------------
	# Load data into table 
	#------------------------------------------------------------------
	def load(self, org_id, name, seqorigin, gtype, start, end, score, strand, phase, attributes, length, proteinId, transcriptId, dbname):
		db = mdb.connect("localhost","asp","1234",dbname) # test of connection has to be performed prior
		cursor = db.cursor()
		values = "\",\"".join(map(str,(org_id ,name, seqorigin, gtype, start, end, score, strand, phase, attributes.replace("\"", ""),length, proteinId,transcriptId)))
		try:
			cursor.execute("INSERT IGNORE INTO gff(org_id, gff_name, gff_seqorigin, gff_type, gff_start, gff_end, gff_score, gff_strand, gff_phase, gff_attributes, gff_length, gff_proteinId, gff_transcriptId) VALUES(\"" + values + "\");")
			db.commit()

		except mdb.Error, e:
			sys.exit("# ERROR %d: %s" % (e.args[0],e.args[1]))
			
		db.close()
"""
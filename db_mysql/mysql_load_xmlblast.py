#!/usr/bin/python

'''##################################################################
# Imports
##################################################################'''
import sys,os
sys.path.append(os.path.join(os.path.dirname(__file__), "../utils/"))
from aspmine_imports import *

'''##################################################################
# ARGUMENTS and setup
##################################################################'''
startTime = datetime.now() # record runtime
parser = CustomArgumentParser(formatter_class=SmartFormatter, 
								description="Get best hits from specific gene name or ID", 
								 usage='%(prog)s -dbname [database name] -blast [BLAST XML output file]\n'
								 "Example: python %(prog)s -dbname aspminedb -out hitlist.csv -gene An15g00560 -strain Aspni_DSM_1")
parser.add_argument("-blast", "-b", required=True, help="Name of BLAST XML file")
parser.add_argument("-dbname", "-d", required=False, default = "aspminedb", help="Database name")
#parser.add_argument("-tab", "-t", required=False, action='store_true', help="Create table new or use existing tmp")
#parser.add_argument("-strain", "-s", required=True, default=[None], action='store', help="JGI organism keys")
#parser.add_argument("-gene", "-g",  required=True, default=[None], action='store', help="List of gene ID or name")

args 	= parser.parse_args()
blast 	= args.blast

'''------------------------------------------------------------------
# Print argument values to screen
------------------------------------------------------------------'''
print "#--------------------------------------------------------------\n\
# ARGUMENTS\n\
# Description\t\t: Script return list of best bidirectional hit for a gene in a specific strain\n\
# Database\t\t: %s\n\
# BLAST XML file\t: %s\n\
#--------------------------------------------------------------" % (args.dbname, blast)

'''------------------------------------------------------------------
# Connect to specific DB
------------------------------------------------------------------'''
try:
    db = mdb.connect("localhost","asp","1234",str(args.dbname))
except mdb.Error, e:
    sys.exit("# ERROR %d: %s" % (e.args[0],e.args[1]))

try:
	cursor = db.cursor()
except mdb.Error, e: 
	sys.exit("# ERROR %d: %s" % (e.args[0],e.args[1]))


'''------------------------------------------------------------------
# Read file
------------------------------------------------------------------'''
from Bio import SearchIO

if not cursor.execute("SHOW TABLES LIKE 'map_niger' ;"):
	print "# Creating table map_niger"
	cursor.execute("CREATE TABLE map_niger ( 		\
	aln_span INT(11), \
	alphabet VARCHAR(11), \
	bitscore FLOAT, \
	bitscore_raw FLOAT, \
	evalue FLOAT, \
	gap_num INT(11), \
	hit_description VARCHAR(100), \
	hit_end INT(11), \
	hit_frame VARCHAR(11), \
	hit_id VARCHAR(100), \
	hit_range VARCHAR(11), \
	hit_span INT(11), \
	hit_start INT(11), \
	hit_strand INT(11), \
	ident_num INT(11), \
	is_fragmented VARCHAR(11), \
	pos_num INT(11), \
	query_description VARCHAR(100), \
	query_end INT(11), \
	query_frame VARCHAR(11), \
	query_id VARCHAR(100), \
	query_range VARCHAR(11), \
	query_span INT(11), \
	query_start INT(11), \
	query_strand INT(11), \
	pident FLOAT \
	)")
	curser.commit()

counter = 0
totalcounter = 0
insert_values = []

for item in SearchIO.parse(blast, 'blast-xml'): #NCBIXML.parse(result):
	for hsp in item.hsps:
		counter += 1

		pident = hsp.ident_num / hsp.aln_span * 100
		values = (	hsp.aln_span, hsp.alphabet, hsp.bitscore, hsp.bitscore_raw, hsp.evalue, hsp.gap_num, \
					hsp.hit_description, hsp.hit_end, hsp.hit_frame, hsp.hit_id, hsp.hit_range,  \
					hsp.hit_span, hsp.hit_start, hsp.hit_strand, hsp.ident_num, hsp.is_fragmented, hsp.pos_num,hsp.query_description, hsp.query_end, \
					hsp.query_frame, hsp.query_id, hsp.query_range, hsp.query_span, hsp.query_start,hsp.query_strand, pident)

		counter_threshold = 1000
		if ( counter == counter_threshold or totalcounter == len(records)):
			if totalcounter % 3000 == 0  or totalcounter == len(records) : print "# INFO: Inserting record number %s" % totalcounter

			try:	
				db = mdb.connect("localhost","asp","1234",dbname)
				cursor = db.cursor()
				cursor.executemany(query,insert_values)
				
				# add the changes to the db, close the db connection 
				db.commit() 
				db.close()	
				
				# reset counter and values to insert ------------------------------------
				counter = 0 
				insert_values = []

			except mdb.Error, e:
				sys.exit( "# ERROR utils %s %d: %s" % (filetype, e.args[0],e.args[1]))		
	sys.exit()

"""
hsp.aln, hsp.aln_all, hsp.aln_annotation, hsp.aln_annotation_all, hsp.aln_span, 
hsp.alphabet, hsp.bitscore, hsp.bitscore_raw, hsp.evalue, hsp.fragment, hsp.fragments, hsp.gap_num, 
hsp.hit, hsp.hit_all, hsp.hit_description, hsp.hit_end, hsp.hit_end_all, hsp.hit_features, hsp.hit_features_all, 
hsp.hit_frame, hsp.hit_frame_all, hsp.hit_id, hsp.hit_inter_ranges, hsp.hit_inter_spans, hsp.hit_range, hsp.hit_range_all, 
hsp.hit_span, hsp.hit_span_all, hsp.hit_start, hsp.hit_start_all, hsp.hit_strand, hsp.hit_strand_all, 
hsp.ident_num, hsp.is_fragmented, hsp.pos_num, hsp.query, hsp.query_all, hsp.query_description, hsp.query_end, hsp.query_end_all, 
hsp.query_features, hsp.query_features_all, hsp.query_frame, hsp.query_frame_all, hsp.query_id, hsp.query_inter_ranges, hsp.query_inter_spans, 
hsp.query_range, hsp.query_range_all, hsp.query_span, hsp.query_span_all, hsp.query_start, hsp.query_start_all, hsp.query_strand, hsp.query_strand_all]


hsp.aln, hsp.aln_annotation, hsp.fragment, hsp.fragments, hsp.hit, hsp.hit_features,hsp.hit_inter_ranges, hsp.hit_inter_spans, hsp.query, 
hsp.query_features, hsp.query_inter_ranges, hsp.query_inter_spans, 

hsp.aln_span, hsp.alphabet, hsp.bitscore, hsp.bitscore_raw, hsp.evalue, hsp.gap_num, 
hsp.hit_description, hsp.hit_end, hsp.hit_frame, hsp.hit_id, hsp.hit_range,  
hsp.hit_span, hsp.hit_start, hsp.hit_strand, hsp.ident_num, hsp.is_fragmented, hsp.pos_num,hsp.query_description, hsp.query_end, 
hsp.query_frame, hsp.query_id, hsp.query_range, hsp.query_span, hsp.query_start,hsp.query_strand, pident



h.aln, h.aln_annotation, h.aln_span, h.alphabet, 
h.hit, h.hit_description, h.hit_end, h.hit_features, h.hit_frame, h.hit_id, h.hit_range, h.hit_span, h.hit_start, h.hit_strand, 
h.query, h.query_description, h.query_end, h.query_features, h.query_frame, h.query_id, h.query_range, h.query_span, h.query_start, h.query_strand']

"""
"""
for alignment in item.alignments:
counter += 1
totalcounter += 1
#print dir(alignment)
#print alignment.accession, alignment.hit_def, alignment.hit_id, alignment.length, alignment.title, alignment.hsps

for hsp in alignment.hsps:
#print hsp.match
#sys.exit()
#print dir(hsp)
# align_length', 'bits', 'expect', 'frame', 'gaps', 'identities', 'match', 'num_alignments', 'positives', 'query', 'query_end', 'query_start', 'sbjct', 'sbjct_end', 'sbjct_start', 'score', 'strand'

scov = hsp.align_length / len(hsp.sbjct) *100
qcov = hsp.align_length / len(hsp.query) *100 #( max(hsp.query_end,hsp.query_start) - min(hsp.query_end,hsp.query_start) ) #/ hsp.align_length *100

values = (str(alignment.accession), str(alignment.hit_id), str(alignment.hit_def), hsp.align_length, \
	hsp.bits, hsp.expect, hsp.frame[0], hsp.gaps, hsp.identities, hsp.num_alignments, hsp.positives, \
	hsp.query_end, hsp.query_start, hsp.sbjct_end, hsp.sbjct_start, hsp.score, hsp.strand[0], qcov, scov)
valuenames = "align_acc, align_hit_id, align_hit_def, align_length, bits, expect, frame, gaps, identities, num_alignments, positives, query_end, query_start, sbjct_end, sbjct_start, score, strand, qcov, scov"
query = ("INSERT INTO map_niger SET (%s) VALUES (%s)" % ( valuenames, ("%s," * len(values)).rstrip(",")))


#if scov != 100:
print (str(alignment.accession), hsp.align_length, len(hsp.sbjct), len(hsp.query), hsp.num_alignments) 

#print hsp.sbjct
#print hsp
#sys.exit()

# Create sets of lists of values ------------------------------------
insert_values.append(values) 

"""
"""
# MySQL will only allow a certain data insert size ------------------------------------
counter_threshold = 1000
if ( counter == counter_threshold or totalcounter == len(records)):
	if totalcounter % 3000 == 0  or totalcounter == len(records) : print "# INFO: Inserting record number %s" % totalcounter

	try:	
		db = mdb.connect("localhost","asp","1234",dbname)
		cursor = db.cursor()
		cursor.executemany(query,insert_values)
		
		# add the changes to the db, close the db connection 
		db.commit() 
		db.close()	
		
		# reset counter and values to insert ------------------------------------
		counter = 0 
		insert_values = []

	except mdb.Error, e:
		sys.exit( "# ERROR utils %s %d: %s" % (filetype, e.args[0],e.args[1]))		
#print ("%s:%s:%s:%s:%s:%s:%s") \
#% 
#print type(hsp)
  	#print hsp
"""  	
"""
 if hsp.expect <0.01:
         print('****Alignment****')
         print('sequence:', alignment.title) 
         print('length:', alignment.length)
         print('score:', alignment.score)
         print('gaps:', alignment.gaps)
         print('e value:', hsp.expect)
         print(hsp.query[0:90] + '...')
         print(hsp.match[0:90] + '...')
         print(hsp.sbjct[0:90] + '...')
"""    
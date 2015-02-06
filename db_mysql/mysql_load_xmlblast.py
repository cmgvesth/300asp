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
	cursor.execute("CREATE TABLE IF NOT EXISTS `map_niger` (\
  `filename` varchar(50) DEFAULT NULL, \
  `hit_id` varchar(100) DEFAULT NULL,\
  `hit_org` varchar(100) DEFAULT NULL,\
  `query_id` varchar(100) DEFAULT NULL,\
  `query_org` varchar(100) DEFAULT NULL,\
  `query_range` varchar(11) DEFAULT NULL,\
  `hit_range` varchar(11) DEFAULT NULL,\
  `aln_span` int(11) DEFAULT NULL,\
  `alphabet` varchar(11) DEFAULT NULL,\
  `bitscore` decimal(10,2) DEFAULT NULL,\
  `bitscore_raw` decimal(10,2) DEFAULT NULL,\
  `evalue` varchar(11) DEFAULT NULL,\
  `gap_num` int(11) DEFAULT NULL,\
  `hit_end` int(11) DEFAULT NULL,\
  `hit_frame` varchar(11) DEFAULT NULL,\
  `hit_span` int(11) DEFAULT NULL,\
  `hit_start` int(11) DEFAULT NULL,\
  `hit_strand` int(11) DEFAULT NULL,\
  `ident_num` int(11) DEFAULT NULL,\
  `is_fragmented` varchar(11) DEFAULT NULL,\
  `pos_num` int(11) DEFAULT NULL,\
  `query_end` int(11) DEFAULT NULL,\
  `query_frame` varchar(11) DEFAULT NULL,\
  `query_span` int(11) DEFAULT NULL,\
  `query_start` int(11) DEFAULT NULL,\
  `query_strand` int(11) DEFAULT NULL,\
  `pident` decimal(10,2) DEFAULT NULL,\
  UNIQUE KEY `filename` (`filename`,`hit_id`,`query_id`,`bitscore_raw`,`pident`),\
  KEY `hit_id` (`hit_id`),\
  KEY `query_id` (`query_id`)\
) ENGINE=InnoDB DEFAULT CHARSET=latin1;")
	db.commit()

counter = 0
totalcounter = 0
insert_values = []

records = SearchIO.parse(blast, 'blast-xml')

for item in records: #NCBIXML.parse(result):
	for hsp in item.hsps:
		#print hsp.evalue 
		#sys.exit()
		filename = (blast.split("/")[-1]).replace(".txt", "")
		qorg = filename.split("_")[0]
		sorg = filename.split("_")[-1]

		pident = (float(hsp.ident_num) / float(hsp.aln_span)) * 100
		values = (	filename, hsp.aln_span, hsp.alphabet, hsp.bitscore, hsp.bitscore_raw, hsp.evalue, hsp.gap_num, \
					hsp.hit_end, hsp.hit_frame, hsp.hit_id, str(hsp.hit_range),  \
					hsp.hit_span, hsp.hit_start, hsp.hit_strand, hsp.ident_num, hsp.is_fragmented, hsp.pos_num, hsp.query_end, \
					hsp.query_frame, hsp.query_id, str(hsp.query_range), hsp.query_span, hsp.query_start,hsp.query_strand, pident, qorg, sorg)
		valuenames = "filename, aln_span, alphabet, bitscore, bitscore_raw, evalue, gap_num, \
					hit_end, hit_frame, hit_id, hit_range,  \
					hit_span, hit_start, hit_strand, ident_num, is_fragmented, pos_num, query_end, \
					query_frame, query_id, query_range, query_span, query_start,query_strand, pident, query_org, hit_org"

		query = "REPLACE INTO map_niger (%s) VALUES (%s)" % (valuenames, ("%s," * len(values)).rstrip(","))
		insert_values.append(values)
		counter_threshold = 100

		
		if ( counter == counter_threshold):
			if totalcounter % 3000 == 0 : print "# INFO: Inserting record number %s" % totalcounter

			try:	
				cursor.executemany(query, insert_values)
				
				# add the changes to the db, close the db connection 
				db.commit() 
				
				# reset counter and values to insert ------------------------------------
				counter = 0 
				insert_values = []

			except mdb.Error, e:
				sys.exit( "# ERROR xmlblast %d: %s" % (e.args[0],e.args[1]))	
		
		counter += 1
		totalcounter += 1

"""
select * from (
select query_id, query_org, 
case when hit_org = "Aniger1" then hit_id end as Aniger1, 
case when hit_org = "Aniger3" then hit_id end as Aniger3, 
case when hit_org = "Aniger7" then hit_id end as Aniger7 
from map_niger) as mapA
group by query_id;

"""

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
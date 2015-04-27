#!/usr/bin/python

'''##################################################################
# Imports
##################################################################'''
import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), "../utils/"))
from aspmine_imports import *

"""
CODEBOOK INFO:

blastp 
-query ./FASTA_prot/A_acidus_CBS_106_47_orf_trans_all.fasta 
-db ./BLAST_db/Aindologenus.aa 
-out ./BLAST_out3/AacidusVsAindologenusTable.txt 
-evalue 1e-100 
-outfmt "6 qseqid sseqid pident qlen qstart qend slen sstart send evalue bitscore"



for x in Az*/*txt
do
python /home/tcve/github/db_mysql/mysql_blast_load.py -filetype blast -action load -source $x
done


blast | CREATE TABLE `blast` (
  `blast_id` int(11) unsigned zerofill NOT NULL AUTO_INCREMENT,
  `blast_filename` varchar(100) NOT NULL,
  `blast_qseq_id` varchar(100) NOT NULL DEFAULT '',
  `blast_sseq_id` varchar(100) NOT NULL DEFAULT '',
  `blast_pident` varchar(50) DEFAULT NULL,
  `blast_qlen` int(11) DEFAULT NULL,
  `blast_qstart` int(11) NOT NULL DEFAULT '0',
  `blast_qend` int(11) DEFAULT NULL,
  `blast_slen` int(11) DEFAULT NULL,
  `blast_sstart` int(11) DEFAULT NULL,
  `blast_send` int(11) DEFAULT NULL,
  `blast_bitscore` int(11) NOT NULL DEFAULT '0',
  `blast_evalue` varchar(10) DEFAULT NULL,
  PRIMARY KEY (`blast_bitscore`,`blast_qseq_id`,`blast_sseq_id`,`blast_qstart`),
  UNIQUE KEY `id` (`blast_id`),
  KEY `blast_qseq_id` (`blast_qseq_id`),
  KEY `blast_sseq_id` (`blast_sseq_id`)
) ENGINE=InnoDB AUTO_INCREMENT=23820677 DEFAULT CHARSET=latin1 |


for x in A*/*txt; do python /home/tcve/Dropbox/aspmine/mysql_blast_load.py -filetype blast -action load -source $x; done

for x in A*cbs*/*euc*txt A*cbs*/*kaw*txt A*cbs*/*cbs*txt A*euc*/*euc*txt A*euc*/*kaw*txt A*euc*/*cbs*txt A*kaw*/*euc*txt A*kaw*/*kaw*txt A*kaw*/*cbs*txt; do python /home/tcve/Dropbox/aspmine/mysql_blast_load.py -filetype blast -action load -source $x; done

select * from proteins where prot_jgi3 LIKE "An08g11070%" or prot_jgi3 LIKE "An15g00320%" or prot_jgi3 LIKE "An06g02420%" 


"""
			
#------------------------------------------------------------------
# Format argument helper
#------------------------------------------------------------------
class SmartFormatter( argparse.HelpFormatter ):
	width = 100
	def _split_lines( self, text, width ):
	#------------------------------------------------------------------
	# this is the RawTextHelpFormatter._split_lines
	#------------------------------------------------------------------
		if text.startswith( 'R|' ):
			return text[2:].splitlines()  
		return argparse.HelpFormatter._split_lines( self, text, width )

#------------------------------------------------------------------
# Get command line arguments
#------------------------------------------------------------------
parser = argparse.ArgumentParser( formatter_class=SmartFormatter, 
								 usage='%(prog)s -filetype filetype -action action -source filename\n'
								 "Example: python %(prog)s -dbname aspminedb -filetype blast -action test -source AacidusVsAindologenusTable.txt\n" )

#------------------------------------------------------------------
# Choose if data should be loaded or only tested
#------------------------------------------------------------------
parser.add_argument( '-action' , required=False, default = "test", choices=['load', 'test'],
					help='R|Test or test and load data\n'
					"-test\t: Read files and verify formats\n"
					"-load\t: Read files and verify formats followed by storing data in MySQL tables\n" )

#------------------------------------------------------------------
# At least one file filetype must be provided, or the option "all" which expects a JGI directory format
#------------------------------------------------------------------
parser.add_argument( '-filetype' ,	required=False,	help="R|BLAST filetype, files in .txt format.", default = "blast", choices=['blast'])
parser.add_argument( "-source",	"-s",	required=True,	help="R|File or directoryname, example: -source AacidusVsAindologenusTable.txt" )
parser.add_argument( "-dbname",		required=False,	help="R|Database name", default = "aspminedb" )

args = parser.parse_args()

#------------------------------------------------------------------
# Print argument values to screen
#------------------------------------------------------------------
print "#--------------------------------------------------------------\
# ARGUMENTS :\n\
# source\t:%s\n\
# filetype\t:%s\n\
# action\t:%s\n\
# Database\t:%s\
#--------------------------------------------------------------" % ( args.source, args.filetype, args.action, args.dbname )

#------------------------------------------------------------------
# Make connection to database
#------------------------------------------------------------------
try:
	db = mdb.connect( "localhost","asp","1234",str( args.dbname ) )
	cursor = db.cursor()
except mdb.Error, e:
	sys.exit( "# ERROR %d: %s" % ( e.args[0],e.args[1] ) )

#------------------------------------------------------------------
# Read BLAST file
#------------------------------------------------------------------
if os.path.isfile( args.source ):
	print "# INFO: file exists: %s" % args.source 
	os.path.isfile( args.source )
else:	
	sys.exit( "# ERROR: file does not exists" )

#------------------------------------------------------------------
# Get lines from file
#------------------------------------------------------------------
file_obj = open( args.source, "r" )
file_lines = file_obj.readlines()
file_obj.close()        

#------------------------------------------------------------------
# Open DB connection
#------------------------------------------------------------------
db = mdb.connect( "localhost","asp","1234",str( args.dbname ) )

#------------------------------------------------------------------
# Get sequence names using BLAST seqids and organism ids
#------------------------------------------------------------------
# Initiate variables
values_to_insert = []
counter = 0
totalcounter = 0

for line in file_lines:
	line = line.rstrip()
	(qjgi1, qjgi2, qjgi3) = ("","","")
	(sjgi1, sjgi2, sjgi3) = ("","","")

	# Always check that the file has the format which is expected - runs for each line
	if not len( line.split( "\t" ) ) == 11:	
		sys.exit( "# ERROR: lines do not have the right number of elements ( 11 ): \
			qseqid, sseqid, pident, qlen, qstart, qend, slen, sstart, send, evalue, bitscore" )

	( qorg_key, qprotein_key, sorg_key, sprotein_key, search_s, search_q ) = ( "","","","","","" )	# Define variables
	
	counter += 1		# Increment counters
	totalcounter += 1	# Increment counters
	values = ()			# Empty list of values

	( qseqid, sseqid, pident, qlen, qstart, qend, slen, sstart, send, evalue, bitscore ) = line.split( "\t" )[0:11]
	 
	if len(sseqid.split("|")) >= 4 : 	(sjgi1, sjgi2, sjgi3) 	= sseqid.split("|")[1:4]
	elif len(sseqid.split("_")) >= 2:	(sjgi1, sjgi2) 			= sseqid.split("_")[0:2]
	elif re.search("Afu", sseqid) :		(sjgi1, sjgi2) 			= ("sApfu1", sseqid)
	else: 								(sjgi1, sjgi2, sjgi3) 	= ( sseqid, "","")

	if len(qseqid.split("|")) >= 4 :	(qjgi1, qjgi2, qjgi3) 	= qseqid.split("|")[1:4]
	elif len(qseqid.split("_")) >= 2:	(qjgi1, qjgi2) 			= qseqid.split("_")[0:2]
	elif re.search("Afu", qseqid) :		(qjgi1, qjgi2) 			= ("Aspfu1", qseqid)
	else: 								(qjgi1, qjgi2, qjgi3) 	= ( qseqid, "","")

	values = (args.source.split("/")[-1] , qseqid, qjgi1, qjgi2, qjgi3 , sseqid, sjgi1, sjgi2, sjgi3, pident, qlen, qstart, qend, slen, sstart, send, evalue, bitscore )
	values = [re.sub( r"^0+", "", i ) for i in values]
	#print values
	values_to_insert.append( values ) # Create sets of lists of values
	# Only load data if action is specified
	if args.action == "load" and ( counter == 500 or totalcounter == len( file_lines ) ):
		if totalcounter % 1000 == 0 : print "# INFO: Inserting record number %s" % totalcounter
		elif totalcounter == len(file_lines): print "# INFO: Inserting record number %s" % totalcounter

	
		try:	
			cursor.executemany( "REPLACE INTO blast ( blast_filename, blast_qseq_id, blast_qseq_jg1,blast_qseq_jg2,blast_qseq_jg3, blast_sseq_id,  blast_sseq_jg1,blast_sseq_jg2,blast_sseq_jg3, blast_pident, blast_qlen, blast_qstart, blast_qend, blast_slen, blast_sstart, blast_send, blast_evalue, blast_bitscore ) values( %s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s ,%s,%s,%s,%s,%s,%s);",values_to_insert )
			db.commit()	# Add changes to database
			
			counter = 0 			# restart counter
			values_to_insert = []	# Empty list of values

		except mdb.Error, e:
			sys.exit( "# ERROR blast load %s %d: %s" % ( args.source, e.args[0],e.args[1] ) )
	
if args.action == "test":
	print "# INFO: testing file, example values:\n#%s" % values
	
print "# FINISHED BLAST file: %s with total number of records %s" % ( args.source , totalcounter )

db.close()	# Close DB connection
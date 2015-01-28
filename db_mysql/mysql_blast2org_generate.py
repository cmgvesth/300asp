#!/usr/bin/python

#------------------------------------------------------------------
# Imports
#------------------------------------------------------------------
'''##################################################################
# Imports
##################################################################'''
import sys 
sys.path.append('../utils/')
from aspmine_imports import *

"""
CODEBOOK INFO:

"""
			
#------------------------------------------------------------------
# Format argument helper
#------------------------------------------------------------------
class SmartFormatter( argparse.HelpFormatter ):
	width = 100
	def _split_lines( self, text, width ):
	# this is the RawTextHelpFormatter._split_lines
		if text.startswith( 'R|' ):
			return text[2:].splitlines()  
		return argparse.HelpFormatter._split_lines( self, text, width )

#------------------------------------------------------------------
# Get command line arguments
#------------------------------------------------------------------
parser = argparse.ArgumentParser( formatter_class=SmartFormatter, 
								 usage='%(prog)s -action connect\n'
								 "Example: python %(prog)s -dbname aspminedb -action connect\n" )

#------------------------------------------------------------------
# Choose if data should be loaded or only tested
#------------------------------------------------------------------
parser.add_argument( '-action' , required=False, default = "test", choices=['connect', 'test'],
					help='R|Test or test and load data\n'
					"-test\t: Read and verify connections between sequence names in blast table and organisms in organism table\n"
					"-connect\t: Read sequence names in MySQL blast table and connect to organisms, Stores connections in blast2org table\n" )

parser.add_argument( "-dbname", required=False, help="R|Database name", default = "aspminedb" )

args = parser.parse_args()

#------------------------------------------------------------------
# Print argument values to screen
#------------------------------------------------------------------
print "#--------------------------------------------------------------"
print "# ARGUMENTS :\n# action\t:%s\n# Database\t:%s" % ( args.action,args.dbname )
print "#--------------------------------------------------------------"

#------------------------------------------------------------------
# Connect to specific DB
#------------------------------------------------------------------
try:
    db = mdb.connect( "localhost","asp","1234",str( args.dbname ) )
    db.close() 
except mdb.Error, e:
    sys.exit( "# ERROR %d: %s" % ( e.args[0],e.args[1] ) )

#------------------------------------------------------------------
# Define functions
#------------------------------------------------------------------
#######################################################################
def orgid_from_seqid ( seqid ):
#######################################################################
	( org_key, protein_key,result, query ) = ( "","","","" )
	
	#print seqid.split( "|" )
	if len( seqid.split( "|" ) ) > 1:
		org_key = seqid.split( "|" )[1]
		protein_key = seqid.split( "|" )[2]

	elif len( seqid.split( "_" ) ) > 1:
		#if re.search("^Asp", seqid.split( "_" )[0]) or re.search("^Aa", seqid.split( "_" )[0]):
		org_key = seqid.split( "_" )[0]
		protein_key = seqid.split( "_" )[1]

		
	if org_key and protein_key:
		cursor.execute( "SELECT org_id, prot_protein_id FROM organism join proteins using ( org_id ) where name = '%s' and prot_protein_id = %s" % ( org_key, protein_key ) )
		( result,query ) = ( cursor.fetchall(), cursor._last_executed )
		
		if len( result ) > 1 or len( result ) < 1:
			protein_key = re.sub( r"^0+", "0?", protein_key )
			org_key = re.match( r"^A([a-zA-Z]+)\d+", org_key ).groups(0)[0]
			
			cursor.execute( "SELECT org_id, prot_protein_id FROM organism join proteins using ( org_id ) where real_name REGEXP '.*%s.*' and prot_protein_id REGEXP '^%s$' " % ( org_key, protein_key ) )
			( result,query ) = ( cursor.fetchall(), cursor._last_executed )
			#print query
		
	"""
	else:
		search = seqid
		search = re.sub( "_0+", "\\\\\\|*_*0*", seqid )
		cursor.execute( "SELECT org_id, prot_protein_id from proteins where prot_seq_name REGEXP '.*%s\\\\|.*';" % search )
		( result,query ) = ( cursor.fetchall(), cursor._last_executed )

		if len( result ) > 1 or len( result ) < 1:
			cursor.execute( "SELECT org_id, prot_protein_id from proteins where prot_seq_name REGEXP '.*%s$';" % search )
			( result,query ) = ( cursor.fetchall(), cursor._last_executed )	

	elif re.search( r"(Asp.+)_", seqid ) and ( len( result ) > 1 or len( result ) < 1 ):
		search = re.search( r"(Asp.+)_", seqid ).group( 1 )
		cursor.execute( "SELECT org_id, prot_protein_id FROM organism join proteins using ( org_id ) where name REGEXP '.*%s.*'" % search )
		( result,query ) = ( cursor.fetchall(), cursor._last_executed )
	"""

	return result, query	
"""
#######################################################################
def get_seqname ( org_key, seqid ):
#######################################################################
	( protein_key,result ) = ( "","" )

	if len( seqid.split( "|" ) ) > 1: 
		protein_key = seqid.split( "|" )[2]
		cursor.execute( "SELECT prot_seq_name FROM organism join proteins using ( org_id ) where org_id = '%s' and prot_protein_id = %s" % ( org_key, protein_key ) )
		( result,query ) = ( cursor.fetchall(), cursor._last_executed )

	else:
		search = re.sub( "_0+", "\\\\\\|*_*0*", seqid )
		cursor.execute( "SELECT prot_seq_name from proteins where prot_seq_name REGEXP '.*%s\\\\|.*' and org_id = '%s';" % ( search, org_key ) )
		( result,query ) = ( cursor.fetchall(), cursor._last_executed )

		if len( result ) > 1 or len( result ) < 1:
			cursor.execute( "SELECT prot_seq_name from proteins where prot_seq_name REGEXP '.*%s$' and org_id = '%s';" % ( search, org_key ) )
			( result,query ) = ( cursor.fetchall(), cursor._last_executed )

	if len( result ) > 1 or len( result ) < 1:
		( result, query ) = ( ("none",), cursor._last_executed )		
		
	return result, query			
"""
#------------------------------------------------------------------
# Open DB connection
#------------------------------------------------------------------
db = mdb.connect( "localhost","asp","1234",str( args.dbname ) )
cursor = db.cursor()

#------------------------------------------------------------------
# Get organism ids from sequence ids
#------------------------------------------------------------------
# Get list of sequence ids
# Taking all values at once is to much so you run through it in batches

cursor.execute( "select max(blast_id) from blast")
max_blast_id = cursor.fetchall()
db.commit

cursor.execute( "select min(blast_id) from blast")
min_blast_id = cursor.fetchall()
db.commit

( seq_count, lower_count) = ( 0, min_blast_id[0][0])	# Define variables

print "# INFO: ID numbers range from %s to %s" % (min_blast_id[0][0], max_blast_id[0][0])

#lower_count = 100000
while lower_count <= max_blast_id[0][0]:

	try:	
		cursor.execute( "select blast_qseq_id, blast_sseq_id from blast where blast_id between %s and %s" % (lower_count, lower_count+1000) )
		#print cursor._last_executed
		results = cursor.fetchall()
		db.commit()
		
	except mdb.Error, e:
		sys.exit( "# ERROR blast2org %d: %s" % ( e.args[0],e.args[1] ) )

	# Veryfy that ids were found
	if len(results) < 1: sys.exit("# ERROR: no sequence ids found in blast table")
	print "# INFO: processing entries with ID numbers between %s and %s, sequence ID list %s" % (lower_count, lower_count+1000, len(results))
	lower_count += 1000
		
	# Loop through each pair and obtain 
	# sequence name from protein table - prot_seq_name
	# organism id from organism table - org_id
	unfound_genes = open('unfound_genes.txt', 'w') 

	values_to_insert = []
	counter = 0
	totalcounter = 0

	for pair in results:
		counter += 1
		totalcounter += 1
		values = ()

		qorg = pair[0]
		sorg = pair[1]
		
		(qresults, qquery) = orgid_from_seqid( qorg )
		(sresults, squery) = orgid_from_seqid( sorg )

		( qorg_id , qprotein_id ) = ( qresults[0][0] ,  qresults[0][1] )
		( sorg_id , sprotein_id ) = ( sresults[0][0] ,  sresults[0][1] )

		values = (qprotein_id,  sprotein_id, qorg_id, sorg_id, qorg, sorg)

		#print "# sequence ids: %s -- %s -- %s -- %s -- %s -- %s" % (qorg,  sorg, qorg_id, qprotein_id, sorg_id,  sprotein_id)

		if not sorg_id or not qorg_id:
			seq_count += 1
			if not sorg_id:	unfound_genes.write("subject:\t" + pair[0] + "\n")
			if not qorg_id:	unfound_genes.write("query:\t" + pair[1] + "\n")
		else:
			values_to_insert.append( values )
			
		counter_threshold = 1000
	
		if args.action == "connect" and ( counter == counter_threshold or totalcounter == len(results)):
			if totalcounter % 3000 == 0  or totalcounter == len(results) : print "# INFO: Inserting record number %s" % totalcounter

			try:	
				cursor = db.cursor()
				cursor.executemany("REPLACE INTO blast2org (qseq_id, sseq_id, qorg_id, sorg_id, qprot_seq_name, sprot_seq_name) VALUES (%s,%s,%s,%s,%s,%s)",values_to_insert)
				db.commit()
				
				counter = 0 
				values_to_insert = []

			except mdb.Error, e:
				sys.exit( "# ERROR utils %d: %s" % (e.args[0],e.args[1]))

		#print "# FINISHED %s file: %s with total number of records %s" % (name,filepath.split("/")[-1], totalcounter)
"""
# Error if organism ids could not be found
if not qorg_id or not sorg_id:
	sys.exit("# ERROR: organism ID could not be obtained for query or subject")
print "# FINISHED BLAST file: %s with total number of records %s" % ( args.source , totalcounter )
"""
unfound_genes.close()
db.close()	# Close DB connection

print "#--------------------------------------------------------------"
print "# INFO: Final error handling"
print "#--------------------------------------------------------------"

print "# INFO total blast entries: ", totalcounter
print "# INFO total genes not connected to organism ids: ", seq_count

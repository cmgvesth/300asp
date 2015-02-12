#!/usr/bin/python

#------------------------------------------------------------------
# Imports
#------------------------------------------------------------------
'''##################################################################
# Imports
##################################################################'''
import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), "../utils/"))
from aspmine_imports import *

"""
CODEBOOK INFO:

for x in $(locate *Asphet1*_sm.txt)
do
python  python /home/tcve/github/db_mysql/mysql_antismash_load.py -filetype anti -action load -source $x
done

for x in *sm.txt
do
python /home/tcve/Dropbox/aspmine/mysql_antismash_load.py -filetype anti -action load -source $x
done

head -n 20 Neofi1_sm.txt 
60 SM backbone genes predicted
PKS 16
HYBRID 1
PKS-Like 2
NRPS-Like 10
DMAT 10
TC 1
NRPS 20
#########################
Cluster size 6  for backbone gene  3784 N_fischeri_1099437636252   72056   79171   NRPS
   16  3780         #  KR adh_short
   17  3781         #  MFS_1
   18  3782         #  Fungal_trans_2
   19  3783         #  
   20  3784  NRPS  #  AMP-binding PP-binding Condensation
   21  3785         #  MFS_1 Sugar_tr
Cluster size 1  for backbone gene  3827 N_fischeri_1099437636252  190763  193849   NRPS-Like
   63  3827  NRPS-Like  #  AMP-binding Condensation
Cluster size 5  for backbone gene  3835 N_fischeri_1099437636252  216976  218337   DMAT
   71  3835  DMAT  #  Trp_DMAT


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
								 "Example: python %(prog)s -dbname aspminedb -filetype anti -action test -source Aspzo1_sm.txt\n" )

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
parser.add_argument( '-filetype', "-f",	required=False,	help="R|Antismash JGI filetype, files in .txt format.", default = "anti", choices=['anti'],  )
parser.add_argument( "-source",	"-s",	required=True,	help="R|File or directoryname, example: -source Aspzo1_sm.txt" )
parser.add_argument( "-dbname",	"-d", 	required=False,	help="R|Database name", default = "aspminedb" )

args = parser.parse_args()

#------------------------------------------------------------------
# Print argument values to screen
#------------------------------------------------------------------
print "#--------------------------------------------------------------"
print "# ARGUMENTS :\n# source\t:%s\n# filetype\t:%s\n# action\t:%s\n# Database\t:%s" % ( args.source, args.filetype, args.action, args.dbname )
print "#--------------------------------------------------------------"

#------------------------------------------------------------------
# Chck connection to database
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
# Get sequence names using BLAST seqids and organism ids
#------------------------------------------------------------------
# Initiate variables
values_to_insert = []
counter = 0
totalcounter = 0
cluster_flag = 0
org_id = ""
filename = args.source.split("/")[-1]

for line in file_lines:
	totalcounter += 1
	
	if re.search("Cluster", line):
		cluster_info = line.rstrip()
		cluster_info = re.sub("\s\s+", " ", cluster_info)
		cluster_flag = 1
		continue

		
	if cluster_flag == 1:
		counter += 1

		line = re.sub("\s\s+", "\t", line)
		line = re.sub("^\s", "", line)
		fields = line.rstrip().split("\t")
		
		cluster_fields = cluster_info.split(" ")
		(cluster, cluster_size, backbone, origin, start, end, cluster_type) = (cluster_info, cluster_fields[2], cluster_fields[6], cluster_fields[7], cluster_fields[8], cluster_fields[9], cluster_fields[10] )
		if len(fields) == 3:
			(sm_id, protein_id, sm, sm_desc) = (fields[0], fields[1], "none", "none")
		elif len(fields) == 4:
			(sm_id, protein_id, sm, sm_desc) = (fields[0], fields[1], "none", fields[3])
		elif len(fields) == 5:
			(sm_id, protein_id, sm, sm_desc) = (fields[0], fields[1], fields[2], fields[4])
			
		if re.search("Cluster", line):
			cluster_info
			cluster_flag = 0

		if not org_id:
			org_key = filename.split("_sm.txt")[0]
			cursor.execute( "SELECT org_id FROM organism where name LIKE '%s' " % org_key )
			org_id = cursor.fetchall()
			
			if not len(org_id) == 1:
				
				sys.exit("# ERROR: data could note be associated with an organism")
				
		values = ( org_id[0][0], args.source, cluster_info, cluster_size, backbone, origin, start, end, cluster_type, sm_id , protein_id, sm, sm_desc)
		values_to_insert.append( values ) # Create sets of lists of values

		# Only load data if action is specified
		if args.action == "load" and ( counter == 100 or totalcounter == len( file_lines ) ):
			if totalcounter % 1000 == 0 : print "# INFO: Inserting record number %s" % totalcounter
			elif totalcounter == len(file_lines): print "# INFO: Inserting record number %s" % totalcounter

			try:	
				cursor.executemany( "REPLACE INTO antismash ( org_id, filename, cluster, clust_size, clust_backbone, clust_origin, clust_start, clust_end, clust_type, sm_id, sm_protein_id, sm_short, sm_desc ) values( %s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s);",values_to_insert )
				db.commit()	# Add changes to database
				
				counter = 0 			# restart counter
				values_to_insert = []	# Empty list of values

			except mdb.Error, e:
				sys.exit( "# ERROR antismash load %s %d: %s" % ( args.source, e.args[0],e.args[1] ) )
		

print "# FINISHED BLAST file: %s with total number of records %s" % ( args.source , totalcounter )

db.close()	# Close DB connection
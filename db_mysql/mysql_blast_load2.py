#!/usr/bin/python


"""
TOtal number of lines:
wc -l `find . -type f`
32.567.719


List all Files: 
ls -1 */* | grep Vs | cut -d "/" -f 2 | sort > allfiles.tab

mysql -u asp -p1234 -B -e 'select distinct(filename) from blast;' aspminedb | grep Vs | sort > allfiles_mysql.tab

Difference in filenames:
tcve@klara:~/ExtraStorage/BLASTOutPuts$ comm -13 allfiles_mysql.tab allfiles.tab | wc = 67
tcve@klara:~/ExtraStorage/BLASTOutPuts$ wc allfiles_mysql.tab 	= 2101  2101 73310 allfiles_mysql.tab
tcve@klara:~/ExtraStorage/BLASTOutPuts$ wc allfiles.tab 		= 2914   2914 101335 allfiles.tab


Aaculeatinus/AaculeatinusVsAacidusTable.txt:jgi|Aspacu1|387311|fgenesh1_kg.1_#_254_#_Locus5009v3rpkm2.98	jgi|Aspfo1|207813|gm1.6437_g	43.61	874	22	874	864	22	864	0.0	  558
Aaculeatinus/AaculeatinusVsAfoetidusTable.txt:jgi|Aspacu1|387311|fgenesh1_kg.1_#_254_#_Locus5009v3rpkm2.98	jgi|Aspfo1|207813|gm1.6437_g	43.61	874	22	874	864	22	864	0.0	  558

Hej. Jeg har et lille problem med nogle af vores BLAST.
Det er noget jeg ahr opdaget pga de test mine scripts laver. 
Naar jeg loader BLAST filerne ind far jeg et problem med duplikerede entries:
Aaculeatinus/AaculeatinusVsAacidusTable.txt:jgi|Aspacu1|387311|fgenesh1_kg.1_#_254_#_Locus5009v3rpkm2.98	jgi|Aspfo1|207813|gm1.6437_g	43.61	874	22	874	864	22	864	0.0	  558
Aaculeatinus/AaculeatinusVsAfoetidusTable.txt:jgi|Aspacu1|387311|fgenesh1_kg.1_#_254_#_Locus5009v3rpkm2.98	jgi|Aspfo1|207813|gm1.6437_g	43.61	874	22	874	864	22	864	0.0	  558

Som det ser ud her er gener of JGI kode for Aacidus og Afoetidus den samme (jgi|Aspfo1|207813|gm1.6437_g)
Naar man slaar Afoetidus op paa JGI faar man kun entriy for Aacidus.

Jeg kan ikke gemme begge disse. De har samme orgkey og sekvens key. 


for d in $(ls -d */); do  python /home/tcve/github/db_mysql/mysql_blast_load2.py -s $d -action load -e $d ; done



select count(distinct(filename)) from blast;
+---------------------------+
| count(distinct(filename)) |
+---------------------------+
|                      2899 |
+---------------------------+

"""
#------------------------------------------------------------------
# Imports
#------------------------------------------------------------------
'''##################################################################
# Imports
##################################################################'''
import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), "../utils/"))
from aspmine_imports import *
	
startTime = datetime.now() # record runtime

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

parser.add_argument( "-source",	"-s",required=True,	help="R|File or directoryname, example: -source AacidusVsAindologenusTable.txt" )
parser.add_argument( "-dbname",		required=False,	help="R|Database name", default = "aspminedb" )
parser.add_argument("-replace", required=False, action='store_true', help="Replace data if file already loaded")
parser.add_argument("-err", "-e", required=False, default = "default", help="Replace data if file already loaded")

args = parser.parse_args()
replace = args.replace

if args.err.endswith("/"):
	errorname = args.err.split( "/" )[0]
else:
	errorname =	args.err
errorlogfile = os.path.basename(sys.argv[0])+ "_" + errorname +".err"
sys.stdout = Logger(errorlogfile)
sys.stderr = Logger(errorlogfile)

#------------------------------------------------------------------
# Print argument values to screen
#------------------------------------------------------------------
print "#--------------------------------------------------------------\n\
# ARGUMENTS :\n\
# Description: loading BLAST tab files, ignores selfhits based on ID\n\
# source\t:%s\n\
# action\t:%s\n\
# Database\t:%s\n\
# Errorlog\t:%s\n\
#--------------------------------------------------------------" % ( args.source, args.action, args.dbname, errorlogfile )

path = args.source

if os.path.isfile(path):
	allpath = [path,]

elif os.path.isdir(path):
	allpath = ["/".join([path, f]) for (root,dirs,files) in os.walk(path) for f in files if f.endswith(".txt")]

else:
	sys.exit()
print "# INFO: Source is a list of %s *.txt files" % len(allpath)
			
#------------------------------------------------------------------
# Make connection to database
#------------------------------------------------------------------
try:
	db = mdb.connect( "localhost","asp","1234",str( args.dbname ) )
	cursor = db.cursor()
except mdb.Error, e:
	sys.exit( "# ERROR %d: %s" % ( e.args[0],e.args[1] ) )

#------------------------------------------------------------------
# Create table if it does not exist
#------------------------------------------------------------------
if not cursor.execute("SHOW TABLES LIKE 'blast';"):
	print "# Creating table blast"
	cursor.execute("CREATE TABLE `blast` (\
	  `q_org` varchar(100) NOT NULL,\
	  `q_seqkey` varchar(100) NOT NULL,\
	  `q_tail` varchar(100) NOT NULL,\
	  `h_org` varchar(100) NOT NULL,\
	  `h_seqkey` varchar(100) NOT NULL,\
	  `h_tail` varchar(100) NOT NULL,\
	  `pident` decimal(10,2) DEFAULT NULL,\
	  `q_len` int(11) NOT NULL,\
	  `q_start` int(11) NOT NULL,\
	  `q_end` int(11) NOT NULL,\
	  `h_len` int(11) NOT NULL,\
	  `h_start` int(11) NOT NULL,\
	  `h_end` int(11) NOT NULL,\
	  `bitscore` int(11) NOT NULL,\
	  `evalue` varchar(10) NOT NULL,\
	  `h_seqid` varchar(100) NOT NULL,\
	  `q_seqid` varchar(100) NOT NULL,\
	  `filename` varchar(100) NOT NULL,\
	  `loadstamp` varchar(100) NOT NULL,\
	  PRIMARY KEY (`bitscore`,q_org,q_seqkey,h_org, h_seqkey,`q_start`, h_start),\
	  KEY `q_org` (`q_org`),\
	  KEY `q_seqkey` (`q_seqkey`),\
	  KEY `h_org` (`h_org`),\
	  KEY `h_seqkey` (`h_seqkey`),\
	  KEY `h_seqid` (`h_seqid`),\
	  KEY `q_seqid` (`q_seqid`),\
	  KEY `filename` (`filename`)\
	) ENGINE=MyISAM DEFAULT CHARSET=latin1")
	db.commit()

#------------------------------------------------------------------
# Get lines from file

def get_lines(filepath):			
	file_obj = open( filepath, "r" )
	file_lines = file_obj.readlines()
	file_obj.close()        

	return file_lines

# END get_lines
#------------------------------------------------------------------

#------------------------------------------------------------------
# Remove rows with same filename

def clear_database_rows(filename):
	#------------------------------------------------------------------
	# Locate rows with same filename
	#------------------------------------------------------------------
	status = 0
	try:
		cursor.execute("SELECT * from blast where filename=\'" + filename +"\'")
		result = cursor.fetchone()
	except mdb.Error, e:
		sys.exit( "# ERROR blast select filename %s %d: %s" % ( filename, e.args[0],e.args[1] ) )
		
	if result:
		status = 1
		print "# WARNING: filename already found in table, run with -replace to overwrite!"

	# Delete any rows if found
	if replace:
		try:
			cursor.execute("DELETE from blast where filename=\'" + filename +"\'")
		except mdb.Error, e:
			sys.exit( "# ERROR blast delete filename %s %d: %s" % ( filename, e.args[0],e.args[1] ) )
		status = 0	

	return status	

# END clear_database_rows
#------------------------------------------------------------------

#------------------------------------------------------------------
# Process lines and load if specified

def process_lines(file_lines, filename):
	#------------------------------------------------------------------
	# Initiate variables
	#------------------------------------------------------------------
	values_to_insert = []
	counter = 0
	totalcounter = 0
	loadstamp = time.strftime("%c")

	#------------------------------------------------------------------
	# Get values for each row
	#------------------------------------------------------------------

	for line in file_lines:
		line = line.rstrip()

		# Always check that the file has the format which is expected - runs for each line
		if not len( line.split( "\t" ) ) == 11:	
			sys.exit( "# ERROR: lines do not have the right number of elements ( 11 ): \
				qseqid, sseqid, pident, qlen, qstart, qend, slen, sstart, send, evalue, bitscore" )

		counter += 1		# Increment counters
		totalcounter += 1	# Increment counters
		values = ()			# Empty list of values

		( q_seqid, h_seqid, pident, q_len, q_start, q_end, h_len, h_start, h_end, evalue, bitscore ) = line.split( "\t" )[0:11]

		# Remove selfhits
		if q_seqid == h_seqid and totalcounter != len( file_lines ):
			continue
		 

		else:		#print q_seqid, h_seqid
			(h_org, h_seqkey, h_tail, q_org, q_seqkey, q_tail) = ("","","","","","")


			if len(h_seqid.split("|")) >= 4 : 	(h_org, h_seqkey, h_tail) 	= h_seqid.split("|")[1:4]
			elif len(h_seqid.split("_")) >= 2:	(h_org, h_seqkey) 			= h_seqid.split("_")[0:2]
			elif re.search("Afu", h_seqid) :	(h_org, h_seqkey) 			= ("sApfu1", h_seqid)
			else: 								(h_org, h_seqkey, h_tail) 	= ( h_seqid, "","")

			if len(q_seqid.split("|")) >= 4 :	(q_org, q_seqkey, q_tail) 	= q_seqid.split("|")[1:4]
			elif len(q_seqid.split("_")) >= 2:	(q_org, q_seqkey) 			= q_seqid.split("_")[0:2]
			elif re.search("Afu", q_seqid) :	(q_org, q_seqkey) 			= ("Aspfu1", q_seqid)
			else: 								(q_org, q_seqkey, q_tail) 	= ( q_seqid, "","")

			if re.search("^Afoetidus", str(filename)):
				q_org = "Afoetidus"
			if re.search("VsAfoetidus", str(filename)):
				h_org = "Afoetidus"
			
			#print filename, q_org, h_org
				
			values = ( filename, q_seqid, q_org, q_seqkey, q_tail , h_seqid, h_org, h_seqkey, h_tail,\
					pident, q_len, q_start, q_end, h_len, h_start, h_end, evalue, bitscore , loadstamp)

			values = [re.sub( r"^0+", "", i ) for i in values]
			values_to_insert.append( values ) # Create sets of lists of values


		# Only load data if action is specified
		if args.action == "load" and ( counter == 6000 or totalcounter == len( file_lines ) ):
			#print values_to_insert[0]
			if totalcounter % 1000 == 0 : print "# INFO: Inserting record number %s" % totalcounter
			elif totalcounter == len(file_lines): print "# INFO: Inserting record number %s" % totalcounter

			try:	
				cursor.executemany( "INSERT INTO blast ( filename, q_seqid, q_org, q_seqkey, q_tail ,\
									h_seqid, h_org, h_seqkey, h_tail, pident, q_len, q_start, q_end,\
									h_len, h_start, h_end, evalue, bitscore , loadstamp )\
									values(%s);" % ("%s," * len(values)).rstrip(","),values_to_insert )
				db.commit()	# Add changes to database
				
				counter = 0 			# restart counter
				values_to_insert = []	# Empty list of values

			except mdb.Error, e:
				sys.exit( "# ERROR blast load %s %d: %s" % ( filename, e.args[0],e.args[1] ) )
		
	if args.action == "test" and values:
		print "# INFO: testing file, example values:\n#%s" % values
	if args.action == "test" and not values:
		print "# WARNING: testing file, example values is empty"

	print "# FINISHED BLAST file: %s with total number of records %s" % ( filename , totalcounter )

# END process_lines
#------------------------------------------------------------------


#------------------------------------------------------------------
# Get lines from file
#------------------------------------------------------------------
status = 0
for source in allpath:
	if not os.path.isfile(source):
		sys.exit("# ERROR: Source is not a file %s" % source)
	print "#---------------------------------------------"	
	print "# INFO: Processing file %s" % source
	print "# INFO: Getting lines from file"
	file_lines = get_lines(source)

	filename = source.split("/")[-1]
	print "# INFO: Specified filename %s" % filename

	print "# INFO: Clearing rows with same filename "
	status = clear_database_rows(filename)

	if status == 0:
		print "# INFO: Processing and loading/testing lines"
		process_lines(file_lines, filename)

db.close()	# Close DB connection

print "# INFO Runtime: ", (datetime.now()-startTime)

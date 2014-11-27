#!/usr/bin/python

#------------------------------------------------------------------
# Imports
#------------------------------------------------------------------
import MySQLdb as mdb
import sys, getopt, argparse, re, glob, os, gzip

from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from mysql_utils_data import * # custom functions
from utilsArgparse import * # custom functions
from datetime import datetime

##############################################################################################################################################
##############################################################################################################################################
####################################################################### ARGUMENTS ############################################################
##############################################################################################################################################
##############################################################################################################################################
startTime = datetime.now() # record runtime


#------------------------------------------------------------------
# Get command line arguments
#------------------------------------------------------------------
parser = CustomArgumentParser(formatter_class=SmartFormatter, 
								description="Load data into local mysql database", 
								 usage='%(prog)s -filetype filetype -action action -source filename\n'
								 "Example: python %(prog)s -dbname aspminedb -filetype uma -action test -source Aspsac1_AssemblyScaffolds.fasta.gz\n"
								 "Example: python %(prog)s -dbname aspminedb -filetype dir -action test -source /home/user/Documents/Aspergillus_flavus_NRRL3357_[Aspfl1]")
#------------------------------------------------------------------
# Choose if data should be loaded or only tested
#------------------------------------------------------------------
parser.add_argument('-action' , "-a", required=True, default = "test", choices=['load', 'test'],
					help='R|Test or test and load data\n'
					"-test\t: Read files and verify formats\n"
					"-load\t: Read files and verify formats followed by storing data in MySQL tables\n")

#------------------------------------------------------------------
# At least one file filetype must be provided, or the option "all" which expects a JGI directory format
#------------------------------------------------------------------
parser.add_argument('-filetype' ,"-f",  required=True, default = "ma", choices=['pf', 'tf', 'cf', 'uma', 'ma', 'gff', 'all', "dir", "kog","kegg","ipr","go","sigp", "pm", "blast", "anti"],
					help="R|Select filetype or complete JGI directory, some files in .gz format (*).\n"
					"-filetype pf\t: Protein FASTA file*\n"
					"-filetype tf\t: Transcript FASTA file*\n"
					"-filetype cf\t: CDS FASTA file*, all\n"
					"-filetype uma\t: Un-masked assembly FASTA file*\n"
					"-filetype ma\t: Masked assembly file*\n"
					"-filetype gff\t: General Feature Format file*\n"
					"-filetype kog\t: KOG annotation file*\n"
					"-filetype kegg\t: KEGG annotation file*\n"
					"-filetype ipr\t: Interpro annotation file*\n"
					"-filetype sigp\t: SignalP annotation file*\n"
					"-filetype go\t: GO annotation file*\n"
					"-filetype dir\t: complete JGI directory structure, see README for further instructions, all files in .gz\n" 
					"-filetype pm\t: primary metabolism model (local production, NOT JGI)\n"
					"-filetype blast\t: BLAST, all against all comparisons (local production, NOT JGI)\n"
					"-filetype anti\t: Antismash annotation, JGI run but not part of directory structure")

parser.add_argument("-source", "-s", required=True, help="R|File or directoryname, example: -source Aspbr1_AssemblyScaffolds.fasta.gz")#, type=argparse.FileType('r'))
parser.add_argument("-dbname", "-d", required=False, default = "aspminedb", help="R|Database name")
parser.add_argument("-prefix", "-p", required=False, help="R|Organism prefix")
parser.add_argument("-orgkey", "-o",required=False, help="R|Organism key/name")
parser.add_argument("-verbose", "-v", required=False, help="R|Increase output verbosity for bug tracking" , action='store_true')

args = parser.parse_args()

#------------------------------------------------------------------
# Check prefix validity 
# user can define a prefix to force the duplicated entering of an organism in the database
#------------------------------------------------------------------
if args.prefix and not re.search(r"^[0-9A-Fa-f]*$", args.prefix):
	sys.exit("# ERROR : prefix contains characters other than letters and numbers")

#------------------------------------------------------------------
# Print argument values to screen
#------------------------------------------------------------------
print "#--------------------------------------------------------------"
print "# ARGUMENTS :\n# source\t:%s\n# filetype\t:%s\n# action\t:%s\n# prefix\t:%s\n# Database\t:%s\n# Organism key\t:%s" % (args.source, args.filetype, args.action, args.prefix, args.dbname, args.orgkey)
print "#--------------------------------------------------------------"

#------------------------------------------------------------------
# Test functions
#------------------------------------------------------------------
if args.verbose:
	test__get_org_id()
	test__fasta()

	#test__gff()
	#test__ipr()
	#test__kegg()
	#test__go()
	#test__sigp()
	#test__kog()

	sys.exit("# FINISHED testing functions")

#------------------------------------------------------------------
# Connect to specific DB
#------------------------------------------------------------------
try:
    db = mdb.connect("localhost","asp","1234",str(args.dbname))
    print "# INFO: success connectiong to existing database %s" % args.dbname

except mdb.Error, e:
	print "# INFO: creating new database %s" % args.dbname
	print "# INFO: loading BLAST table will take hours!"
	try:
		db = mdb.connect("localhost","root","1234")
		cursor = db.cursor()
		sql = 'CREATE DATABASE %s' % str(args.dbname)
		cursor.execute(sql)	
	except mdb.Error, e:
		sys.exit("# ERROR %d: %s" % (e.args[0],e.args[1]))

	try:
		db = mdb.connect("localhost","asp","1234",str(args.dbname))
	except mdb.Error, e:
		sys.exit("# ERROR %d: %s" % (e.args[0],e.args[1]))

	cursor = db.cursor()
	myfile = open ("/home/tcve/github/db_mysql/aspminedb.sql", "r") 
	data = ''.join( [line.replace('\n', '') for line in myfile.readlines()] )
	data = data.split(";")

	cursor.execute("SET FOREIGN_KEY_CHECKS = 0;")
	for d in data:
		cursor.execute(d)
		db.commit()
	cursor.execute("SET FOREIGN_KEY_CHECKS = 1;")
	db.commit()
#sys.exit()
##############################################################################################################################################
##############################################################################################################################################
####################################################################### Main #################################################################
##############################################################################################################################################
##############################################################################################################################################
source = args.source
filetype = args.filetype
action = args.action
dbname = args.dbname


org_name = ""
org_id = ""

# Define organism key/name
if args.orgkey:
	org_name = args.orgkey
	print "# INFO: organism key/name was determined from argument: ", org_name
elif re.search("\[(.+)\]", str(args.source)):
	org_name = (re.search("\[(.+)\]", str(args.source))).group(1)
	print "# INFO: organism key/name was determined from path: ", org_name
else:
	parser.print_help()
	sys.exit("# ERROR: organism key was not defined and could not be determined from path")


# Load entire JGI directory structure
if filetype == "dir":
	
	# Verify source 
	source_exists("dir", source)
	
	print "#--------------------------------------------------------------"
	print "# INFO: reading sequence data files"
	print "#--------------------------------------------------------------"

	# Create file of complete file list for the project
	filelist = open(org_name+'_filelist.txt', 'w') 

	for root, subfold, files in os.walk(source):
		for filename in files: 
			path = root + "/" + filename

			# Write all file names to file list for future reference	
			filelist.write(path)
			filelist.write("\n")
			
			# Read genecatalog files, fasta and gff
			# These structures will try to automatically deside what data is found in each file
			# assumption: all protein fasta files countain the word "proteins" and "fasta" in the file name
			# The same sort of assumptions are made for transcripts, cds and assemblies
			if not re.search(r"fasta", filename) and re.search(r"gz$", filename): continue
			
			if re.search("_(GeneCatalog|filtered)_transcripts", filename):
				org_id = fasta(org_name,path, action, "tf", dbname, args.prefix)

			elif re.search("_(GeneCatalog|filtered)_proteins", filename):
				org_id = fasta(org_name,path, action, "pf", dbname, args.prefix)

			elif re.search("_(GeneCatalog|filtered)_CDS", filename):
				org_id = fasta(org_name,path, action, "cf", dbname, args.prefix)

			elif re.search(r"Assembly", filename) or re.search(r"_masked_scaffold", filename):
				if not re.search(r"unmasked", str(path)):
					org_id = fasta(org_name,path, action, "ma", dbname, args.prefix)

	filelist.close()
	
	print "#--------------------------------------------------------------"
	print "# INFO: reading annotation files"
	print "#--------------------------------------------------------------"

	for root, subfold, files in os.walk(source):
		for filename in files: 
			path = root + "/" + filename
			if not re.search("(GeneCatalog|filtered)", filename) and re.search(r"gz$", filename): continue
			
			# Read genecatalog files, gff and annotations
			if re.search(r"gff\.", filename):
				org_id = gff(org_name,path, action, dbname)

			elif re.search(r"_KOG", filename):
				org_id = kog(org_name,path, action, dbname)

			elif re.search(r"_SigP", filename):
				org_id = sigp(org_name,path, action, dbname)

			elif re.search(r"_GO", filename):
				org_id = go(org_name,path, action, dbname)

			elif re.search(r"_KEGG", filename):
				org_id = kegg(org_name,path, action, dbname)

			elif re.search(r"_IPR", filename):
				org_id = ipr(org_name,path, action, dbname)

#------------------------------------------------------------------
# Load individual files
#------------------------------------------------------------------

else: 	
	# Verify source 
	source_exists("file", source)
		
	if filetype == "uma" or filetype == "ma" or filetype == "pf" or filetype == "tf" or filetype == "cf":
		org_id = fasta(org_name,source, action, filetype, dbname, args.prefix)
		
	elif filetype == "gff":
		org_id = gff(org_name,source, action, dbname)

	elif filetype == "kog":
		org_id = kog(org_name,source, action, dbname)

	elif filetype == "sigp":
		org_id = sigp(org_name,source, action, dbname)

	elif filetype == "go":
		org_id = go(org_name,source, action, dbname)

	elif filetype == "kegg":
		org_id = kegg(org_name,source, action, dbname)

	elif filetype == "ipr":
		org_id = ipr(org_name,source, action, dbname)
	else:
		sys.exit("# ERROR: filetype has no defined action")

#------------------------------------------------------------------
# At this point the org_id should have been defined by one or more functions
# If this has not happened something went wrong
#------------------------------------------------------------------
print "#--------------------------------------------------------------"
print "# INFO: Final error handling"
print "#--------------------------------------------------------------"

if org_id:
	print "# INFO Organism id: ", org_id
else:
	sys.exit("# ERROR: something went wrong, organism ID was not returned, filetype = %s" % source)
	
#------------------------------------------------------------------
# Print database content for organism
#------------------------------------------------------------------
db = mdb.connect("localhost","asp","1234",dbname)
cursor = db.cursor()
cursor.execute("SELECT organism.org_id, name, ( SELECT count(*) from proteins where org_id=%s), ( SELECT count(*) from transcripts where org_id=%s ), ( SELECT count(*) from assembly where org_id=%s ), ( SELECT count(*) from cds where org_id=%s ) from organism where org_id = %s" % (org_id,org_id,org_id,org_id,org_id))
results = cursor.fetchall()
print "# INFO database data: " , (results)

if action == "load":
	if results[0][2] == 0:
		print "# WARNING: no protein data has been loaded for this project"
	if results[0][3] == 0:
		print "# WARNING: no transcript data has been loaded for this project"
	if results[0][4] == 0:
		print "# WARNING: no assembly data has been loaded for this project"
	if results[0][5] == 0:
		print "# WARNING: no cds data has been loaded for this project"

#------------------------------------------------------------------
# Print program runtime
#------------------------------------------------------------------
print "# INFO Runtime: ", (datetime.now()-startTime)
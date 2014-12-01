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
#from mysql_classes import Fasta	# custom class
#from mysql_classes import Gff	# custom class
from mysql_utils import * # custom functions
import pprint

##############################################################################################################################################
##############################################################################################################################################
####################################################################### ARGUMENTS ############################################################
##############################################################################################################################################
##############################################################################################################################################
from datetime import datetime

startTime = datetime.now()

#------------------------------------------------------------------
# Format argument helper
#------------------------------------------------------------------
class SmartFormatter(argparse.HelpFormatter):
	width = 100
	def _split_lines(self, text, width):
	#------------------------------------------------------------------
	# this is the RawTextHelpFormatter._split_lines
	#------------------------------------------------------------------
		if text.startswith('R|'):
			return text[2:].splitlines()  
		return argparse.HelpFormatter._split_lines(self, text, width)

#------------------------------------------------------------------
# Get command line arguments
#------------------------------------------------------------------

parser = argparse.ArgumentParser(formatter_class=SmartFormatter, 
								 usage='%(prog)s -filetype filetype -action action -source filename\n'
								 "Example: python %(prog)s -dbname aspminedb -filetype uma -action test -source Aspsac1_AssemblyScaffolds.fasta.gz\n"
								 "Example: python %(prog)s -dbname aspminedb -filetype dir -action test -source /home/user/Documents/Aspergillus_flavus_NRRL3357_[Aspfl1]")

#------------------------------------------------------------------
# Choose if data should be loaded or only tested
#------------------------------------------------------------------
parser.add_argument('-action' , required=True, default = "test", choices=['load', 'test'],
					help='R|Test or test and load data\n'
					"-test\t: Read files and verify formats\n"
					"-load\t: Read files and verify formats followed by storing data in MySQL tables\n")

#------------------------------------------------------------------
# At least one file filetype must be provided, or the option "all" which expects a JGI directory format
#------------------------------------------------------------------
parser.add_argument('-filetype' , required=True, default = "ma", choices=['pf', 'tf', 'cf', 'uma', 'ma', 'gff', 'all', "dir", "kog","kegg","ipr","go","sigp"],
					help="R|Select filetype or complete JGI directory, all files in .gz format.\n"
					"-filetype pf\t: Protein FASTA file\n"
					"-filetype tf\t: Transcript FASTA file\n"
					"-filetype cf\t: CDS FASTA file, all\n"
					"-filetype uma\t: Un-masked assembly FASTA file\n"
					"-filetype ma\t: Masked assembly file\n"
					"-filetype gff\t: General Feature Format file\n"
					"-filetype kog\t: KOG annotation file\n"
					"-filetype kegg\t: KEGG annotation file\n"
					"-filetype ipr\t: Interpro annotation file\n"
					"-filetype sigp\t: SignalP annotation file\n"
					"-filetype go\t: GO annotation file\n"
					"-filetype dir\t: complete JGI directory structure, see README for further instructions ")

parser.add_argument("-source", required=True, help="R|File or directoryname, example: -source Aspbr1_AssemblyScaffolds.fasta.gz")#, type=argparse.FileType('r'))
parser.add_argument("-dbname", required=False, default = "aspminedb", help="R|Database name")
parser.add_argument("-prefix", required=False, help="R|Organism prefix")
parser.add_argument("-orgkey", required=False, help="R|Organism key/name")

args = parser.parse_args()

#------------------------------------------------------------------
# Check prefix validity
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
# Connect to specific DB
#------------------------------------------------------------------
try:
    db = mdb.connect("localhost","asp","1234",str(args.dbname))

except mdb.Error, e:
    sys.exit("# ERROR %d: %s" % (e.args[0],e.args[1]))

	
##############################################################################################################################################
##############################################################################################################################################
####################################################################### Main #################################################################
##############################################################################################################################################
##############################################################################################################################################
org_id = ''
source = args.source
filetype = args.filetype
action = args.action
dbname = args.dbname
allfiles = []
org_name = ""

if args.orgkey:
	org_name = arg.orgkey
elif re.search("\[(.+)\]", str(args.source)):
	org_name = (re.search("\[(.+)\]", str(args.source))).group(1)
else:
	sys.exit("# ERROR: organism key was not defined and could not be determined from path")

print "# INFO: organism key/name was determined from path: ", org_name

#------------------------------------------------------------------
# Load entire JGI directory structure
#------------------------------------------------------------------
if filetype == "dir":
	
	#------------------------------------------------------------------
	# Make sure that the path exists
	#------------------------------------------------------------------
	if os.path.exists(source):	print "# INFO: path exists: %s" % source 
	else:	sys.exit("# ERROR: path does not exists" )
	
	#------------------------------------------------------------------
	# Loop through all files in directory
	#------------------------------------------------------------------
	for root, subfold, files in os.walk(source):
		allfiles.append(files)
		
		#------------------------------------------------------------------
		# Exit if no files are found
		#------------------------------------------------------------------
		if allfiles == []:
			sys.exit("# ERROR: directory is empty")
		
	print "#--------------------------------------------------------------"
	print "# INFO: reading basic data files"
	print "#--------------------------------------------------------------"

	#------------------------------------------------------------------
	# Create file of complete file list for the project
	#------------------------------------------------------------------
	filelist = open(org_name+'_filelist.txt', 'w') 
	
	for root, subfold, files in os.walk(source):
		for f in files: 
			path = root + "/" + f
			filelist.write(path)
			filelist.write("\n")
			
			if not re.search(r"fasta", str(f)) and re.search(r"gz$", str(f)): continue
			
			#------------------------------------------------------------------
			# Read genecatalog files, fasta and gff
			# These structures will try to automatically deside what data is found in each file
			# assumption: all protein fasta files countain the word "proteins" and "fasta" in the file name
			# The same sort of assumptions are made for trasncripts, cds and assemblies
			#------------------------------------------------------------------
			if re.search("_(GeneCatalog|filtered)_transcripts", str(f)):
				org_id = fasta(org_name,path, action, "tf", dbname, args.prefix)

			elif re.search("_(GeneCatalog|filtered)_proteins", str(f)):
				org_id = fasta(org_name,path, action, "pf", dbname, args.prefix)

			elif re.search("_(GeneCatalog|filtered)_CDS", str(f)):
				org_id = fasta(org_name,path, action, "cf", dbname, args.prefix)

			elif re.search(r"Assembly", str(f)) or re.search(r"_masked_scaffold", str(f)):
				if not re.search(r"unmasked", str(path)):
					org_id = fasta(org_name,path, action, "ma", dbname, args.prefix)
	filelist.close()
	
	print "#--------------------------------------------------------------"
	print "# INFO: reading annotation files"
	print "#--------------------------------------------------------------"

	for root, subfold, files in os.walk(source):
		for f in files: 
			path = root + "/" + f
			if not re.search("(GeneCatalog|filtered)", str(f)) and re.search(r"gz$", str(f)): continue
			
			#------------------------------------------------------------------
			# Read genecatalog files, fasta and gff
			#------------------------------------------------------------------
			if re.search(r"gff\.", str(f)):
				org_id = gff(org_name,path, action, dbname)

			elif re.search(r"_KOG", str(f)):
				org_id = kog(org_name,path, action, dbname)

			elif re.search(r"_SigP", str(f)):
				org_id = sigp(org_name,path, action, dbname)

			elif re.search(r"_GO", str(f)):
				org_id = go(org_name,path, action, dbname)

			elif re.search(r"_KEGG", str(f)):
				org_id = kegg(org_name,path, action, dbname)

			elif re.search(r"_IPR", str(f)):
				org_id = ipr(org_name,path, action, dbname)

#------------------------------------------------------------------
# Load individual files
#------------------------------------------------------------------

else: 	
	if os.path.isfile(source):
		print "# INFO: file exists: %s" % source 
	else:
		sys.exit("# ERROR: file does not exists" )
		
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
		sys.exit("ERROR: filetype has no defined action")

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

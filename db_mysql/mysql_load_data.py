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
from mysql_utils_data import * # custom functions
import pprint
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


if args.verbose:
	test()
	sys.exit

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
# Connect to specific DB
#------------------------------------------------------------------
try:
    db = mdb.connect("localhost","asp","1234",str(args.dbname))

except mdb.Error, e:
    sys.exit("# ERROR %d: %s" % (e.args[0],e.args[1]))

sys.exit()	
##############################################################################################################################################
##############################################################################################################################################
####################################################################### Main #################################################################
##############################################################################################################################################
##############################################################################################################################################
source = args.source
filetype = args.filetype
action = args.action
dbname = args.dbname

allfiles = []
org_name = ''
org_id = ''

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
	print "# INFO: reading sequence data files"
	print "#--------------------------------------------------------------"

	#------------------------------------------------------------------
	# Create file of complete file list for the project
	#------------------------------------------------------------------
	filelist = open(org_name+'_filelist.txt', 'w') 
	"""
	for root, subfold, files in os.walk(source):
		for filename in files: 
			path = root + "/" + f
			filelist.write(path)
			filelist.write("\n")
			
			if not re.search(r"fasta", filename) and re.search(r"gz$", filename): continue
			
			#------------------------------------------------------------------
			# Read genecatalog files, fasta and gff
			# These structures will try to automatically deside what data is found in each file
			# assumption: all protein fasta files countain the word "proteins" and "fasta" in the file name
			# The same sort of assumptions are made for trasncripts, cds and assemblies
			#------------------------------------------------------------------
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
	"""
	print "#--------------------------------------------------------------"
	print "# INFO: reading annotation files"
	print "#--------------------------------------------------------------"

	for root, subfold, files in os.walk(source):
		for filename in files: 
			path = root + "/" + filename
			if not re.search("(GeneCatalog|filtered)", filename) and re.search(r"gz$", filename): continue
			
			#------------------------------------------------------------------
			# Read genecatalog files, fasta and gff
			#------------------------------------------------------------------
			"""
			if re.search(r"gff\.", filename):
				org_id = gff(org_name,path, action, dbname)

			elif re.search(r"_KOG", filename):
				org_id = kog(org_name,path, action, dbname)

			elif re.search(r"_SigP", filename):
				org_id = sigp(org_name,path, action, dbname)
			"""
			if re.search(r"_GO", filename):
				org_id = go(org_name,path, action, dbname)
			"""
			elif re.search(r"_KEGG", filename):
				org_id = kegg(org_name,path, action, dbname)

			elif re.search(r"_IPR", filename):
				org_id = ipr(org_name,path, action, dbname)
			"""
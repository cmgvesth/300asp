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
####################################################################### TEST ############################################################
##############################################################################################################################################
##############################################################################################################################################
def test():
	(org_name, dbname, prot_id) = ("Aspni1", "aspminedb", "19949")

	# Test get_org_id
	org_id = get_org_id(org_name, dbname, 1)
	if not (isinstance( org_id, long ) or  isinstance( org_id, int )):
		sys.exit("# ERROR: get_org_id does not return an integer number")

	# Create temporary dna FASTA file
	fasta_dna_entry = ">scaffold_1\nTTTTTAATTTTCGAACTAACTAAAAAAACTTTTATTTAATCTGATAAAAA\nAAACCAACtattattatattaaattctaattatatattttatttagttta\nttaattttttaataatattaCAGAGAAGAGCAGTAATTACTGTAGACtta\n>scaffold_7\ntaatgtaaaccctaatgtaaaccctaatgtaaaccctaatgtaaacccta\natgtaaaccctaatgtaaccctaacccccctaatctgacaaccctaaTGT\nAATACCTTCGGCGGAAGGTTTATCCGGTGCTGACGGGCCGAGGTTTCTCC"
	#fasta_dna_entry = ">jgi|Aspac1|48858|fgenesh1_pm.1_#_1\nMSPAFRIQAISIPNNYLPPRKSFTPLLSPLSPSLRLFGRYSLHRHLLVELGEWEAAYEKASAFVADLTTA\nRGFERVEIDSGGHSYVIFDLRRRDLSYWDVMAQEWLVAPGSYNVFIGASSRDLRLNGTFSLTVTA*\n>jgi|Aspac1|37826|Genemark1.2_g\nMDCDETKPECSNCVNHSVRCIYDSPAPKTKKPVTLPAQTHKTENKVQVQFVEFNFLHVQPAPAPVPATPA\nRLEQPFWEALQRKQPVALMLVAHFAALMRWMDGTWVSVGWPEHILRGAWMFLEEEERGLVRWLIDDVGVL\nGLSG*\n>jgi|Aspac1|37827|Genemark1.3_g\nMVRGIIFQLVFVGLVVDFVLRLSKRGRPQALWSNRPLLLLGGATALSLLLIYIRSVYRTIELLHGWTSST\nMHNEMLLIGLDGAIMVPAFSVYNLLHPGYLLPKVQREVGYLDARGLQMEMEAVKDGRYQIVEVEGGKDDS\nVTLLDRS*"
	tmp_dna_filename = 'tmp_' + org_name + '_dna.fsa.gz'
	tmp_dna_file = gzip.open(tmp_dna_filename, 'wb') 
	tmp_dna_file.write(fasta_dna_entry)
	tmp_dna_file.close()

	# Create temporary gene dna FASTA file
	fasta_gene_entry = ">jgi|Aspac1|48858|fgenesh1_pm.1_#_1\nATGTCCCCTGCCTTTCGCATTCAAGCTATCTCGATACCCAATAATTACTTGCCTCCACGCAAGTCGTTCA\nATCTCTCATATTGGGACGTCATGGCTCAGGAGTGGCTCGTCGCTCCCGGCTCATATAATGTCTTCATTGG\nCGCGAGCTCACGCGACCTAAGACTGAATGGAACGTTTTCTTTGACTGTTACAGCCTAA\n>jgi|Aspac1|37826|Genemark1.2_g\nATGGATTGCGACGAGACAAAGCCCGAGTGCTCGAACTGCGTCAACCATTCTGTGCGCTGTATCTACGATT\nTGCGTGGATGTTTTTGGAGGAGGAGGAGCGGGGGTTGGTGAGATGGTTGATTGACGATGTCGGGGTGTTG\nGGGTTGTCTGGTTGA\n>jgi|Aspac1|37827|Genemark1.3_g\nATGGTCCGGGGTATTATCTTCCAGCTCGTATTCGTGGGACTGGTCGTCGACTTCGTGCTGCGCCTGTCCA\nGTTACGTTGCTGGACAGGAGCTGA\n"
	#fasta_dna_entry = ">jgi|Aspac1|48858|fgenesh1_pm.1_#_1\nMSPAFRIQAISIPNNYLPPRKSFTPLLSPLSPSLRLFGRYSLHRHLLVELGEWEAAYEKASAFVADLTTA\nRGFERVEIDSGGHSYVIFDLRRRDLSYWDVMAQEWLVAPGSYNVFIGASSRDLRLNGTFSLTVTA*\n>jgi|Aspac1|37826|Genemark1.2_g\nMDCDETKPECSNCVNHSVRCIYDSPAPKTKKPVTLPAQTHKTENKVQVQFVEFNFLHVQPAPAPVPATPA\nRLEQPFWEALQRKQPVALMLVAHFAALMRWMDGTWVSVGWPEHILRGAWMFLEEEERGLVRWLIDDVGVL\nGLSG*\n>jgi|Aspac1|37827|Genemark1.3_g\nMVRGIIFQLVFVGLVVDFVLRLSKRGRPQALWSNRPLLLLGGATALSLLLIYIRSVYRTIELLHGWTSST\nMHNEMLLIGLDGAIMVPAFSVYNLLHPGYLLPKVQREVGYLDARGLQMEMEAVKDGRYQIVEVEGGKDDS\nVTLLDRS*"
	tmp_gene_filename = 'tmp_' + org_name + '_gene.fsa.gz'
	tmp_gene_file = gzip.open(tmp_gene_filename, 'wb') 
	tmp_gene_file.write(fasta_gene_entry)
	tmp_gene_file.close()

	# Create temporary protein FASTA file
	fasta_protein_entry = ">jgi|Aspac1|43358|fgenesh1_pm.1_#_1\nMSPAFRIQAISIPNNYLPPRKSFTPLLSPLSPSLRLFGRYSLHRHLLVELGEWEAAYEKASAFVADLTTA\nRGFERVEIDSGGHSYVIFDLRRRDLSYWDVMAQEWLVAPGSYNVFIGASSRDLRLNGTFSLTVTA*\n>jgi|Aspac1|37826|Genemark1.2_g\nMDCDETKPECSNCVNHSVRCIYDSPAPKTKKPVTLPAQTHKTENKVQVQFVEFNFLHVQPAPAPVPATPA\nRLEQPFWEALQRKQPVALMLVAHFAALMRWMDGTWVSVGWPEHILRGAWMFLEEEERGLVRWLIDDVGVL\nGLSG*\n>jgi|Aspac1|37827|Genemark1.3_g\nMVRGIIFQLVFVGLVVDFVLRLSKRGRPQALWSNRPLLLLGGATALSLLLIYIRSVYRTIELLHGWTSST\nMHNEMLLIGLDGAIMVPAFSVYNLLHPGYLLPKVQREVGYLDARGLQMEMEAVKDGRYQIVEVEGGKDDS\nVTLLDRS*"
	#fasta_protein_entry = ">scaffold_1\nTTTTTAATTTTCGAACTAACTAAAAAAACTTTTATTTAATCTGATAAAAA\nAAACCAACtattattatattaaattctaattatatattttatttagttta\nttaattttttaataatattaCAGAGAAGAGCAGTAATTACTGTAGACtta\n>scaffold_7\ntaatgtaaaccctaatgtaaaccctaatgtaaaccctaatgtaaacccta\natgtaaaccctaatgtaaccctaacccccctaatctgacaaccctaaTGT\nAATACCTTCGGCGGAAGGTTTATCCGGTGCTGACGGGCCGAGGTTTCTCC"
	tmp_protein_filename = 'tmp_' + org_name + '_protein.fsa.gz'
	tmp_protein_file = gzip.open(tmp_protein_filename, 'wb') 
	tmp_protein_file.write(fasta_protein_entry)
	tmp_protein_file.close()

	# Test different versions of FASTA function
	action = "test"
	org_id = fasta(org_name, tmp_gene_filename, action, "tf", dbname, args.prefix)
	org_id = fasta(org_name, tmp_protein_filename, action, "pf", dbname, args.prefix)
	org_id = fasta(org_name, tmp_gene_filename, action, "cf", dbname, args.prefix)
	org_id = fasta(org_name, tmp_dna_filename, action, "ma", dbname, args.prefix)

	# Test annotation functions
	# Create temporary GFF TAB file
	# Create temporary KOG TAB file
	# Create temporary SignalP TAB file
	# Create temporary GO TAB file
	# Create temporary KEGG TAB file
	# Create temporary Interpro TAB file
	"""
	org_id = gff 	(org_name, path, action, dbname)
	org_id = kog 	(org_name, path, action, dbname)
	org_id = sigp 	(org_name, path, action, dbname)
	org_id = go 	(org_name, path, action, dbname)
	org_id = kegg 	(org_name, path, action, dbname)
	org_id = ipr 	(org_name, path, action, dbname)
	"""
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
	test()
	sys.exit("# FINISHED testing functions")

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
source = args.source
filetype = args.filetype
action = args.action
dbname = args.dbname

allfiles = []
org_name = ""
org_id = ""

# Define organism key/name
if args.orgkey:
	org_name = arg.orgkey
elif re.search("\[(.+)\]", str(args.source)):
	org_name = (re.search("\[(.+)\]", str(args.source))).group(1)
else:
	sys.exit("# ERROR: organism key was not defined and could not be determined from path")

print "# INFO: organism key/name was determined from path: ", org_name


# Load entire JGI directory structure
if filetype == "dir":
	
	# Make sure that the path exists
	if os.path.exists(source):	print "# INFO: path exists: %s" % source 
	else:	sys.exit("# ERROR: path does not exists" )
	
	# Loop through all files in directory and check if it is empty
	for root, subfold, files in os.walk(source):
		allfiles.append(files)
		
	if allfiles == []:
		sys.exit("# ERROR: directory is empty")
		
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
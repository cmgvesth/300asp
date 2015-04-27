#!/usr/bin/python

'''------------------------------------------------------------------
# Imports
------------------------------------------------------------------'''
import sys,os
root = os.path.abspath(__file__).split("github")[0]
sys.path.append(root+"github/utils/")
from aspmine_imports import *

'''------------------------------------------------------------------
# Code uses local functions: 
------------------------------------------------------------------

CustomArgumentParser
DBconnect
executeQuery

------------------------------------------------------------------'''

'''------------------------------------------------------------------
# ARGUMENTS and setup
------------------------------------------------------------------'''
startTime = datetime.now() # record runtime
parser = CustomArgumentParser(formatter_class=SmartFormatter, usage='%(prog)s -dbname [database name]')
parser.add_argument("-dbname", "-d", required=False, default = "aspminedb", help="Database name")
parser.add_argument("-afolder", "-af",  required=False, default="afolder_clusterPresenceAbsence", help="Name of folder for analysis results")

""" Parameters """
parser.add_argument("-anntype", "-at", required=False, default = "gene", choices=["gene","antismash","interpro", "go", "kegg", "kog"], help="Annotation type")
parser.add_argument("-annpattern", "-ap", required=False, default = "100693", help="Annotation pattern, gene name, interpro domain ")
parser.add_argument("-annorg", "-ao", required=False, default = "Aspacu1", help="Annotation organism, what organism this pattern comes from, Aspacu1")

""" ANALYSIS """
parser.add_argument("-a1", required=False, action='store_true', help="Homologs from gene name")

""" PARSE ARGUMENTS """
args 	= parser.parse_args()
anntype = args.anntype
annpattern = args.annpattern
annorg = args.annorg

startTime = datetime.now() # record runtime

Rpath = "/" + "/".join( os.path.abspath(__file__).split("/")[1:-1] ) + "/"
afolder = os.getcwd() + "/" + args.afolder
afolderString = afolder.replace("/", "\\/")
os.system("mkdir -p " + afolder + "" )
	
'''------------------------------------------------------------------
# Print argument values to screen
------------------------------------------------------------------'''
print "#--------------------------------------------------------------\n\
# ARGUMENTS\n\
# Database\t: %s\n\
# Analysis folder\t: %s\n\
# Rpath\t: %s\n\
# Annotation type, -at [default = 'gene']\t: %s\n\
# Annotation pattern, -ap [default = '100693']\t: %s\n\
# Annotation organism, -ao [default = 'Aspacu1']\t: %s\n\
#--------------------------------------------------------------" % (args.dbname, afolder, Rpath, anntype, annpattern, annorg)

'''------------------------------------------------------------------
# Connection to database
------------------------------------------------------------------'''
cursor = DBconnect("localhost", "asp", "1234", args.dbname) # costum function

'''------------------------------------------------------------------
# Make sure data can be found in database
------------------------------------------------------------------'''
if anntype == "gene":
	queryroot = "SELECT *, prot_seqkey as seqkey, prot_tailkey as id from proteins join organism using (org_id)"
	querytail = "where prot_seqkey=\'%s\' and prot_orgkey=\'%s\'" % (annpattern, annorg)

if anntype == "antismash":
	queryroot = "SELECT *, sm_protein_id as seqkey, sm_short as id from antismash join organism using (org_id)"
	querytail = " where sm_short=\'%s\' and name=\'%s\'" % (annpattern, annorg)

if anntype == "interpro":
	queryroot = "SELECT *, protein_has_ipr.protein_id as seqkey, ipr_desc as id from ipr join protein_has_ipr using (ipr_id) join organism using (org_id)"
	querytail = " where ipr_id=\'%s\' and name=\'%s\'" % (annpattern, annorg)

if anntype == "go":
	queryroot = "SELECT *, protein_has_go.protein_id as seqkey, go_name as id from go join protein_has_go using (go_term_id) join organism using (org_id)"
	querytail = " where go_term_id=\'%s\' and name=\'%s\'" % (annpattern, annorg)

if anntype == "kog":
	queryroot = "SELECT *, protein_has_go.protein_id as seqkey, kog_defline as id from kog join protein_has_go using (kog_id) join organism using (org_id)"
	querytail = " where kog_id=\'%s\' and name=\'%s\'" % (annpattern, annorg)

if anntype == "kegg":
	queryroot = "SELECT *, protein_has_kegg.protein_id as seqkey, kegg_definition as id from kegg join protein_has_kegg using (kegg_id) join organism using (org_id)" 
	querytail = " where kegg_id=\'%s\' and name=\'%s\'" % (annpattern, annorg)

'''------------------------------------------------------------------
# Test basic data
------------------------------------------------------------------'''
print "# INFO: Verifying basic data from gene and annotation tables"

(columns, values) = executeQuery(cursor, queryroot + querytail)

if len(values) > 0:
	print "# DATA: rows were found to match the pattern, row count: %s" % len(values)

else:
	query = "SELECT * from organism where name=\'%s\'" % (annorg)
	(columns, values) = executeQuery(cursor, query)
	if len(values) <= 0:
		sys.exit("# ERROR: no data found for organism" )		
	sys.exit("# ERROR: no data found for pattern")	

print "-----------------------------"
'''------------------------------------------------------------------
# Combine basic data to BLAST table and get homologs
------------------------------------------------------------------'''
print "# INFO: Combining basic data with BLAST table"

startTime_a = datetime.now()

queryblast = "SELECT q_org, q_seqkey, h_org, h_seqkey,pident, ta.id  from (" + queryroot + querytail + ") as ta  \
JOIN blast on ((ta.name=q_org and q_seqkey=ta.seqkey) or (ta.name=h_org and h_seqkey=ta.seqkey)) where q_seqkey > h_seqkey"

(columns, values) = executeQuery(cursor, queryblast)

print "# RUNTIME: queryblast " + str( datetime.now() - startTime_a )

if len(values) > 0:

	print "# DATA: homologs were found to match the pattern, homolog count: %s" % len(values)
	filename = afolder + "/homologs_ap%s_ao%s_at%s_c%s.csv" % (annpattern, annorg, anntype, len(values))

	cursor2csv(columns, values, filename)
	print "# INFO: wrote results to %s" % filename

if len(values) > 2000:
	sys.exit("# ERROR: to many homologs were found to match the pattern, homolog count: %s" % len(values) )

print "-----------------------------"

'''------------------------------------------------------------------
# Plot data in R
------------------------------------------------------------------'''
print "# INFO: running Rscript"
filenameString = filename.replace("/", "\\/")

os.system("sed \"s/replace_folderpath/\'" + afolderString + "\'/g\" " + Rpath + "sequence_blast_clusters.R > " + afolder + "/tmp.R")
os.system("sed -i \"s/replace_infile/\'" + filenameString + "\'/g\" " + afolder + "/tmp.R")
os.system("R CMD BATCH " + afolder + "/tmp.R " + afolder + "/test.out ")

'''------------------------------------------------------------------
------------------------------------------------------------------

ADD: R script
read in data
cast to wide format
make clustering
- complete
- single 
- kmer

plot clustering
- heatmap
- clusplot



------------------------------------------------------------------
------------------------------------------------------------------'''


'''------------------------------------------------------------------
# Finish
------------------------------------------------------------------'''

print "# RUNTIME: " + str( datetime.now() - startTime )

print "-----------------------------"

#!/usr/bin/python

'''##################################################################
# Imports
##################################################################'''
import sys 
root = os.path.abspath(__file__).split("github")[0]
sys.path.append(root+"github/utils/")
from aspmine_imports import *

'''##################################################################
# ARGUMENTS and setup
##################################################################'''
startTime = datetime.now() # record runtime
parser = CustomArgumentParser(formatter_class=SmartFormatter, 
								description="Get database statistics fro each organism", 
								 usage='%(prog)s -dbname [database name] -out [output filename] \n'
								 "Example: python %(prog)s -dbname aspminedb -out dbstats.csv")
parser.add_argument("-out", "-o", required=False, default = "dbstats.csv", help="Name of output file, default hitlist.csv")
parser.add_argument("-dbname", "-d", required=False, default = "aspminedb", help="Database name")

args 	= parser.parse_args()
outfile = args.out

'''------------------------------------------------------------------
# Print argument values to screen
------------------------------------------------------------------'''
print "#--------------------------------------------------------------\n\
# ARGUMENTS\n\
# Database\t: %s\n\
# Outfile\t: %s\n\
#--------------------------------------------------------------" % (args.dbname, outfile)

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
# Get data
------------------------------------------------------------------'''
f = open(outfile,'wb')
writer = csv.writer(f, dialect = 'excel', delimiter=";")
dbstats = ()	

try:
	cursor.execute("SELECT org_id from organism ;")
except mdb.Error, e:
	sys.exit("# ERROR: executing following query failed:\n%s\n # ERROR %d: %s" % (cursor._last_executed, e.args[0],e.args[1]))

org_ids = cursor.fetchall()

for org_id in org_ids:
	#print org_id
	org_id = str(org_id[0])  
  	try:
		cursor.execute("drop table IF EXISTS tmp ;")
		cursor.execute("select * from \
		(select name, real_name, org_id from organism where org_id=" + org_id + ") as o,\
		(select count(distinct(ipr_id)) as ipr from protein_has_ipr where org_id=" + org_id + ") as a ,\
		(select count(distinct(kegg_id)) as kegg from protein_has_kegg where org_id=" + org_id + ") as b ,\
		(select count(distinct(kog_id)) as kog from protein_has_kog where org_id=" + org_id + ") as c, \
		(select count(distinct(go_term_id)) as go from protein_has_go where org_id=" + org_id + ") as d, \
		(select sum(length(assembly_seq)) as dna_seq from assembly where org_id=" + org_id + ") as f, \
		(select count(gff_type) as exon_count from gff where org_id=" + org_id + " and gff_type='exon') as g, \
		(select count(distinct(prot_seqkey)) as prots from proteins where org_id=" + org_id + ") as e;")

	except mdb.Error, e:
		sys.exit("# ERROR: executing following query failed:\n%s\n # ERROR %d: %s" % (cursor._last_executed, e.args[0],e.args[1]))
	columns = map(lambda x:x[0], cursor.description) 	
	dbstats = dbstats + cursor.fetchall()

writer.writerow(columns)

writer.writerows(dbstats)
f.close()

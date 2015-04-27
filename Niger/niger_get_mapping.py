#!/usr/bin/python

'''##################################################################
# Imports
##################################################################'''
import sys,os
root = os.path.abspath(__file__).split("github")[0]
sys.path.append(root+"github/utils/")
from aspmine_imports import *

'''##################################################################
# ARGUMENTS and setup
##################################################################'''
startTime = datetime.now() # record runtime
parser = CustomArgumentParser(formatter_class=SmartFormatter, 
								description="Get mapping for specific set of identifiers.", 
								 usage='%(prog)s -dbname [database name] -id [transciprt or protein ids] -type [prot or trans]\n'
								 "Example: python %(prog)s -dbname aspminedb -id 52305, 534305 -type trans")
parser.add_argument("-dbname", "-d", required=False, default = "aspminedb", help="Database name")
parser.add_argument("-id", required=True, action='store',  default=[None], nargs='*', help="Protein or transcipt identifiers")
parser.add_argument("-type", "-t", required=True, choices=['prot', 'trans'], help="Protein or Transcripts")
parser.add_argument("-out", "-o", required=False, default = "default_mapping.csv", help="Name of output file, default default_mapping.csv")
#parser.add_argument("-strain", "-s", required=True, default=[None], action='store', help="JGI organism keys")
#parser.add_argument("-gene", "-g",  required=True, default=[None], action='store', help="List of gene ID or name")

args 	= parser.parse_args()
ids 	= args.id
dtype 	= args.type
outfile = args.out
'''------------------------------------------------------------------
# Print argument values to screen
------------------------------------------------------------------'''
print "#--------------------------------------------------------------\n\
# ARGUMENTS\n\
# Description\t\t: Script returns a tabble orf mapping of ID numbers between Aniger1, 3 and 7\n\
# Database\t\t: %s\n\
# Type\t\t\t: %s\n\
# IDs\t\t\t: %s\n\
# Outfile\t\t: %s\n\
#--------------------------------------------------------------" % (args.dbname, dtype, ids, outfile)

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



gene_string = ""
for i in range(0,len(ids)):
	gene_string += " query_id like \'%" + ids[i] + "%\' "
	if i != len(ids)-1:
		gene_string += "or"

#print gene_string

try:
	cursor.execute("SELECT * FROM map_niger_" + dtype + "_complete WHERE " + gene_string)
	result = cursor.fetchall()
except mdb.Error, e: 
	sys.exit("# ERROR %d: %s" % (e.args[0],e.args[1]))


f = open(outfile,'wb')
writer = csv.writer(f, dialect = 'excel', delimiter=";")
columns = map(lambda x:x[0], cursor.description) 	
writer.writerow(columns)
writer.writerows(result)
f.close()

"""
(strain_string, gene_string) = ("","")

strain_string = "q_org = \'" + strain + "\' or h_org = \'" + strain + "\'"
gene_string = "q_seqid like \'%" + gene + "%\' or h_seqid like \'%" + gene + "%\'"



	for s in range(0,len(strains)):
		strain_string += " b.q_org = \'" + strains[s] + "\' or b.h_org = \'" + strains[s] + "\'"
		if s != len(strains)-1:
			strain_string += " or "


	for i in range(0,len(genes)):
		gene_string += " a.q_seqid like \'%" + genes[i] + "%\' or a.h_seqid like \'%" + genes[i] + "%\'"
		if i != len(genes)-1:
			gene_string += " or "

"""			
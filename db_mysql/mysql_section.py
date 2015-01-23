#!/usr/bin/python

#------------------------------------------------------------------
# Imports
#------------------------------------------------------------------
import MySQLdb as mdb
import sys, getopt, argparse, re, glob, os
import gzip
import xlrd

"""
CODEBOOK INFO:

Program reads a EXCEL file as exported from Podio, from the Aspergillus genome sequencing workspace.
Export Excel from the Aspergillus species app.
Needed columns include "JGI code", "" and "section"

Created on	Created by	Species	Species alias(es)	Botanical Code name	Strain CBS Acc. No	Strain IBT No	NCBI Organism ID	MycoBank	JGI code	Picture	Section	Project status	DNA sequencing status	Dropbox folder with pictures	Link to genome portal	Reference on species	Genome publication	Description	Tags
2015-01-15 11:55	Mikael Rordam Andersen	acidus	Aspergillus foetidus		CBS 106.47				Aspfo1		Nigri	Previous project/Completed by another lab	Annotation complete		http://genome.jgi.doe.gov/Aspfo1/Aspfo1.home.html				
2015-01-15 10:57	Mikael Rordam Andersen	acristatulus		Emericella acristata	CBS 119.55	11111	41739					Material shipped	Sequencing started						
2015-01-15 9:33	Mikael Rordam Andersen	aculeatinus			121060	29077	487661		Aspacu1	aculeatinus.jpg	Nigri	Project completed	Annotation complete						
2015-01-22 14:52	Mikael Rordam Andersen	aculeatus							Aspac1		Nigri	Previous project/Completed by another lab	Annotation complete						
2015-01-15 10:57	Mikael Rordam Andersen	affinis			CBS 129190	32310	1070780					Material shipped	Sequencing started			http://www.ncbi.nlm.nih.gov/pubmed/21788229			
2015-01-15 10:57	Mikael Rordam Andersen	alabamensis			IBT 12702	12702	657433					Material shipped	Sequencing started						
2015-01-15 10:57	Mikael Rordam Andersen	albertensis		Petromyces albertensis	IBT 14317	14317	41046					Material shipped	Sequencing started						
2015-01-15 10:57	Mikael Rordam Andersen	allahabadii			CBS 164.63	23179	1131641					Material shipped	Sequencing started						
2015-01-15 10:57	Mikael Rordam Andersen	alliaceus		Petromyces alliaceus	CBS 536.65	13376	209559					Material shipped	Sequencing started						

........
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
								 usage='%(prog)s -dbname [database name] -source [Podio export Excel file]\n'
"# Example: python %(prog)s -dbname aspminedb -source /home/tcve/Downloads/Aspergillus_species.xlsx\n\
# Program reads a EXCEL file as exported from Podio, from the Aspergillus genome sequencing workspace.\n\
# Export Excel from the Aspergillus species app.\n\
# Created on	Created by	Species	......	Strain CBS Acc. No	Strain IBT No	NCBI Organism ID	MycoBank	JGI code	Picture	Section......\n\
# 2015-01-15 11:55	MR	acidus	Aspergillus foetidus		CBS 106.47				Aspfo1		Nigri...." )

#------------------------------------------------------------------
# Choose if data should be loaded or only tested
#------------------------------------------------------------------

parser.add_argument( "-source",		required=True,	help="R|File name, example: -source Aspergillus_species.xlsx" )
parser.add_argument( "-dbname",		required=False, help="R|Database name", default = "aspminedb" )

args = parser.parse_args()

#------------------------------------------------------------------
# Print argument values to screen
#------------------------------------------------------------------
print "#--------------------------------------------------------------"
print "# ARGUMENTS :\n# source\t:%s\n# database\t:%s" % ( args.source,args.dbname )
print "#--------------------------------------------------------------"

#------------------------------------------------------------------
# Connect to specific DB
#------------------------------------------------------------------
try:
    db = mdb.connect( "localhost","asp","1234",str( args.dbname ) )
    db.close() 
except mdb.Error, e:
    sys.exit( "# ERROR %d: %s" % ( e.args[0],e.args[1] ) )

db = mdb.connect( "localhost","asp","1234",str( args.dbname ) )
cursor = db.cursor()

#------------------------------------------------------------------
# Get lines from file
#------------------------------------------------------------------
workbook = xlrd.open_workbook(args.source)
worksheets	= workbook.sheet_names()
worksheet	= workbook.sheet_by_name(worksheets[0])
headers		= worksheet.row(0) # get first row as headers

for i in range(0,len(headers)): # loop through each element in header
	if re.match("^Section$",	headers[i].value):	section_index	= i
	if re.match("^Species$",	headers[i].value):	species_index	= i
	if re.match("^JGI code$",	headers[i].value):	jgi_index		= i 

if not section_index or not species_index or not jgi_index:
	sys.exit("# ERROR: one of the required fields were not found, Species, JGI code, Section")

updated = 0
for r in range(0, worksheet.nrows): # run through each row in sheet
	section_value	= worksheet.cell_value(r, section_index)
	species_value	= worksheet.cell_value(r, species_index)
	jgi_value 		= worksheet.cell_value(r, jgi_index)

	if not section_value or not species_value or not jgi_value:
		print ("# ERROR: one of the required value were not found, Species %s, JGI code %s, Section %s" % (species_value, section_value, jgi_value ))
	
	else:	
		try:	
			cursor.execute( "UPDATE organism SET real_name = '%s', section = '%s' WHERE name = '%s';" % (species_value, section_value, jgi_value ))
			db.commit()	# Add changes to database
			updated += 1
			
		except mdb.Error, e:
			sys.exit( "# ERROR organism name %s %d: %s" % ( args.source, e.args[0],e.args[1] ) )

print "# FINISHED: updated %s organism rows" % updated

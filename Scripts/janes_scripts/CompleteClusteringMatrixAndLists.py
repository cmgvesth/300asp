
import os, datetime, getpass
from os.path import isfile, join
import sys
from sys import argv
import argparse
from argparse import ArgumentParser
import re
import pandas as pd
from pandas import Series, DataFrame
import operator
import errno
import csv
import MySQLdb as mdb


#===================================================================
# INPUT ARGUMENTS DESCRIPTIONS
#===================================================================
def ParseArguments():
	""" Handles program parameters and returns an argument class containing 
	all parameters """
	parser = argparse.ArgumentParser(description="Preprocessing databases for gene clustering", 
									usage="%(prog)s -out filename")
	parser.add_argument("--outfile", "-o", type = str, required=False, 
									default = "initial.csv", help="Name of output file")
	parser.add_argument("--user", "-u", type = str, choices=["setd", "jlnr"], default = "jlnr", dest="user",
									required=True, help="Specify user")
	parser.add_argument("--database", "-d", type = str, required=False, default = "testasp", help="Database in MySQL",
									dest="database", choices=["testasp", "aspminedb", "joomla", "phpmyadmin", "asp"])
	parser.add_argument("--coverage", "-c", type = int, required=False, default = 90, dest= "coverage",
									help="Minimum alignment coverage")
	parser.add_argument("--identity", "-i", type = int, required=False, default = 50, dest= "identity",
									help="Minimum alignment identity")
	parser.add_argument("--species", "-s", type = str, required=False, default = "Aspnid1", dest= "species",
									help="Name of the species genome")

	args = parser.parse_args()

	# PRINT RUNNING DESCRIPTION 
	now = datetime.datetime.now()
	print '# ' + ' '.join(argv)
	print '# ' + now.strftime("%a %b %d %Y %H:%M")
	print '# USER: ' + getpass.getuser()
	print '# CWD: ' + os.getcwd()
	if os.name != "nt":
		print '# HOST: ' + ' '.join([ os.uname()[0] , os.uname()[4] , 
	       	                         os.uname()[2] , os.uname()[1] ])							
	return args


#===================================================================
# Create hit matrix from database
#===================================================================

def make_table(args):
	""" Checking program parameters and printing the arguments for 
	visual validation """
	data = []
	index=[]
	

	# Connect to specific database
	try:
		db = mdb.connect(host="192.38.13.9", user=args.user ,passwd="1234",db=args.database)
	except mdb.Error, e:
		sys.exit("Cannot connect to database")#"# ERROR %d: %s" % (e.args[0],e.args[1]))
	try:
		cursor = db.cursor()
	except mdb.Error, e:
		sys.exit("# ERROR %d: %s" % (e.args[0],e.args[1]))


	# Extraxt index from CDS table
	try:
		# Build index_query
		index_query = """SELECT `q_tail` FROM `blast` WHERE (`h_org` LIKE %s AND 
					`q_org` LIKE %s) AND (`q_tail` = `h_tail`) ORDER BY `q_tail`;""" 
		cursor.execute(index_query, (str(args.species),str(args.species)))
		indexes = cursor.fetchall()
	except mdb.Error, e:
		sys.exit('Setting up initial table failed')


	# EXTRACT DATA FROM DATABASE
	# When creating data_query. Generates two columns containing shortest alignment and 
	# longest alignment pct. coverages ((seq. end- seq. start)/(seq. length)).
	# The alignment length is chosen as the shortest alignment length, because it resembles 
	# the amount of aligned nucleotides the most (Adjusting for BLAST GABS - shortes length = fewer gabs).
	# Included selfhits to be able to cluster on patterns (therefore are selfhits needed).
	# Retriveing only hits with a coverage (either in short or long coverage) >= "predecided number"
	try:
		# Build data_query
		data_query = """SELECT * FROM (SELECT `q_tail`, `h_tail`,
				ROUND(((((LEAST((`h_end`-`h_start`), (`q_end`-`q_start`)))+1)/GREATEST(`q_len`,`h_len`))*100),2) AS 'pct_longCov',
				ROUND(((((LEAST((`h_end`-`h_start`), (`q_end`-`q_start`)))+1)/LEAST(`q_len`,`h_len`))*100),2) AS 'pct_shortCov',
				`pident` FROM `blast`
				WHERE `h_org` LIKE %s AND `q_org` LIKE %s) AS C
				WHERE (GREATEST(C.`pct_longCov`, C.`pct_shortCov`) >= %s AND C.`pident` >= %s);"""

		#Execute the data_query and save in variable
		cursor.execute(data_query, (str(args.species), str(args.species), int(args.coverage), int(args.identity)))
		data = cursor.fetchall()
	except mdb.Error, e:
		sys.exit('Setting up initial table failed')

	db.close()


	# Create table
	# BUG in .loc (the iteration app for name-location instead of integer-location) so the iteration is very slow.
	# Fast fix - changed to .iloc and found the indexed location (which is fixed in the matrix) in the indexing (protein_names) list
	protein_names = []
	for item in indexes:
		protein_names.append(item[0])

	hit_matrix = pd.DataFrame(index=protein_names, columns=protein_names)
	# Add 1's to the right cells
	for row in data:
		hit_matrix.iloc[protein_names.index(row[0]),protein_names.index(row[1])] = 1
	# Check for the right dimension
	print "\nMatrix dimensions: ",hit_matrix.shape

	# Write to csv file
	#hit_matrix.to_csv('hit_matrix.csv', sep=';',  encoding='utf-8') # Problems in outprint

	#Sort hit matrix
	#hit_matrix_sorted = hit_matrix.sort(axis=0, ascending=False)
	#hit_matrix_sorted.to_csv('hit_matrix_sorted.csv', sep=';',  encoding='utf-8') # Problems in outprint

	#hit_matrix['compressed'] = hit_matrix.apply(lambda x: ''.join([ str(v) for v in x ]),1)
	hit_matrix_filled = hit_matrix.fillna(0)
	
# 1. Create families version 1 (grouping on label hits - series of tuples)
	concat_matrix = pd.concat([ Series([ tuple(x[x.astype(bool)].index.tolist()) ], index=[row]) for (row,x) in hit_matrix_filled.iterrows() ])
	concat_matrix['compressed'] = concat_matrix.apply(lambda x: ''.join([ str(v) for v in x ]),1)
	protein_groups = concat_matrix.groupby('compressed').apply(lambda x: x.index.tolist())
	#print blabla[10:20]
	singleton_list = list()
	family_list = list()
	singleton_count = 0
	family_count = 0
	

	# Test if any proteins are located in multiple groups or are registred both in a group and as a selfhit
	# This test might be redundant.
	now = datetime.datetime.now()
	f_familylist = open(args.species+"_COV%s_ID%s_family_%s.txt"%(args.coverage,args.identity,now.strftime("%Y-%m-%d")), "w")
	f_familyshort = open(args.species+"_COV%s_ID%s_family_WO_singletons_%s.txt"%(args.coverage,args.identity,now.strftime("%Y-%m-%d")), "w")
	for element in protein_groups:
		if len(element) > 1:
			family_list.append(element)
			family_count += 1
			f_familylist.write('%s\t%s\n' %(element,len(element)))
			f_familyshort.write('%s\t%s\n' %(element,len(element)))
			#print element
		else:
			singleton_list.append(element)
			f_familylist.write('%s\t%s\n' %(element,len(element)))
			singleton_count += 1
	f_familylist.close()
	f_familyshort.close()
	print "# of singletons: ", singleton_count
	print "# of families: ", family_count
	print "Total families and singletons formed: ", len(protein_groups)
	# Check for location in both a group and in selfhits
	for singleton_element in singleton_list:
		for family in family_list:
			for family_member in family:
				if singleton_element == family_member:
					print "There are proteins in multiple locations:\n%s\t%s" % (singleton_element, family_list)
	# Check for location in in the different groups
	for n in range(len(family_list)-1):
		for family_member in family_list[n]:
			for check_member in family_list[n+1]:
				if family_member == check_member:
					print "There are proteins in multiple locations:\n%s\t%s" % (family[n], family_list[n+1])

	

# 2. Create families version 2 (grouping on label hits - dataframe)
	# test_matrix = DataFrame(dict([ (row,x[x.astype(bool)].index.tolist()) for (row,x) in hit_matrix_filled.iterrows() ])).T
	# test_matrix['compressed'] = test_matrix.apply(lambda x: ''.join([ str(v) for v in x ]),1)
	# blabla = test_matrix.groupby('compressed').apply(lambda x: x.index.tolist())
	# print blabla[0:6]
	# for element in blabla:
	# 	print element

	
# 3. Create families version 3 (grouping on 1 and 0 pattern - slow)
	# hit_matrix_filled['compressed'] = hit_matrix_filled.apply(lambda x: ''.join([ str(v) for v in x ]),1)
	# blabla = hit_matrix_filled.groupby('compressed').apply(lambda x: x.index.tolist())
	# print len(blabla)
	# print blabla[0:6]
	# for element in blabla:
	# 	print element
	# real	6m17.373s
	# user	6m5.319s
	# sys	0m5.036s




#===================================================================
# MAIN PROGRAM
#===================================================================
def main(argv):
	args = ParseArguments()
	make_table(args)			



# Running the main program
if __name__ == '__main__':
	main(argv)
	print "\n### Your program has finished! ###\n"



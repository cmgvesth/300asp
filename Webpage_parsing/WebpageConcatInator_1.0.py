#!/usr/bin/python

""" DESCRIPTION:
	The program takes a predefined list of webpages (with no pictures) and "concatenates" them into a single HTML-file.
	0. Input arguments and help program
	1. Printing overview of the arguments onto screen
	2. Opening and concatenating the webpages
 
INPUT FILE FORMATS:
	Input file: .txt file with list of webpages
  
Help:
	python ./WebpageConcatInator_1.0.py -h
  
Run example:
	python ./WebpageConcatInator_1.0.py -i ./WebpageListTestSet.txt -o Out.html
	
"""
import os, datetime, getpass
from os.path import isfile, join
import sys
from sys import argv
from argparse import ArgumentParser
from bs4 import BeautifulSoup # HTML handler
import requests # HTTP parser
import csv
# import pandas as pd
# from pandas import Series, DataFrame


# INPUT ARGUMENTS DESCRIPTION
# 0. Saving input arguments + help program
def ParseArguments():
	""" Handles program parameters and returns an argument class containing 
	all parameters """
	parser = ArgumentParser(description='Description of program: The program takes a predefined list of webpages (with no pictures) and "concatenates" them into a single HTML-file.')
	# Input parameters
	parser.add_argument("-i", type=str, dest="input_file", 
	default="", help="List of webpages. Required input", required=True)
	parser.add_argument("-o", type=str, dest="output_file", 
	default="OutputFile.html", help="Name of the output file. Required input", required=True)
						
	args = parser.parse_args()
	
	# PRINT RUNNING DESCRIPTION 
	now = datetime.datetime.now()
	print '# ' + ' '.join(argv)
	print '# ' + now.strftime("%a %b %d %Y %H:%M")							
	return args

# 1. Checking input
def InputChecker(args):
	if os.path.os.path.isfile(args.input_file) == False:
		return sys.exit("\nThe input file do not exist: %s\nRerun the program with a new file entry.\n" %args.input_file)
	elif os.path.os.path.isfile(args.output_file) == True:
		return sys.exit("\nThe output file already exists: %s\nRerun the program with a new output file.\n" %args.output_file)
	else:
		print '\n# COMMAND LINE PARAMETERS SET TO:'
		print "#\t Input file including path: ", args.input_file
		print "#\t Output file  including path:  ", args.output_file
		print ""









# Running the program
# 0. Importing arguments
args = ParseArguments()
# 1. Checking arguments
InputChecker(args)

# Setting up for a single abstract:
fh_website_list = open(args.input_file, "r")
website_list = list(csv.reader(fh_website_list, delimiter='\t'))
fh_website_list.close()

fh_output = open(args.output_file, 'w')

#Set flag for first run
flag = 1


for website in website_list:
	print website
	# This bit uses requests to get the HTML-part of the websites
	abstract_object = requests.get(''.join(website))
	HTML = abstract_object.text
	# Here Beautiful Soup takes over handling the HTML
	soup = BeautifulSoup(HTML)
	if flag == 1:
		fh_output.write('<html>')
		fh_output.write(soup.head.prettify("latin-1"))
		fh_output.write(soup.style.prettify("latin-1"))
	fh_output.write(soup.body.prettify("latin-1"))
	# After first run flag is set to 0
	flag = 0
fh_output.write('</html>')



fh_output.close()



'''
# Defining set of organisms (Test phase only: To be substituted with calling of arguments.)
JGI_ID_Tuple = ("Aspni7","Asplac1")
organism_Dict = dict()


# Connecting to database
db = MySQLdb.connect(host="192.38.13.9",user="mr",passwd="1234",db="aspminedb" )
cursor = db.cursor()

# Build a JGI_name-to-org_id dictionary
for id in JGI_ID_Tuple:
	# Building Query
	query = "SELECT org_id, name FROM organism WHERE name = %s"
	#Execute the query to identify lines in the table of the correct organism ID
	cursor.execute(query, id)
	data = cursor.fetchone()
	organism_Dict[data[0]] = data[1]

# To build list of unique clusters, cycle through data frame (sigh...)

ClusterGeneList = list()

for org_id in organism_Dict.keys():
	
	# Building Query
	query = "SELECT org_id, clust_backbone, clust_size, sm_protein_id FROM antismash WHERE org_id = %s"
	#Execute the query to identify lines in the table of the correct organism ID
	cursor.execute(query, org_id)
	data = cursor.fetchall()
	for row in data:
		ClusterGeneList.append( ( str(row[0]) + "_" + str(row[1]) + "_" + str(row[2]) , str(row[3]) ) ) 


# Now ClusterGeneList holds all protein IDs from the organisms in organism_tuple (index 1), and a unique Species_BackboneProtein_ClusterSize Identifier

UniqueClusters = set()
for gene in ClusterGeneList:
	UniqueClusters.add(gene[0])

UniqueClusters = list(UniqueClusters)
UniqueClusters.sort()

Results_df = pd.DataFrame(index=UniqueClusters, columns=UniqueClusters)


# Start building the results:

for Cluster in UniqueClusters:
	ClusterSet = list()
	Organism = int(Cluster.split('_')[0])
	for line in ClusterGeneList:
		if line[0] == Cluster:
			ClusterSet.append(list(line))
	# Now ClusterSet contains all members of the cluster with the geneID present in column index 1 in each item of the list
	# Examine this cluster in all genomes
	for org_id2 in organism_Dict.keys():
		# 1. look up gene in org_id1
		# Generate placeholder for gene-links
		ExpandedClusterSet = list()

		for n in range(len(ClusterSet)):
			query = "SELECT q_seqkey, h_seqkey, blast_bitscore FROM blast WHERE q_org = %s and h_org = %s and q_seqkey = %s"
			cursor.execute(query, (organism_Dict[Organism], organism_Dict[org_id2], ClusterSet[n][1] ) )
			# 1b. pick best hit
			# To find the best hit, loop through results and save the single best hit in BestHit_List
			hit = cursor.fetchone()
			# Check for absence of hit and assign None if none found:
			if hit is None:
				BestHit_List = [ClusterSet[n][1], 'None', 'None']
			else:
				BestHit_List = list(hit)
			
			while hit is not None: # Clunky loop, fix later!
				if hit[2] > BestHit_List[2]:
					BestHit_List = list(hit)
				hit = cursor.fetchone()
			# Now gene is identified and index 2 (bit score) is not necessary any longer and can be replaced by name of cluster
			
			# 2. look up geneID in antismash considering org_id2
			if BestHit_List[1] is not 'None':
				query = "SELECT org_id, clust_backbone, clust_size, sm_protein_id FROM antismash WHERE org_id = %s and sm_protein_id = %s"
				cursor.execute(query, (org_id2, BestHit_List[1]))
				org2_cluster_hit = cursor.fetchone()
				if org2_cluster_hit is not None: 
					BestHit_List[2] = str(org2_cluster_hit[0]) + "_" + str(org2_cluster_hit[1]) + "_" + str(org2_cluster_hit[2])
				else:
					BestHit_List[2] = 'None'
			else:
				BestHit_List[2] = 'None'
			# 3. build new table with orthologous gene cluster identifiers
			ExpandedClusterSet.append([ClusterSet[n][0], BestHit_List[0], BestHit_List[1], BestHit_List[2]])
			
						
		# 4. count numbers of genes in each cluster
		ClusterCounting_dict = dict()
		ClusterNames_Set = set()
		0
			if ExpandedClusterSet[n][3] is not 'None':
				ClusterNames_Set.add(ExpandedClusterSet[n][3])
		# print ClusterNames_Set
		for key in ClusterNames_Set:
			ClusterCounting_dict[key] = 0
		# Untested from here on, pending rebuild of BLAST table
		for n in range(len(ExpandedClusterSet)):
			if ExpandedClusterSet[n][3] is not 'None':
				ClusterCounting_dict[ExpandedClusterSet[n][3]] = ClusterCounting_dict[ExpandedClusterSet[n][3]] + 1
		# Remember to handle empty sets
		for keys in ClusterCounting_dict:
			# 5. Write results to dataframe based on Query_Cluster versus Hit_Cluster
			Results_df.loc[str(ExpandedClusterSet[0][0]),str(keys)] =  ClusterCounting_dict[keys]
		
db.close()

Results_df.to_csv('Results.csv')

'''


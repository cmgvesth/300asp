import sys,os
sys.path.append("../utils/")
from aspmine_imports import *
from seth_imports import *
import pandas as pd

parser = argparse.ArgumentParser(description="Preprocessing databases for gene clustering", usage="%(prog)s --out filename")
parser.add_argument("--met_classes", "-mc", required=False, action='store_true', help="Get a list of compound classes")
parser.add_argument("--met_compounds", "-ms", required=False, action='store_true', help="get a list of single compounds")
parser.add_argument("--get_list", "-gl", required=False, action='store_true', help="Print a list of class/uniques")

parser.add_argument("--orgs", "-og", required=False, action='store_true', help="Screen for compounds in Organisms")
parser.add_argument("--seq", "-s", required=False, action='store_true', help="Create genebank files from our data")

args = parser.parse_args()
met_classes = args.met_classes
met_compounds = args.met_compounds
get_list = args.get_list
orgs = args.orgs
seq = args.seq
cursor = asp_con('192.38.13.9', 'setd', '1234')



if met_classes:
	search_term = "metabolite1"

if met_compounds:
	search_term = "metabolite2"

if get_list:
	query ="SELECT DISTINCT(%s) FROM metaboliteMap" % search_term
	# print query
	cursor.execute(query)
	handle = cursor.fetchall()

	query ="SELECT %s,count(*) from metaboliteMap where status = '1' group by metabolite1;" % search_term
	# print query
	cursor.execute(query)
	overview = cursor.fetchall()
	print "Summary of metabolites found in organisms:\n# Metabolite: Number of Orgs"
	counter = 0
	for i in overview:
		print "%s %s : %s" % (counter, i[0], i[1])
		counter += 1


if orgs:
	print "Creating lists of organisms..."
	metabolite_org = {}
	m_org_counts = {}
	for i in handle:
		sub_handle = i[0]
		query = "SELECT metabolite1, real_name FROM metaboliteMap where metabolite1 = '%s' and status ='1';" % sub_handle
		# print query
		cursor.execute(query)
		handle = cursor.fetchall()

		for line in handle:
			a,b = line
			if a in metabolite_org:
				metabolite_org[a].append(b)
			else:
				metabolite_org[a] = list()
				metabolite_org[a].append(b)
	with open('met_in_orgs.csv', 'w') as f:  # Just use 'w' mode in 3.x
	    w = csv.writer(f)
	    w.writerow(metabolite_org.items())
	org_index = raw_input("Please choose a class of metabolits you want to work on: ")
	
	print "Processing %s" % org_index
	sub_query = "('"+"','".join(metabolite_org[org_index])+"','.*')" # TODO implement option to choose from metabolite
	minimum_orgs = len(metabolite_org[org_index])
	size = '2'
	query = "SELECT q_clustid, q_clust_size, candidateMembers ,orgs_repr from (\
		SELECT *, count(*) as orgs_repr from t_antismashLoopAntismashCandidates\
		where h_realname IN %s and clustCov >= 0.5 group by q_clustid) ta\
		where orgs_repr >= %s and q_clust_size > %s;" % (sub_query, minimum_orgs, size)
	cursor.execute(query)
	testing = cursor.fetchall()
	print testing


if seq:
	pass






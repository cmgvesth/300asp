import sys,os
sys.path.append("../utils/")
from aspmine_imports import *
import pandas as pd

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

def fetch_seq(infile,outfile):
	with open('best.csv') as f:
	reader = csv.reader(f)

	for i in bestHits:
		clust, score = i
		print(clust)
		query = "select * from t_antismash2blast_reduced where clust_id = '%s' and h_org = 'Aspeuc1'; where ;" % (clust)
		print(query)
		try:
			cursor.execute(query)
			result = cursor.fetchall()
		except:
			print "Error in Analysis"

		genome = {}

		for line in result:
			a,b = line
			if a in genome:
				genome[a].append(int(b))
			else:
				genome[a] = list()
				genome[a].append(int(b))

	#f = open(outfile,'wb')
	#writer = csv.writer(f, dialect = 'excel')
	#writer.writerows(result)
	#f.close()


# Now writing objects
	for i in proteins:
		record = SeqRecord(Seq("MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF",
		                   IUPAC.protein),
		                   id="YP_025292.1", name="HokC",
		                   description="toxic membrane protein, small")
		print record
from Bio import SeqIO
from reportlab.lib.colors import red, grey, orange, green, brown, blue, lightblue, purple, fade, PCMYKColor
from Bio.Graphics import GenomeDiagram
import sys,os
sys.path.append("../utils/")
from aspmine_imports import *
from seth_imports import *
import pandas as pd

cursor = asp_con('192.38.13.9', 'setd', '1234')

recA = SeqIO.read("GC_11_382281_6.gb", "gb")
recB = SeqIO.read("GC_28_10278_8.gb", "gb")
recC = SeqIO.read("GC_28_3654_7.gb", "gb")

blue  = PCMYKColor(91.0,  43.0,  0.0, 0.0)
black = PCMYKColor(0, 0, 0, 100)

pal = fade(blue,[100,60,40,20,5,0]) # TODO better coloring!

records = [recA, recB, recC]

hits = []
tagcol = []
#all_h_tags=[]
rank_col = []
for i in records[1:]: # Leaving query cluster out # Exchange hardcoded query cluster later!
	tags = []
	for n in i.features:
		tags.append(n.qualifiers['locus_tag'][0][2:])
	tags_formatted = "("+str(tags)[1:-1]+")"
	print tags_formatted
	#all_h_tags.append(str(tags)[1:-1])
	query ="SELECT  h_seqkey,CAST(pident AS UNSIGNED) from t_antismash2blast_reduced where clust_id = '11_382281_6' and h_seqkey IN %s;" % tags_formatted
	cursor.execute(query)

	result = list(cursor.fetchall())
	query =" SELECT ta.q_seqkey, tx.h_seqkey  from t_antismash2blast_reduced as ta  left join (SELECT * from t_antismash2blast_reduced as tb where tb.clust_id = '11_382281_6' and tb.h_seqkey IN %s ) tx on ta.q_seqkey = tx.q_seqkey where ta.clust_id = '11_382281_6' and ta.h_seqkey = ta.q_seqkey; "% tags_formatted
	cursor.execute(query)
	ranks_raw = cursor.fetchall()
	#print ranks_raw


	ranks = {}

	rank_counter = 1

	for line in ranks_raw:
		a,b = line
		if b in ranks:
			if b!='NULL':
				ranks[b].append(rank_counter)
				rank_counter+=1
		else:
			if b!='NULL':
				ranks[b] = list()
				ranks[b].append(rank_counter)
				rank_counter+=1

	rank_col.append(ranks)



	tag_pid = {}

	for line in result:
		a,b = line
		if a in tag_pid:
			tag_pid[a].append(int(b))
		else:
			tag_pid[a] = list()
			tag_pid[a].append(int(b))
	hits.append(tag_pid)
	tagcol.append(tags)

#f_all_h_tags= ','.join(all_h_tags)

#formatted = "('"+str(f_all_h_tags)[1:-1]+"')"

print tagcol
print '\n'
print hits

rankA=range(7)[1:]
rankB=[]
rankC=[]

ranklist = [rankB, rankC]

colB=[]
colC=[]
cols=[colB,colC]

counter=-1
for i in tagcol:
	counter+=1
	for n in i:
		if n in rank_col[counter]:
			ranklist[counter].append(rank_col[counter][n][0])
		else:
			ranklist[counter].append(0)
		if n in hits[counter]:
			calc = round(hits[counter][n][0])
			if calc > 80:
				index=0
			elif calc >60:
				index = 1
			elif calc > 40:
				index = 2
			elif calc > 20:
				index = 3
			else:
				index = 4
			cols[counter].append(pal[index])
		else:
			cols[counter].append(black)



print ranklist

colA=[blue]*6

def set_color():
	pass


name = "Gene Cluster comparison"
gd_diagram = GenomeDiagram.Diagram(name)
max_len = 0

for record, gene_colors, rank in zip([recA, recB, recC], [colA, colB, colC], [rankA]+ranklist):
	max_len = max(max_len, len(record))
	gd_track_for_features = gd_diagram.new_track(1,
							name=record.name,
							greytrack=True,
							start=0, end=len(record))
	gd_feature_set = gd_track_for_features.new_set()

	i = 0
	print len(rank)
	print len(record.features)
	for feature, single_rank in zip(record.features, rank):
		#if feature.type != "gene": # Doesn't work with my records because they don't have genes, only cds.... should work but doesn't....
			#Exclude this feature
		#    continue
		gd_feature_set.add_feature(feature, sigil="BIGARROW",
								   color=gene_colors[i], label=True,
								   name = str(single_rank), # this is too much, maybe implement if other label is found: str(feature.qualifiers['locus_tag'])+'M'+str(single_rank),
								   label_position="start",
								   label_size = 10, label_angle=0)
		i+=1

gd_diagram.draw(format="linear", pagesize='A4', fragments=1,
				start=0, end=max_len)
gd_diagram.write(name + ".pdf", "PDF")
gd_diagram.write(name + ".eps", "EPS")
gd_diagram.write(name + ".svg", "SVG")

from Bio import SeqIO
from reportlab.lib.colors import grey, orange, green, brown, blue, lightblue, purple, fade, PCMYKColor, Color
from Bio.Graphics import GenomeDiagram
import sys,os
sys.path.append("../utils/")
from aspmine_imports import *
from seth_imports import *
import pandas as pd


cursor = asp_con(path='192.38.13.196', user='setd', pw='1234')

#recA = SeqIO.read("GC_11_382281_6.gb", "gb")
#recB = SeqIO.read("GC_28_10278_8.gb", "gb")
#recC = SeqIO.read("GC_28_3654_7.gb", "gb")

red = PCMYKColor(0,  100.0,  0.0, 0.0)
blue  = PCMYKColor(91.0,  43.0,  0.0, 0.0)
black = PCMYKColor(0, 0, 20.0, 0)

pal = fade(blue,[100,60,40,20,5,0]) # TODO better coloring!

#records = [recA, recB, recC] # for testing
records = seq_builder(q_cluster = '11_382281_6', save = False)

# records = [recA, recB, recC]
##########################################
# Doing some black magic on rank numbers #
##########################################

hits = []
tagcol = []
#all_h_tags=[]
rank_col = []
for i in records[1:]: # Leaving query cluster out # Exchange hardcoded query cluster later!
	tags = []
	for n in i.features:
		if n.qualifiers['locus_tag'][2:] is None:
			print "cannot get rank for this protein"
			print n 
		else:
			tags.append(n.qualifiers['locus_tag'][2:]) # Hod to put a 0 here before... don't know why it doesn't need that now...
	tags_formatted = "("+str(tags)[1:-1]+")"
	#all_h_tags.append(str(tags)[1:-1])
	query ="SELECT  h_seqkey,CAST(pident AS UNSIGNED) from t_antismash2blast_reduced where clust_id = '11_382281_6' and h_seqkey IN %s;" % tags_formatted
	cursor.execute(query)

	result = list(cursor.fetchall())
	query =" SELECT ta.q_seqkey, tx.h_seqkey  from t_antismash2blast_reduced as ta  left join (SELECT * from t_antismash2blast_reduced as tb where tb.clust_id = '11_382281_6' and tb.h_seqkey IN %s ) tx on ta.q_seqkey = tx.q_seqkey where ta.clust_id = '11_382281_6' and ta.h_seqkey = ta.q_seqkey; "% tags_formatted
	cursor.execute(query)
	ranks_raw = cursor.fetchall()

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

#print rank_col
#f_all_h_tags= ','.join(all_h_tags)

#formatted = "('"+str(f_all_h_tags)[1:-1]+"')"

#print tagcol
#print '\n'
#print hits

rankA=range(7)[1:]

# Let's see if [[] for i in range(len(records))] works because some records from org 3 show empty features

ranklist = [[] for i in range(len(records))]
#print range(len(records))

cols=[[] for i in range(len(records))]

counter= -1
for i in tagcol:
	counter+=1
	for n in i:
		#print counter
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



#print ranklist

colA=[blue]*6

def set_color():
	pass


name = "testing"
gd_diagram = GenomeDiagram.Diagram(name)
max_len = 0

#print type([rankA]+ranklist)
#print [rankA]+ranklist

for record, gene_colors, rank in zip(records, [colA]+cols, [rankA]+ranklist):
	max_len = max(max_len, len(record))
	gd_track_for_features = gd_diagram.new_track(1,
							name=record.description + record.id,
							greytrack=True,
							greytrack_labels = 1,
							greytrack_font_rotation = 0,
							greytrack_font_colour = Color(1,0,0),
							greytrack_fontsize = 5,
							axis_labels = True,
							scale_smallticks = 0.6,
							start=0, end=len(record))
	gd_feature_set = gd_track_for_features.new_set()

	i = 0
	#print len(rank)
	#print len(record.features)
	#print rank
	for feature, single_rank in zip(record.features, rank):
		#if feature.type != "gene": # Doesn't work with my records because they don't have genes, only cds.... should work but doesn't....
			#Exclude this feature
		#    continue
		if feature.strand == -1:
			temp_angle = 180
		else:
			temp_angle = 1
		if single_rank == 0:
			temp_name = ''
		else:
			temp_name = str(feature.qualifiers['locus_tag'])+'M'+str(single_rank)
		gd_feature_set.add_feature(feature, sigil="BIGARROW",
								   color=gene_colors[i], label=True,
								   name = temp_name, # this is too much, maybe implement if other label is found: str(feature.qualifiers['locus_tag'])+'M'+str(single_rank),
								   label_position="middle",
								   label_size = 6, label_angle= temp_angle)
		i+=1

gd_diagram.draw(format="linear", pagesize='A4', fragments=1, orientation = 'portrait', track_size = 0.6, xl = 0.15,
				start=0, end=max_len)
gd_diagram.write(name + ".pdf", "PDF")
#gd_diagram.write(name + ".eps", "EPS")
gd_diagram.write(name + ".svg", "SVG")

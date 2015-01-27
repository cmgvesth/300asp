# Loading packages and creating a dictionary out of the table. I will append values to a key if that key is already existing and create a new key, if not.

import csv
import numpy as np
second= dict()
with open('smallertable.tsv') as tsv:
    for line in csv.reader(tsv, dialect='excel-tab'):
        second.setdefault(line[0],[]).append(line[1])

del second['blast_qseq_jg3']

names = set(sum(second.values(), []))
netLen = len(names)
network = np.zeros(shape=(len(second)+1,netLen+1), dtype='object')
counter = 0


# This will assign the first matrix column all dictionary keys.
for i in second:
	counter+=1
	try:
		network[counter, 0] = i 
	except:
		print"This entry could not be appended to the array %s" % i

counter=0


# This will assign all values to the first row
for i in names:
	counter += 1
	try:
		network[0, counter] = i 
	except:
		print"This entry could not be appended to the array %s" % i

# Changing matrix value from 0 to 1 where there's a hit
counter = 0
for i in second.keys(): 
	counter+=1
	try:
		for n in second[i]:
			try:
				matchIndex = np.where(network[0] == n)[0][0]
			except:
				print'Cannot find gene hits in array'
			try:
				network[counter,matchIndex] = 1
			except:
				print'Cannot add network count'
	except:
		'Something is wrong'

# Here comes the tricky part. The script will access every element in the network and will compare its matches to 
# the matches of the other elements. If they have a common match, they will be fused (this needs some fixing, sometimes it works
# but there should be more cases)
counter = 0
for i in network[1:]:
	counter+=1
	
	for n in network[1:]:
		if n[0]!=i[0]:
			if sum(i[1:]*n[1:]) >= 1:
				i[1:] = i[1:]+n[1:] #-----------------------------------> This has to be changed because it changes the value of matches
				print'Found single linkage!'
				if (type(i[0]) == str and type(n[0]) == str):
					i[0] = [i[0],n[0]]
					matchIndex = np.where(network[:,0] == n[0])[0][0]
					network = np.delete(network, matchIndex, axis = 0)
				elif type(i[0]) == list:
					i[0].append[n[0]]
					matchIndex = np.where(network[:,0] == n[0])[0][0]
					network = np.delete(network, matchIndex, axis = 0)
				elif type(i[0]) == str and type(n[0]) == list:
					i[0]=[i[0]]
					i[0].extend[n[0]]
					matchIndex = np.where(network[:,0] == n[0])[0][0]
					network = np.delete(network, matchIndex, axis = 0)
				else:
					print"Cannot add id to id list!"


family=list()
family_count=list()

for i in network[1:,1:]:
	try:
		family_count.append(sum(i))
	except:
		print "Cannot calculate family count"


import matplotlib.pyplot as plt
from numpy.random import normal
plt.hist(family_count)
plt.title("NidulansVsNidulans Single linkage")
plt.xlabel("Number of genes in family")
plt.ylabel("Frequency")
# If we were to simply plot pts, we'd lose most of the interesting
# details due to the outliers. So let's 'break' or 'cut-out' the y-axis
# into two portions - use the top (ax) for the outliers, and the bottom
# (ax2) for the details of the majority of our data
f, (ax, ax2) = plt.subplots(2, 1, sharex=True)
# plot the same data on both axes
ax.hist(family_count)
ax2.hist(family_count)
# zoom-in / limit the view to different portions of the data
ax.set_ylim(120, 150) # outliers only
ax2.set_ylim(0, 20) # most of the data
# hide the spines between ax and ax2
ax.spines['bottom'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax.xaxis.tick_top()
ax.tick_params(labeltop='off') # don't put tick labels at the top
ax2.xaxis.tick_bottom()
# This looks pretty good, and was fairly painless, but you can get that
# cut-out diagonal lines look with just a bit more work. The important
# thing to know here is that in axes coordinates, which are always
# between 0-1, spine endpoints are at these locations (0,0), (0,1),
# (1,0), and (1,1). Thus, we just need to put the diagonals in the
# appropriate corners of each of our axes, and so long as we use the
# right transform and disable clipping.
d = .015 # how big to make the diagonal lines in axes coordinates
# arguments to pass plot, just so we don't keep repeating them
kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
ax.plot((-d, +d), (-d, +d), **kwargs) # top-left diagonal
ax.plot((1 - d, 1 + d), (-d, +d), **kwargs) # top-right diagonal
kwargs.update(transform=ax2.transAxes) # switch to the bottom axes
ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs) # bottom-left diagonal
ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs) # bottom-right diagonal
# What's cool about this is that now if we vary the distance between
# ax and ax2 via f.subplots_adjust(hspace=...) or plt.subplot_tool(),
# the diagonal lines will move accordingly, and stay right at the tips
# of the spines they are 'breaking'

plt.show()
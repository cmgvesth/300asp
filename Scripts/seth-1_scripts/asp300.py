##########################################################################
# Loading data into dictionary
##########################################################################

def loadNidulans(data = 'table.tsv'):
	import csv
	import numpy as np
	second= dict()
	
	# Return a dictionary made out of the table-----------------------> I will append values to a key if that key is already existing and create a new key, if not.
	
	with open(data) as tsv:
	    for line in csv.reader(tsv, dialect='excel-tab'):
	        second.setdefault(line[0],[]).append(line[1])

	del second['blast_qseq_jg3']

	return second


##########################################################################
# Create a network matrix with hits shown as 1; takes dictionary as input
##########################################################################

def netW(data):
	import numpy as np
	names = set(sum(data.values(), [])) # Returns unique list of all values inside the dictionary
	netLen = len(names)
	network = np.zeros(shape=(len(data)+1,netLen+1), dtype='object')
	
	counter = 0


# Assigning gene identifiers to matrix ---------------------------------> All keys assigned to first column and values assigned to first row
	for i in data:
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

	counter=0

	# Finding gene -------------------------------------------------------> matchindex function returns an array, I subset it to get the coordinates for the column where the hit occured
	for i in data.keys(): 
		counter+=1
		try:
			for n in data[i]:
				try:
					matchIndex = np.where(network[0] == n)[0][0]
				except:
					print'Cannot find gene hits in array'
				try:
					network[counter,matchIndex] = 1 	# Changing matrix value from 0 to 1 if there's a hit
				except:
					print'Cannot add network count'
		except:
			'Something is wrong'

###########################################################################
# Here comes the tricky part. The script will access every element in the network and will compare its matches to 
# the matches of the other elements. If they have a common match, they will be fused (this needs some fixing, sometimes it works
# but there should be more cases)
###########################################################################
	counter = 0
	for i in network[1:]:
		counter+=1
		
		for n in network[1:]:
			if n[0]!=i[0]:
				if sum(i[1:]*n[1:]) >= 1:
					i[1:] = i[1:]+n[1:] #-----------------------------------> This has to be changed because it changes the value of matches, hm... really?
					
					# Single linkage families are created by lists and this part just serves to create only one list and not list of lists

					if (type(i[0]) == str and type(n[0]) == str):
						i[0] = [i[0],n[0]]
						matchIndex = np.where(network[:,0] == n[0])[0][0]
						network = np.delete(network, matchIndex, axis = 0)

					elif type(i[0]) == list:
						i[0].append(n[0])
						matchIndex = np.where(network[:,0] == n[0])[0][0]
						network = np.delete(network, matchIndex, axis = 0)

					elif type(i[0]) == str and type(n[0]) == list:
						i[0]=[i[0]]
						i[0].extend(n[0])
						try:
							matchIndex = np.where(network[:,0] == n[0])[0][0]
						except:
							print "Cannot find entry to delete %s" % n
						try:
							network = np.delete(network, matchIndex, axis = 0)
						except:
							print "Cannot delete entry from network %s" % n

					else:
						print"Cannot add id to id list!"
	return network

################################################################################
# Plotting function with cut y-axis for better visibility
################################################################################

def plotNidulans(data):
	import numpy as np
	family_count=list()

	for i in data[1:,1:]:
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
	ax.set_ylim(90, 500) # outliers only
	ax2.set_ylim(0, 30) # most of the data
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
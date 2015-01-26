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
				i[1:] = i[1:]+n[1:]
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

#!/usr/bin/python

'''##################################################################
# Imports
##################################################################'''
import sys,os
sys.path.append(os.path.join(os.path.dirname(__file__), "/home/tcve/github/utils/"))
from aspmine_imports import *
import numpy as np
import scipy.cluster.hierarchy as hac
import matplotlib.pyplot as plt
import scipy.spatial.distance as ssd
from matplotlib import cm as cm

'''##################################################################
# ARGUMENTS and setup
##################################################################'''
startTime = datetime.now() # record runtime
parser = CustomArgumentParser(formatter_class=SmartFormatter, usage='%(prog)s -fasta [FASTA file name]')
parser.add_argument("-fasta", "-f", required=False, help="File containing multiple sequences, FASTA format")
parser.add_argument("-dist", "-d", required=False, help="File containing distance matrix from sequence alignment")

# clustalw -infile=all.aln -tree -pim

args 		= parser.parse_args()
fastafile 	= args.fasta
distfile 	= args.dist

if not distfile and not fastafile:
	sys.exit("# ERROR: neither FASTA or distance specified, please choose only one")

if distfile and fastafile:
	sys.exit("# ERROR: both FASTA and distance specified, please choose only one")

'''------------------------------------------------------------------
# Print argument values to screen
------------------------------------------------------------------'''
print "#--------------------------------------------------------------\n\
# ARGUMENTS\t\t:\n\
# FASTA query file\t: %s\n\
# Distance matrix file\t: %s\n\
#--------------------------------------------------------------" % (fastafile, distfile)

'''------------------------------------------------------------------
# Do alignment from FASTA file
# AND
# Get distance matrix
------------------------------------------------------------------'''
if fastafile and not distfile:
	source_exists("file", fastafile)
	#os.system("clustalw -infile=%s" % fastafile)
	os.system("clustalw -infile=%s -tree -pim" % fastafile.replace("fsa", "aln"))
	infile = open(fastafile.replace("fsa", "pim"), 'r')

elif distfile and not fastafile:
	source_exists("file", distfile)
	infile = open(distfile, 'r')


values = []
names = []

for line in infile.readlines():
	# skip comment and empty lines
	if re.search("#", line) or re.search("^$", line):
		continue
	line = re.sub('^\s+', "", line)
	line = re.sub('\s+', "\t", line)
	names.append( line.split("\t")[1] )
	values.append( [int(x) if x else 0 for x in line.split("\t")[2:]] )
	#distMatrix.append( [ int(x) if x else 0 for x in line.split("\t")[2:]] )
	#line = re.sub('', "0", line)
	#line = re.sub('\t\t', "\t0\t", line)

#print names

out = open(distfile.replace("pim", "tab"),'wb')

for x in range(0,len(names),1):
	for y in range(0,len(names),1):
		line = str(names[x]) + ";" + str(names[y]) + ";" + str(values[x][y]) 
		out.write(line + "\n")
out.close()		

#np_distMatrix = np.array(distMatrix)


"""
plt.pcolormesh(np_distMatrix)
plt.colorbar()
curves = 10
m = max([max(row) for row in np_distMatrix])
levels = np.arange(0, m, (1 / float(curves)) * m)
plt.contour(np_distMatrix, colors="white", levels=levels)
plt.show()
"""

#print np_distMatrix
#Y = hac.pdist(np_distMatrix, 'euclidean')
#print np_distMatrix

#sys.exit()
"""
fig, axes23 = plt.subplots(2, 4)

for method, axes in zip(['single', 'complete'], axes23):
	z = hac.linkage(np_distMatrix, method=method)

	# The knee estimates how many clusters to expect
	knee = np.diff(z[::-1, 2], 2)

	# Set number of clusters
	num_clust1 = knee.argmax() + 2
	knee[knee.argmax()] = 0
	num_clust2 = knee.argmax() + 2

	#( num_clust1,num_clust2 ) = (4,8)

	# Plotting
	axes[0].plot(range(1, len(z)+1), z[::-1, 2])
	axes[0].plot(range(2, len(z)), knee)

	# Add writing to plot
	axes[0].text(num_clust1, z[::-1, 2][num_clust1-1], 'possible\n<- knee point')

	part1 = hac.fcluster(z, num_clust1, 'maxclust')
	part2 = hac.fcluster(z, num_clust2, 'maxclust')

	dd = hac.dendrogram(hac.linkage(z))
	
	i = 0
	for p in part1:
		#print p, rowNames[i]
		i+=1

	clr = ['#2200CC' ,'#D9007E' ,'#FF6600' ,'#FFCC00' ,'#ACE600' ,'#0099CC' ,'#8900CC' ,'#FF0000' ,'#FF9900' ,'#FFFF00' ,'#00CC01' ,'#0055CC']
	
	for part, ax in zip([part1, part2], axes[1:]):
		for cluster in set(part):
			#print cluster, clr[cluster]
			#print np_distMatrix[part == cluster, 1], clr[cluster]
			ax.scatter(np_distMatrix[part == cluster, 0], np_distMatrix[part == cluster, 1], color=clr[cluster])
			
	m = '\n(method: {})'.format(method)
	plt.setp(axes[0], title='Screeplot{}'.format(m), xlabel='partition',ylabel='{}\ncluster distance'.format(m))
	plt.setp(axes[1], title='{} Clusters'.format(num_clust1))
	plt.setp(axes[2], title='{} Clusters'.format(num_clust2))
	
plt.tight_layout()
plt.show()
"""

'''------------------------------------------------------------------
# Write result to file
------------------------------------------------------------------'''
#out = open(fastafile+".fsa",'wb')
#out.write(">" + comment + "\n")
#out.write(seq + "\n")
#out.close()		

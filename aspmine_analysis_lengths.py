#!/usr/bin/python

#------------------------------------------------------------------
# Imports
#------------------------------------------------------------------
import MySQLdb as mdb
import sys, getopt, argparse, re, glob, os, gzip

from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

import numpy
import scipy
import scipy.stats
from datetime import datetime

sys.path.append('/home/tcve/github/utils/')
from utilsArgparse import * # custom functions

# API plotting
import plotly.plotly as py
from plotly.graph_objs import *

# Local plots
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

##############################################################################################################################################
##############################################################################################################################################
####################################################################### ARGUMENTS ############################################################
##############################################################################################################################################
##############################################################################################################################################
startTime = datetime.now() # record runtime

#------------------------------------------------------------------
# Get command line arguments
#------------------------------------------------------------------
parser = CustomArgumentParser(formatter_class=SmartFormatter, 
								description="Load primary metabolism model data into local mysql database", 
								 usage='%(prog)s -out filename\n'
								 "Example: python %(prog)s -dbname aspminedb -filetype model -action test -source NigerModel2013.csv")
parser.add_argument("-limit", "-l", required=False, default = 100, help="limit length selects, dafault 100, to select all set to 'all' ")
parser.add_argument("-out", "-o", required=False, default = "lengths.tab", help="Name of output file")
parser.add_argument("-dbname", "-d", required=False, default = "aspminedb", help="Database name")
parser.add_argument("-verbose", "-v", required=False, help="R|Increase output verbosity for bug tracking" , action='store_true')

args = parser.parse_args()
limit = args.limit
#------------------------------------------------------------------
# Print argument values to screen
#------------------------------------------------------------------
print "#--------------------------------------------------------------"
print "# ARGUMENTS :\n# Database\t:%s\n# Outfile\t:%s\n# Limit\t:%s" % (args.dbname, args.out, limit)
print "#--------------------------------------------------------------"

#------------------------------------------------------------------
# Connect to specific DB
#------------------------------------------------------------------
outfile = open(args.out, 'w') 

try:
    db = mdb.connect("localhost","asp","1234",str(args.dbname))

except mdb.Error, e:
    sys.exit("# ERROR %d: %s" % (e.args[0],e.args[1]))

#------------------------------------------------------------------
# Get protein lengths
#------------------------------------------------------------------
try:
	db = mdb.connect("localhost","asp","1234",str(args.dbname))
	cursor = db.cursor()
	if limit != "all":
		cursor.execute("SELECT org_id,prot_orgkey,prot_seqkey,prot_tailkey,prot_seqname,prot_orgseq,char_length(prot_seq) FROM `proteins` limit %s" % limit)
	else:
		cursor.execute("SELECT org_id,prot_orgkey,prot_seqkey,prot_tailkey,prot_seqname,prot_orgseq,char_length(prot_seq) FROM `proteins` ")
	
	#cursor.execute("SELECT char_length(prot_seq) FROM `proteins`")
	result = cursor.fetchall() 
except mdb.Error, e:
	sys.exit("# ERROR %d: %s" % (e.args[0],e.args[1]))


#------------------------------------------------------------------
# Display histogram of lengths
#------------------------------------------------------------------
lengths = []
for r in result:
	lengths.append(r[-1])

x = numpy.asarray( lengths )
#py.sign_in('Python-Demo-Account', 'gwt101uhh0')
#data = Data([ Histogram( x=x ) ])
#plot_url = py.plot(data, filename='basic-histogram')

print len(x), x.mean(), x.max()
print "data size\tmin value\tmax value\tmean\tunbiased variance\tvariance\tskewness\tkurtoises (0 if normal dist.)"
print scipy.stats.describe(x)
"""
(min, max): tuple of ndarrays or floats
minimum and maximum value of data array
arithmetic mean : ndarray or float
mean of data along axis
unbiased variance : ndarray or float
variance of the data along axis, denominator is number of observations minus one.
biased skewness : ndarray or float
skewness, based on moment calculations with denominator equal to the number of observations, i.e. no degrees of freedom correction
biased kurtosis : ndarray or float
kurtosis (Fisher), the kurtosis is normalized so that it is zero for the normal distribution. No degrees of freedom or bias correction is used.
"""

# The ``normed`` flag, which normalizes bin heights so that the integral of the histogram is 1. The resulting histogram is a probability density.
#plt.hist(x, 20, normed=1, histtype='bar', facecolor='g', alpha=0.75)
#plt.hist(x, 20, histtype='bar', facecolor='g', alpha=0.75, orientation="horizontal")
#plt.hist(x, 100, histtype='bar', facecolor='g', alpha=0.75)
#binwidth = 50
#plt.hist(x, bins=range(min(x), max(x) + binwidth, binwidth), histtype='bar', facecolor='g', alpha=0.75)
#plt.hist(x, bins=(range(0,2000, 100) + range(2000,max(x), 500)), histtype='bar', facecolor='g', alpha=0.75)
plt.style.use('bmh')

plt.hist(x, bins=(range(0,max(x),100)), histtype='stepfilled', facecolor='g', alpha=0.75)

plt.plot()
plt.xlabel('Protein lengths (amino acids)')
plt.ylabel('Number of proteins')
plt.title(r'Histogram of protein lengths')

#plt.subplots_adjust(left=0.15)
plt.show()

"""
# example data
mu = 100 # mean of distribution
sigma = 15 # standard deviation of distribution
#x = mu + sigma * numpy.random.randn(10000)

num_bins = 50
# the histogram of the data
n, bins, patches = plt.hist(x, num_bins, normed=1, facecolor='green', alpha=0.5)
# add a 'best fit' line
#y = mlab.normpdf(bins, mu, sigma)
plt.plot(bins, 'r--')
plt.xlabel('Smarts')
plt.ylabel('Probability')
plt.title(r'Histogram of IQ: $\mu=100$, $\sigma=15$')

# Tweak spacing to prevent clipping of ylabel
plt.subplots_adjust(left=0.15)
plt.show()
"""
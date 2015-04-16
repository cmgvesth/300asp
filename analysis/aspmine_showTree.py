#!/usr/bin/python

'''##################################################################
# Imports
##################################################################'''
import sys,os
root = os.path.abspath(__file__).split("github")[0]
sys.path.append(root+"github/utils/")
from aspmine_imports import *
from Bio.Align.Applications import ClustalwCommandline
from Bio import Phylo
import pylab
import matplotlib.pyplot as plt
import numpy as np

'''##################################################################
# ARGUMENTS and setup
##################################################################'''
startTime = datetime.now() # record runtime
parser = CustomArgumentParser(formatter_class=SmartFormatter, 
								description="Create distance tree from FASTA file", 
								 usage='%(prog)s -f [Tree file, newick] \n'
								 "Example: python %(prog)s -f hitseqs.fsa ")
parser.add_argument("-file", "-f", required=True, help="Name of input file")

args 	= parser.parse_args()
treefile = args.file
assert os.path.isfile(treefile), sys.exit("# ERROR File " + treefile + " does not exists")

'''------------------------------------------------------------------
# Print argument values to screen
------------------------------------------------------------------'''
print "#--------------------------------------------------------------\n\
# ARGUMENTS\n\
# Description\t: Create distance tree from FASTA file\n\
# Tree file\t: %s\n\
#--------------------------------------------------------------" % (treefile)

'''------------------------------------------------------------------
# Draw tree
------------------------------------------------------------------'''

# Open file and remove negatives
newlines = ""
phbfile = open(treefile, "r")

# Create tree from phb file
tree = Phylo.read(treefile, "newick", rooted=False)

# Pretty up branch names
"""
for clade in tree.get_terminals():
	cname = str(clade.name)
	if re.match("jgi|", cname):
		clade.name=(cname).replace("jgi|", "")

	if re.match("|", cname):
		clade.name=("_").join((cname).split("|")[0:2])
"""
# Draw ascii tree
fhascii = open(treefile + "_asciitree.txt", "wb")
Phylo.draw_ascii(tree, file=fhascii)
print "# INFO: created ascii tree in file %s" % treefile + "_asciitree.txt"

# Fix size of plot
nrclades = len(tree.get_terminals())
fig = plt.figure(figsize=(20,30+nrclades/10))
ax = fig.add_subplot(111)

# Plot tree
tree_draw = Phylo.draw(tree, do_show=False, branch_labels=(lambda c: c.branch_length > 0.02 and format( c.branch_length, '.2f') or None), axes=ax)
pylab.savefig(treefile + "_tree.pdf",  format="pdf", orientation="landscape", dpi=200)
print "# INFO: created PDF tree in file %s" % treefile + "_tree.pdf"
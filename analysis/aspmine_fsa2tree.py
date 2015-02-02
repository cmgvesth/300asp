#!/usr/bin/python

'''##################################################################
# Imports
##################################################################'''
import sys,os
sys.path.append(os.path.join(os.path.dirname(__file__), "../utils/"))
from aspmine_imports import *

'''##################################################################
# ARGUMENTS and setup
##################################################################'''
startTime = datetime.now() # record runtime
parser = CustomArgumentParser(formatter_class=SmartFormatter, 
								description="Create distance tree from FASTA file", 
								 usage='%(prog)s -f [FASTA file with nucleotides or proteins] -align [perform multiple alignment]\n'
								 "Example: python %(prog)s -f hitseqs.fsa -align")
parser.add_argument("-file", "-f", required=True, help="Name of input file")
parser.add_argument("-align", "-a", required=False, action='store_true', help="Perform multiple alignment")

args 	= parser.parse_args()
fastafile = args.file
align = args.align
assert os.path.isfile(fastafile), sys.exit("# ERROR File " + fastafile + " does not exists")

'''------------------------------------------------------------------
# Print argument values to screen
------------------------------------------------------------------'''
print "#--------------------------------------------------------------\n\
# ARGUMENTS\n\
# Description\t: Script return list of best bidirectional hit for a gene in a specific strain\n\
# FASTA file\t: %s\n\
#--------------------------------------------------------------" % (fastafile)

'''------------------------------------------------------------------
# Clean up fasta file
------------------------------------------------------------------'''




'''------------------------------------------------------------------
# Run alignment
------------------------------------------------------------------'''
from Bio.Align.Applications import ClustalwCommandline
from Bio import Phylo
import pylab
import matplotlib.pyplot as plt
import numpy as np

"""
Changed initial branch length value in class:
sudo gedit /usr/lib/python2.7/dist-packages/Bio/Phylo/Newick.py
"""

if not os.path.isfile(fastafile.replace("fsa", "aln")) or align:
	print "# INFO: alignment file not found, doing alignment"
	cline = ClustalwCommandline("clustalw", infile=fastafile, outfile=fastafile.replace("fsa", "aln"), type="protein", align=True)
	stdout, stderr = cline()
	cline = ClustalwCommandline("clustalw", infile=fastafile.replace("fsa", "aln"),  bootstrap=1000, output="PHYLIP")
	stdout, stderr = cline() # hitseqs.fsa -tree -type="protein" -align -bootstrap=1000

if not os.path.isfile(fastafile.replace("fsa", "aln")):
	sys.exit("# ERROR: alignment not found, something went wrong ")

#if not os.path.isfile(fastafile.replace("fsa", "ph")):
#	sys.exit("# ERROR: tree file ph not found, something went wrong ")

'''------------------------------------------------------------------
# Draw tree
------------------------------------------------------------------'''

tree = Phylo.read(fastafile.replace("fsa", "phb"), "newick")
#print tree

for clade in tree.get_terminals():
	if re.match("jgi|", clade.name):
		clade.name=(clade.name).replace("jgi|", "")

	if re.match("|", clade.name):
		clade.name=("_").join((clade.name).split("|")[0:2])

fhascii = open(fastafile.replace("fsa", "asciitree.txt"), "wb")
#Phylo.draw_ascii(tree)
Phylo.draw_ascii(tree, file=fhascii)
print "# INFO: created ascii tree in file %s" % fastafile.replace("fsa", "asciitree.txt")


tree_draw = Phylo.draw(tree, do_show=False, branch_labels=(lambda c: c.branch_length > 0.02 and format( c.branch_length, '.2f') or None))
pylab.savefig(fastafile.replace(".fsa", "_tree.pdf"),  format="pdf", orientation="landscape", dpi=72.72)

"""
x = np.linspace(0, 2 * np.pi, 100)
y = 2 * np.sin(x)
fig, (ax0, ax1) = plt.subplots(nrows=2)

ax0.plot(x, y)
ax0.set_title('normal spines')

ax1.plot(x, y)
ax1.set_title('bottom-left spines')

# Hide the right and top spines
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)
# Only show ticks on the left and bottom spines
ax1.yaxis.set_ticks_position('left')
ax1.xaxis.set_ticks_position('bottom')

# Tweak spacing between subplots to prevent labels from overlapping
plt.subplots_adjust(hspace=0.5)

plt.show()
"""
"""
ax = pylab.subplot(111)
#print type(ax)
ax.set_title('bottom-left spines')

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')

ax.plot(tree_draw)
pylab.show()
#print tree_draw
#axhspan=((0.25, 7.75), {'facecolor':'0.5'}))
#tree_draw.axis([0, 6, 0, 20])
pylab.savefig(fastafile.replace(".fsa", "_tree.pdf"),  format="pdf", orientation="landscape", dpi=72.72)
#print "# INFO: created pdf tree in file %s" % fastafile.replace("fsa", "tree.txt")

#Phylo.draw(tree, branch_labels=lambda c: str(c.branch_length)+"_"+str(c.comment))
"""
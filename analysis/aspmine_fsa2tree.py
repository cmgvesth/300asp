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
# Description\t: Create distance tree from FASTA file\n\
# FASTA file\t: %s\n\
#--------------------------------------------------------------" % (fastafile)

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

# CLEAN UP negative branch lengths
# Open file and remove negatives
newlines = ""
phbfile = open(fastafile.replace("fsa", "phb"), "r")

for line in phbfile.readlines():
	newlines += line.replace(":-0.", ":0.")
phbfile.close()

# Open file and write new content
phbfile = open(fastafile.replace("fsa", "phb"), "w")
phbfile.write(newlines)
phbfile.close()

# Create tree from phb file
tree = Phylo.read(fastafile.replace("fsa", "phb"), "newick", rooted=True)

# Pretty up branch names
for clade in tree.get_terminals():
	cname = str(clade.name)
	if re.match("jgi|", cname):
		clade.name=(cname).replace("jgi|", "")

	if re.match("|", cname):
		clade.name=("_").join((cname).split("|")[0:2])

# Draw ascii tree
fhascii = open(fastafile.replace(".fsa", "_asciitree.txt"), "wb")
Phylo.draw_ascii(tree, file=fhascii)
print "# INFO: created ascii tree in file %s" % fastafile.replace(".fsa", "_asciitree.txt")

# Fix size of plot
nrclades = len(tree.get_terminals())
fig = plt.figure(figsize=(20,30+nrclades/10))
ax = fig.add_subplot(111)

# Plot tree
tree_draw = Phylo.draw(tree, do_show=False, branch_labels=(lambda c: c.branch_length > 0.02 and format( c.branch_length, '.2f') or None), axes=ax)
pylab.savefig(fastafile.replace(".fsa", "_tree.pdf"),  format="pdf", orientation="landscape", dpi=200)
print "# INFO: created PDF tree in file %s" % fastafile.replace("fsa", "_tree.pdf")
#!/usr/bin/python2.7 -w

""" DESCRIPTION:
	Produce a matrix output of genes/proteins (=hits) in common between species.
	0.  Input argument and help program
	1.  Checking arguments and printing an overview of the input onto screen
	2.  Creating output directories
	3.  Generating lists containing hits shared between "speciesX vs speciesY" and "speciesY vs speciesX" 
			a. Checking BLAST input files
			b. Extracting top hits from input
			c. Generating lists containing shared hits by both species
	4. 	Producing error message if BLAST files are missing
	5.  Transforming the hit lists into matrices for each species
	6.  Generating a collected matrix containing hits patterns and the amounts of hits with these patterns 
			a. Assigning all hits a pattern
			b. Grouping  hits with same pattern
			c. Collecting patterns for each species into matrix, showing the sum of hits with these pattern
			d. Adding statistical calculations to the values in the matrix
	    Pattern: Hits shared between species = 1, not shared = 0
  
REMEMBER:
	To install pandas module before running this program.
	http://pandas.pydata.org/pandas-docs/stable/install.html
 
INPUT FILE FORMATS:
	Input file: tab separated .txt file. Preferably BLAST output file.
	The file NEEDS to be named: [speciesX]vs[speciesY]Table.txt
  
Help:
	python ./BLASTprofile.py -h
  
Run example:
	python ./scripts/BLASTprofile.py -i /Users/Jane/Dropbox/A.Niger/BLAST_out/ 
	-o /Users/Jane/Dropbox/A.Niger/
"""

#-------------------------------------------------------------------
# IMPORT MODULES
#-------------------------------------------------------------------

try:
	import os, datetime, getpass
	from os.path import isfile, join
	import sys
	from sys import argv
	from argparse import ArgumentParser
	import re
	import pandas as pd
	from pandas import Series, DataFrame
	import operator
	import errno
	import csv
except ImportError:
	print "One of your modules can't be imported"



#-------------------------------------------------------------------
# INPUT ARGUMENTS DESCRIPTIONS
#-------------------------------------------------------------------

def ParseArguments():
	""" Handles program parameters and returns an argument class containing 
	all parameters """
	parser = ArgumentParser(description='Description of program: The programs main function is to convert BLAST alignment files between two species <spXvsspY> and <spYvsspX> into a matrix containing the hits shared between these species. This is possible to do for multiple species at once. It also generates a collected matrix containing the amount of shared hits across all the input species.')
	parser.add_argument("-V", "--version", action="version", version="%(prog)s (version 0.1)")
	# Input parameters
	parser.add_argument("-i", type=str, dest="input_dir", 
						default="", help="directory of your input files. Required input.", required=True)
	parser.add_argument("-o", type=str, dest="output_dir", 
						default="", help="directory of your output files. Required input", required=True)
	parser.add_argument("-s", choices=['mean', "median", "sd", 'standard', 'pct', 'all', 'none'], type=str, dest="select",
					default="none", help="selection of which columns you would include in your sum matrix. 'mean' = mean, 'median' = median, 'sd' = standard deviation (SD), 'standard' = SD, mean and median, 'pct' = pct SD from median and pct RSD, 'all' = all of the above, 'none' = none. Default = none.")
						
	args = parser.parse_args()
	# PRINT RUNNING DESCRIPTION 
	now = datetime.datetime.now()
	print '# ' + ' '.join(argv)
	print '# ' + now.strftime("%a %b %d %Y %H:%M")
	print '# USER: ' + getpass.getuser()
	print '# CWD: ' + os.getcwd()
	if os.name != "nt":
		print '# HOST: ' + ' '.join([ os.uname()[0] , os.uname()[4] , 
	       	                         os.uname()[2] , os.uname()[1] ])							
	return args


#-------------------------------------------------------------------
# SUB PROGRAMS
#-------------------------------------------------------------------
#===================================================================
# 1. Checking and printing an overview of input arguments
#===================================================================

def my_function(args):
	# Checking input arguments
	if os.path.exists(args.input_dir) == False:
		return sys.exit("\nThe input directory do not exist: %s\nRerun the program with a new directory.\n" %args.input_dir)
	elif os.listdir(args.input_dir) == []:
		return sys.exit("\nThe input directory do not contain any files: %s\nRerun the program with a new directory.\n" %args.input_dir)
	elif os.path.exists(args.output_dir) == False:
		return sys.exit("\nThe output directory do not exist: %s\nRerun the program with a new directory.\n" %args.output_dir)
	else:
		# Printing arguments for visual validation
		print '\n# COMMAND LINE PARAMETERS SET TO:'
		print "#\t Directory containing input files:\t", args.input_dir
		print "#\t Directory for saving output files:\t", args.output_dir
		print "#\t Selected option for protein count matrix:\t", args.select
		print ""



#===================================================================
# 2. Creating output directories
#===================================================================

def mkdir(args):
	root=args.output_dir+'/BLAST_profile' ; hits=root+'/BLASThits' ; matrix=root+'/matrix'; hitcounts=root+'/counts_patterns'
	path = [root, hits, hits+'/BLASTinfo_tophits', hits+'/shared_hits',matrix,matrix+'/unsorted',matrix+'/sorted',
	hitcounts, hitcounts+'/hitPattern',hitcounts+'/hitcounts']
	for folder in path:
		try:
			os.makedirs(folder)
		except OSError as exc:
			# If the folders already exist - then pass
			if exc.errno == errno.EEXIST and os.path.isdir(folder):
				pass
			else: raise
	return args.output_dir + '/BLAST_profile'



#===================================================================
# 3. Extracting top hints for each BLASTfile spXvsspY
#    and generating a list containing hits that are only present in both species tophits: sp1vssp2 = sp2vssp1
#===================================================================

def shared(args, root): 
	""" Output:
		1. Files containing shared hits between species 
		2. File containing # of top hits for one species
		3. File containing # shared hits between species 
	"""
	## 1. Checking input files and renaming them
	specieslist = [] ; BLASTlist = [] ; errorlist = []
	for files in os.listdir(args.input_dir):
		# Extracting species names
		sp = re.split('^(.+)[Vv][Ss](.+)[Tt].*.txt$', files)
		# Renaming files
		# And adjusting for human capitalization errors
		try:
			os.rename(join(args.input_dir,files), join(args.input_dir,"%sVs%sTable.txt" %(sp[1].capitalize(), sp[2].capitalize())))
		except IndexError:
			if files != ".DS_Store":
				errorlist.append(files)
			continue
		# Checking that both files are present sp1vssp2 and sp2vssp1
		if os.path.isfile(join(args.input_dir,"%sVs%sTable.txt" %(sp[1], sp[2]))) and os.path.isfile(join(args.input_dir,"%sVs%sTable.txt" %(sp[2], sp[1]))):
			if "%sVs%sTable.txt" %(sp[1], sp[2]) and "%sVs%sTable.txt" %(sp[2], sp[1]) not in BLASTlist:
				BLASTlist.append(files)
				if sp[1] not in specieslist:
					specieslist.append(sp[1])
		else:
				errorlist.append(files)
	# If there is a missing file run missing_files()
	if errorlist != []:
		missing_files(errorlist)
	print '\nTop BLAST hits have been gathered and shared hits are generated.\nThis may take a while ..Please wait..'
	## 2. Extracting top hits for each species
	# File containing the count of tophits of each species
	f_count = open(root + '/BLASThits/count_tophits.txt', "w")
	# File containing the counts of shared tophit between two species
	f_count_shared= open(root + '/BLASThits/count_sharedhits.txt', "w")
	# Saving top hits into dictionaries																	
	for files in BLASTlist:
		sp = re.split('^(.+)[Vv][Ss](.+)[Tt].*.txt$', files)
		# For spXvsspX: output top hits to as shared hits
		if sp[1] == sp[2]:
			n_hits = 0
			hit = ''
			with open(args.input_dir+ '/' + files, "r") as f_BLAST:
				# File containing BLAST information on top hits
				f_info = open(root + '/BLASThits/BLASTinfo_tophits/%sVs%sTable.txt' %(sp[1], sp[2]), "w")
				# File containing BLAST information on top hits
				f_same = open(root + '/BLASThits/shared_hits/%sVs%sTable.txt' %(sp[1], sp[2]), "w")
				line = f_BLAST.readline()
				# Extracting each new hit and print to file
				while line:
					prot = re.split('\t', line)
					# If hit already exist - read new line
					if prot[0] == hit:
						line = f_BLAST.readline()
					# If hit doesn't exist - save hit to file
					else:
						n_hits += 1
						hit = prot[0]
						f_info.write(line)
						f_same.write('%s\t%s\n' %(prot[0], prot[1]))
						line = f_BLAST.readline()
				print '  %s vs %s' %(sp[1], sp[2])
				f_count_shared.write('%s\t%s\t%s\n' % (sp[1], sp[2], n_hits))
				f_count.write('%s\t%s\t%s\n' % (sp[1], sp[2], n_hits))
				f_info.close()
				f_same.close()
		# For spXvsspX: collect top hits into dictionaries and output shared hits
		else:
			tophit = dict()
			tophit_inverse = dict()
			n_hits = 0
			hit = ''
			# File spXvsspY, extracts tophits to "tophit" dictionary
			with open(args.input_dir+ '/' + files, "r") as f_BLAST_for:
				f_info = open(root + '/BLASThits/BLASTinfo_tophits/%sVs%sTable.txt' %(sp[1], sp[2]), "w")
				line = f_BLAST_for.readline()
				while line:
					prot = re.split('\t', line)
					# If hit already exist - read new line
					if prot[0] == hit:
						line = f_BLAST_for.readline()
					# If hit doesn't exist - save hit to dictionary
					else:																													
						n_hits += 1
						hit = prot[0]
						f_info.write(line)
						# Save to dict. Column 1 = key, column 2 = value
						tophit[prot[0]] = prot[1]
						line = f_BLAST_for.readline()
				f_count.write('%s\t%s\t%s\n' % (sp[1], sp[2], n_hits))
				n_hits = 0
				hit = ''
			# File spYvsspX (INVERSE), extraxt tophits to "tophit_inverse" dictionary
			with open(args.input_dir+"/%sVs%sTable.txt" %(sp[2], sp[1]), "r") as f_BLAST_rev:
				f_info = open(root + '/BLASThits/BLASTinfo_tophits/%sVs%sTable.txt' %(sp[2], sp[1]), "w")
				line = f_BLAST_rev.readline()
				while line:
					prot = re.split('\t', line)
					# If hit already exist - read new line
					if prot[0] == hit:
						line = f_BLAST_rev.readline()
					# If hit doesn't exist - save hit to dictionary
					else:
						n_hits += 1
						hit = prot[0]
						f_info.write(line)
						# key:value variables are switched
						# Save to dict. Column 2 = key, column 1 = value
						tophit_inverse[prot[1]] = prot[0]
						line = f_BLAST_rev.readline()
				f_count.write('%s\t%s\t%s\n' % (sp[2], sp[1], n_hits))
			n_hits = 0
			# File containing the shared hits between two species
			f_shared = open(root + '/BLASThits/shared_hits/%sVs%sTable.txt' %(sp[1], sp[2]), "w")
			## 3. Create shared hit files by extracting identical key:value entries between species
			for key in tophit.keys():
				# Find all key:value entries that exists in "inverse"																				
				if key in tophit_inverse:
					if tophit[key] == tophit_inverse[key]:
						f_shared.write('%s\t%s\n' %(key, tophit[key]))
						n_hits += 1	
			f_count_shared.write('%s\t%s\t%s\n' % (sp[1], sp[2], n_hits))
			print '  %s vs %s' %(sp[1], sp[2])
			f_shared.close()
	f_count.close()
	f_count_shared.close()
	# Remove empty spaces and sort count_tophits.txt alphabetically
	with open(root + '/BLASThits/count_tophits.txt', "r") as f_count:
		# Omit empty lines and lines containing only whitespace
		lines = [line for line in f_count if line.strip()]
		f_count.close()
	# Sort lines alphabetically
	lines.sort()
	f_count_sort = open(root + '/BLASThits/count_tophits.txt', "w")
	f_count_sort.writelines(lines)
	f_count_sort.close()
	specieslist.sort()
	return specieslist, join(root,'/BLASThits/shared_hits/')



#===================================================================
# 4. Checking for missing files
#===================================================================

def missing_files(errorlist):
		# Accepted inputs
		yes = set(['yes','y', 'ye', ''])
		no = set(['no','n'])
		print "\nThese files do not have an opposite BLAST file or aren't named <speciesA>Vs<speciesB>Table.txt:"
		i = 0
		for item in errorlist:
			i += 1
			print "   %s. %s" %(i, item)
		missing = 0
		while missing == 0:
			choice = raw_input("\nDo you want to continue without these files? [Yes/No]: ").lower()
			if choice in no:
				return sys.exit("\nRerun this program, after editing the files")
			elif choice in yes:
				missing = 1
				pass
			else:
				missing = 0
				sys.stdout.write("\nPlease respond with 'yes' or 'no'")



#===================================================================
# 5. Collecting the hits into a matrix
#===================================================================

def matrix(root, shared_folder, specieslist):
	""" Output:
		1. Matrices containing shared hits between species. One for each species.
		2. Matrices containing alphabetic sorted columns
	"""
	title = 'test'
	flag = '0'
	# Generate combinations of species from specieslist
	for sp1 in specieslist:
		for sp2 in specieslist:
			# Filename: spXvsspY
			files = "%sVs%sTable.txt" %(sp1, sp2)
			# Filename: spYvsspX
			files_inverse = "%sVs%sTable.txt" %(sp2, sp1)
			## 1. Save files into dictionaries, checking if file exists else reverse dictionary keys
			# If file was generated in the step before
			# Read tab-delimited file into a table
			# Save to sp_dict. Normal key:values
			if isfile(join(shared_folder,files)):																							
				sp_file = csv.reader(open(shared_folder+files, 'r'), delimiter='\t')						
				sp_dict = dict((rows[0],rows[1]) for rows in sp_file if rows != "")							
			else:
			# If not, find opposite file and save hits into a dictionary
			# Read tab-delimited file into a table
			# Save to sp_dict. Reverse key:values
				if isfile(join(shared_folder,files_inverse)):
					sp_file = csv.reader(open(shared_folder+files_inverse, 'r'), delimiter='\t')
					# Reverse key:value variables 
					sp_dict = dict((rows[1],rows[0]) for rows in sp_file if rows != "")							
				else:
					sys.exit("\nExpected file not present: %s" %files_inverse)
			## 2. Collect the dictionaries into matrices containing shared hits between all species
			# When encountering a new species in sp1 filelist
			if sp1 != title:
				# Retrieve spXvsspX and use as start matrix										    																				
				sp_file_start = csv.reader(open(shared_folder+"%sVs%sTable.txt" %(sp1, sp1), 'r'), delimiter='\t')
				sp_dict_start = dict((rows[0],rows[1]) for rows in sp_file_start if rows != "")
				# Create new matrix
				matrix = DataFrame.from_dict(sp_dict_start, orient='index', dtype=None)
				# Rename first column to start species
				matrix.columns = [sp1]
				title = sp1
			# For every new species
			if sp1 != sp2:
				# Add new species column to matrix
				matrix[sp2] = Series(sp_dict, index = matrix.index.values)
		# Output matrix
		matrix.to_csv('%s/matrix/unsorted/%s.txt'%(root,title), sep="\t")
		# Sort matrix after alphabetic order
		matrix_sorted = matrix.sort_index(axis=1)
		# Output sorted matrix
		matrix_sorted.to_csv('%s/matrix/sorted/%s.txt'%(root,title), sep="\t")
		print "  %s" %title
	return root + '/matrix/sorted/'		



#===================================================================
# 6. Assign pattern to hits, group hits by pattern and count the groups
#===================================================================

def hit(root, sorted_folder, args):
	""" Output:
		1. Files containing hit list with assigned pattern. 
			 Pattern reprecents hits in other species; hit = 1, no hit = 0.
		2. Files containing counted pattern lists.
		3. File containing a matrix over counted pattern lists.
	"""
	flag = 0
	# Iterating through all files in sorted_folder
	for files in os.listdir(sorted_folder):
		pattern = dict()											
		# File containing assign pattern to each AA hit - pattern correlates to existance in a species
		f_pattern=open(root+'/counts_patterns/hitPattern/'+files, 'w')
		with open(sorted_folder+files, 'r') as f:
			line = f.readline()
			# Saving header containing species in the right order
			header = line
			f_pattern.write(header)
			line = f.readline()
			# Count total amount of hits per species
			count = 0
			# Extracting hit name 
			while line:
				# For every indexed hit
				# Output hit to column 1
				count += 1
				line = re.split('\s',line)									
				f_pattern.write('%s\t'%line[0])
				i = 1
				hit = ''
				# Transforming shared hits into 1 or 0
				# make pattern over column 2 to end column (exclude index column)
				while i < len(line)-1:										
					# If no hit - append 0
					if line[i] == '':
						hit = hit + str(0)
						i += 1
					# If hit - append 1
					else:
						hit = hit + str(1)
						i += 1
				f_pattern.write(hit)
				f_pattern.write('\n')
				## 2. Counting how often there is the same hit pattern
				# If hit is already present - continue counting
				if hit in pattern:											
					pattern[hit] += 1
				# Else create new pattern										
				else:
					pattern[hit] = 1
				line = f.readline()
			f_pattern.close()
		# Output pattern count
		# Sort on pattern: itemgetter(0) OR sort on count: itemgetter(1)
		sorted_pattern = sorted(pattern.iteritems(), key=operator.itemgetter(1), reverse=True)
		f_count = open(root+'/counts_patterns/hitcounts/'+files, 'w')
		f_count.write(header)
		f_count.write('total amount:\t%s\n'%int(count))
		# Output pattern vs. hit amounts
		f_count.write('\n'.join('%s\t%s' % x for x in sorted_pattern))
		f_count.close()
		# Create pattern matrix
		if flag == 0:
			# Create new matrix from dictionaries
			matrix = DataFrame.from_dict(pattern, orient='index', dtype=int)
			# Rename columns
			matrix.columns = [re.findall('(.*).txt',files)[0]]
			flag = 1
		else:
			# Create temporary matrix
			matrix_temp = DataFrame.from_dict(pattern, orient='index', dtype=int)
			# Rename columns
			matrix_temp.columns = [re.findall('(.*).txt',files)[0]]
			# Join the two matrices and combine/add the pattern indexing
			matrix = matrix.join(matrix_temp, on=None, how='outer', lsuffix='', rsuffix='', sort=False)
	# Select options to include to the pattern matrix
	column = len(matrix.columns)
	print "column length: %s" %column
	# Statistical calculations over specified columns
	if args.select == "median" or args.select == "standard" or args.select == "all":
		# median
		matrix["median"] = matrix.iloc[:,0:column].median(axis=1, skipna=True)
	if args.select == "mean" or args.select == "standard" or args.select == "all":
		# mean
		matrix["mean"] = matrix.iloc[:,0:column].mean(axis=1, skipna=True)
	if args.select == "sd" or args.select == "standard" or args.select == "all":
		# standard deviation
		matrix["sd"] = matrix.iloc[:,0:column].std(axis=1, skipna=True)
	if args.select == "pct" or args.select == "all":
		# pct SD over median
		matrix["%std_median"] = (matrix.iloc[:,0:column].std(axis=1, skipna=True)/matrix.iloc[:,0:column].median(axis=1, skipna=True))*100
		# pct SD over mean
		matrix["%RSD"] = (matrix.iloc[:,0:column].std(axis=1, skipna=True)/matrix.iloc[:,0:column].mean(axis=1, skipna=True))*100
	# Counting the species with the same pattern and adding it to the matrix to sort it
	count_dict = dict()
	for pattern in matrix.index.tolist():
		count_dict[pattern] = sum(map(int,pattern))
	matrix_hit = DataFrame.from_dict(count_dict, orient='index', dtype=None)
	matrix.insert(0, "hits", matrix_hit)
	# Assigning specific rows floats
	if args.select == "median" or args.select == "standard" or args.select == "all":
		matrix["median"] = matrix["median"].map(lambda x: '%2.1f' % x)
	if args.select == "mean" or args.select == "standard" or args.select == "all":
		matrix["mean"] = matrix["mean"].map(lambda x: '%2.2f' % x)
	if args.select == "sd" or args.select == "standard" or args.select == "all":
		matrix["sd"] = matrix["sd"].map(lambda x: '%2.2f' % x)
	if args.select == "pct" or args.select == "all":
		matrix['%std_median'] = matrix['%std_median'].map(lambda x: '%2.2f' % x)
		matrix['%RSD'] = matrix['%RSD'].map(lambda x: '%2.2f' % x)
	# Sort matrix
	matrix_sorted = matrix.sort(columns='hits', ascending=False)
	matrix_sorted.to_csv('%s/counts_patterns/countMatrix_sorted.txt'%root, sep="\t", float_format='%.0f')



#-------------------------------------------------------------------
# MAIN PROGRAM
#-------------------------------------------------------------------
def main(argv):
	args = ParseArguments()
	my_function(args)
	root = mkdir(args)
	specieslist, shared_folder = shared(args, root)
	print 'Succesfull comparison.\n\nCollecting the shared hits into matrices. .. Please wait..'
	sorted_folder = matrix(root, root+'/BLASThits/shared_hits/', specieslist)
	print 'Succesfull matrix generation.\n\nGenerating shared-hit counts matrices. .. Please wait..'
	hit(root,sorted_folder, args)
	print 'Succesfull matrix generation.'


if __name__ == '__main__':
	main(argv)
	print "\n### Your program has finished! ###\n"

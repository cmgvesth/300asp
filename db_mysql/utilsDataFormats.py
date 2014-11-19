
def check_seq_protein(seq):
	# This function takes a sequences and tests if the alphabet matches a amino acid string
	# Function does not return any values but prints warning to the screen 
	symbol_flag = 0

	if len(set(str(seq.lower())) - set("actg")) == 0	:
		print "# WARNING: found only A, C, T, and G - possible DNA sequence?"

		non_ExtendedIUPACProtein = set(str(seq.lower())) - set("ACDEFGHIKLMNPQRSTVWYBXZJUO".lower()) 
		# ExtendedIUPACProtein: ACDEFGHIKLMNPQRSTVWYBXZJUO http://biopython.org/DIST/docs/api/Bio.Alphabet.IUPAC.ExtendedIUPACProtein-class.html
		if len(non_ExtendedIUPACProtein) != 0:
			if len(non_ExtendedIUPACProtein) == 1 and "*" in non_ExtendedIUPACProtein:
				if symbol_flag != 1:
					# Only print the warning once for each file, symbol_flag
					print "# WARNING: Found ambiguouse symbol: *"
					symbol_flag = 1
			else:
				sys.exit("# ERROR : Found amino acids other than ACDEFGHIKLMNPQRSTVWYBXZJUO = %s" % ",".join(non_ExtendedIUPACProtein))



def check_seq_dna(seq):
	# This function takes a sequences and tests if the alphabet matches a DNA string
	# Function does not return any values but pritns warning to the screen 

	# Make "sure" sequence is DNA ------------------------------------ 
	non_IUPACUnambiguousDNA = set(str(seq.lower())) - set("atcg") 
	# ExtendedIUPACDNA : GATCBDSW http://biopython.org/DIST/docs/api/Bio.Alphabet.IUPAC.ExtendedIUPACDNA-class.html
	# IUPACAmbiguousDNA : GATCRYWSMKHBVDN http://biopython.org/DIST/docs/api/Bio.Alphabet.IUPAC.IUPACAmbiguousDNA-class.html
	# IUPACUnambiguousDNA : GATC http://biopython.org/DIST/docs/api/Bio.Alphabet.IUPAC.IUPACUnambiguousDNA-class.html

	if len(non_IUPACUnambiguousDNA) > 0:
		if len(non_IUPACUnambiguousDNA) == 1 and "n" in non_IUPACUnambiguousDNA:
			print "# WARNING: Found ambiguouse base: N"
		else:
			sys.exit("# ERROR : Found bases other than A, T, G, C or N: " + ",".join(non_IUPACUnambiguousDNA))

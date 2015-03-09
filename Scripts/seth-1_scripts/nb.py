from Bio import Entrez
from Bio import SeqIO
Entrez.email="setd@bio.dtu.dk"

def search_entrez(db, term):
	handle=Entrez.esearch(db, term) # pPCP1
	record = Entrez.read(handle)
	record["IdList"]
	temp_list = record["IdList"]
	print temp_list
	print Entrez.epost("nucleotide", id=','.join(temp_list)).read()
	search_result = Entrez.read(Entrez.epost("nucleotide", id=','.join(temp_list)))
	webenv = search_result["WebEnv"]
	query_key = search_result["QueryKey"]
	#return webenv, query_key

# Getting a summary of search results
# def overview(session_variables):
	#webenv, query_key = session_variables
	summary_handle = Entrez.esummary(db="nuclteotide", webenv=webenv, query_key = query_key)
	record = Entrez.read(summary_handle)
	record[0]
	counter = 1
	for n in record:
		print "Entry #" , counter
		counter +=1

		print "Title:", n["Title"]
		print "Date:" , n["CreateDate"]
		print "Caption:", n["Caption"]
		print "Length:", n["Length"]

	print "Which Entries shall be saved?"
	input_entries = input('Please provide a list of entries you want to save')

	print "Converting to IDs"
	if type(input_entries)==list:
		id_list = []
		for n in input_entries:
			id_list.append(record[n]["Id"])
		
	else:
		id_list = record[input_entries]["Id"]
	#return id_list

# EFetch part
	print "Downloading files"
	handle = Entrez.efetch(db="nucleotide", id=id_list, rettype = "gb", retmode = "text")
#	result = handle.read()

	print "Do you want to convert the file into a sequence object?(y/n)"
	answer = raw_input(">")

	if answer == "y":
		if type(input_entries)==list:
			for seq_record in SeqIO.parse(handle, "gb"):
				print seq_record.id, seq_record.description[:50] + "..."
				print "Sequence length %i," % len(seq_record),
				print "%i features," % len(seq_record.features),
				print "from: %s" % seq_record.annotations["source"]
		else:
			seq_record = SeqIO.read(handle, "gb")
			handle.close()
			print "%s with %i features" % (seq_record.id, len(seq_record.features))


	return result

	#print "Saving data"
	#save_file = open("work_in_progress.gb", "w")
	#save_file.write(handle.read())
	#save_file.close()
	#handle.close()




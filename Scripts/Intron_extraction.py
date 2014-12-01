from Bio import  SeqIO
from Bio.Seq import Seq
import os
import matplotlib.pyplot as plt

assembly_file  = ""
gff_file = ""
gene_catalog = ""
output_file = ""



#===============================================================================
# Gene class collecting the information in the gff file
#===============================================================================

class gene:
    # Assign boolean values to strand description
    lookup_reverse = {'+':False,'-':True}
    
    def __init__(self, name, supercontig, strand):
        self.name = name
        self.cds = []
        self.supercontig = supercontig
        self.reverse = self.lookup_reverse[strand]

    def add_cds(self, start, end):
        self.cds.append((start, end))
    
    def calculate_intron_positions(self):
        # Sort the coding sequences list by start position
        sorted_cds = sorted(self.cds, key=lambda x: x[0])
        intron_positions = []
        # Get the start and end position of every intron i.e. the sequence between two coding sequences
        for index in range(1,len(sorted_cds)):
            # The start position of the intron corresponds to the end position plus 1 of the previous coding sequence
            intron_start = sorted_cds[index-1][1]+1
            # The end position of the intron corresponds to the start position minus 1 of the next coding sequence
            intron_end = sorted_cds[index][0]-1
            intron_positions.append((intron_start,intron_end))
        return intron_positions
    
    def get_reverse_complement(self, seq_string):
        return str(Seq(seq_string).reverse_complement())
            
    def get_sequences(self, contigs, introns=0):
        """ Return the sequences of either introns or exons.
            If not specified CDS are returned.
        """
        sequences = []
        contig_sequence = contigs[self.supercontig]
        if introns:
            intron_positions = self.calculate_intron_positions()
            for intron in intron_positions:
                # Get DNA sequence of the introns. Note: Positions are shifted by one, because of the zero-based numbering in python,
                # as opposed to numbering starting from 1 in the gff-files! For the end position the shift cancels out the necessary 
                #increment of 1 due to the non-inclusiveness of the end-position in string slices
                sequences.append(contig_sequence[intron[0]-1:intron[1]])
        else:
            for cds in sorted(self.cds, key=lambda x: x[0]):
                sequences.append(contig_sequence[cds[0]-1:cds[1]])

        if self.reverse:
            # If the elements are located on the minus strand return the reverse complement of the sequences in inverse order
            return [self.get_reverse_complement(sequence) for sequence in reversed(sequences)]
        else:
            return sequences
            
        
        
            

#===============================================================================
# Function definitions
#===============================================================================

def parse_gff(ggf_file):
    """ Input: name of a gff-file
    Output: Dictionary of gene objects from the information provided in the file.
    The names given in the gff-file are used as keys.
    """
    genes = {}
    with open(ggf_file,'r') as file:
        for line in file:
            split = line.split('\t')
            name = split[8].strip('"')[6:16]
            supercontig = split[0]
            line_type = split[2]
            start = int(split[3])
            end = int(split[4])
            strand = split[6]

            if name not in genes:
                genes[name] = gene(name, supercontig, strand)
                
            if line_type == 'CDS':
                genes[name].add_cds(start, end)         
    return genes

def read_fasta(fasta_file):
    """ Input: Name of a fasta-file
    Output: Dictionary of the sequences converted to upper case.
    The identifiers in the fasta sequence are used as keys.    
    """
    result = {}
    for entry in SeqIO.parse(fasta_file,"fasta"):
        if entry.name.find('|') != -1:
            name = entry.name.split('|')[-1][:-4]
        else:
            name = entry.name
        result[name]=str(entry.seq).upper()
    return result

def quality_check(gene, contigs, official_genes):
    """ Input: Gene object, dictionary of super contig sequences and a
    dictionary containing the official gene predictions (without introns)
    Output: True if gene is contained in official gene list and 
    the reconstructed gene sequence matches the official one. False otherwise.
    """
    n = 0
    if gene.name in official_genes:
        if "".join(gene.get_sequences(contigs)) == official_genes[gene.name]:
            return True
        else:
            print("".join(gene.get_sequences(contigs)), official_genes[gene.name])
            for element in gene.cds:
                n = n + element[0]
            if n == len(official_genes[gene.name]):
                print("Length does not match for "+gene.name)
                print()
            return False
    else:
        return False


#===============================================================================
# Load files
#===============================================================================

genes = parse_gff(gff_file)
contigs = read_fasta(assembly_file)
official_genes = read_fasta(gene_catalog)

intron_list = []


for entry in genes:
    gene_entry = genes[entry]
    if quality_check(gene_entry, contigs, official_genes):
        for i, intron in enumerate(gene_entry.get_sequences(contigs,introns=1)):
            intron_list.append((gene_entry.name+'_'+str(i+1), intron))
    else:
        print()
        print(entry+' is not matching!')


#===============================================================================
# Print result file
#===============================================================================

with open(output_file,'w') as result_file:
    for element in sorted(intron_list):
        # Write the output file in a fasta compliant way
        # Write identifier
        result_file.write('>'+element[0]+'\n')
        # Write Sequence with 80 characters per line
        result_file.write('\n'.join(element[1][j:j+80] for j in range(0,len(element[1]), 80)))
        result_file.write('\n')

#===============================================================================
# End of execution
#===============================================================================
print("Done.")
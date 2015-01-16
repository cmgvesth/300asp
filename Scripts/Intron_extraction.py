from Bio import  SeqIO
from Bio.Seq import Seq
import sys, argparse
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA

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
        """ Return the sequences of either introns or CDS.
            If not specified CDS are returned.
        """
        sequences = []
        contig_sequence = contigs[self.supercontig]
        if introns:
            intron_positions = self.calculate_intron_positions()
            for intron in intron_positions:
                # Get DNA sequence of the introns. Note: Positions are shifted by one, because of the zero-based numbering in python,
                # as opposed to numbering starting from 1 in the gff-files! For the end position the shift cancels out the necessary 
                # increment of 1 due to the non-inclusiveness of the end-position in string slices
                sequences.append(contig_sequence[intron[0]-1:intron[1]])
        else:
            for cds in sorted(self.cds, key=lambda x: x[0]):
                sequences.append(contig_sequence[cds[0]-1:cds[1]])

        if self.reverse:
            # If the elements are located on the minus strand return the reverse complement of the sequences in inverse order
            return [self.get_reverse_complement(sequence) for sequence in reversed(sequences)]
        else:
            return sequences
            
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


def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('--gff', '-g', required=True)
    parser.add_argument('--contigs', '-c', required=True)
    parser.add_argument('--output', '-o', required=True)
    args = parser.parse_args()
    output_file = args.output
    assembly_file = args.contigs
    gff_file = args.gff
      
    # Read files
    genes = parse_gff(gff_file)
    contigs = read_fasta(assembly_file)
    # Get introns
    intron_list = []
    for entry in sorted(genes.keys()):
        gene_entry = genes[entry]
        for i, intron in enumerate(gene_entry.get_sequences(contigs,introns=1)):
            intron_list.append(SeqRecord(Seq(intron, IUPACAmbiguousDNA), id=entry+'.'+str(i+1), description=""))
              
    n = SeqIO.write(intron_list, output_file, "fasta")
    print("{0} introns printed in {1}".format(str(n), output_file))


if __name__ == "__main__":
    main(sys.argv[1:])
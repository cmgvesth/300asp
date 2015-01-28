#!/usr/bin/python

'''##################################################################
# Imports
##################################################################'''
import sys 
sys.path.append('../utils/')
from aspmine_imports import *


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


"""
Basic data table query followed by sub queries to get sequences or annotation 

Types of options:

Search pattern:
-- gene name (from what table?)
-- protein id
-- transcipt id

Method of search:
-- simple extract DNA, amino acids, annotations
-- one directional best matches to search (from BLAST table)
-- bi-directional best matches to search (from BLAST table)
-- all matches with same annotaton (interpro, go, PFAM)

Output data:
-- FSA
-- FNA
-- FSA, with BLAST score to search
-- FNA, with BLAST score to search
-- Annotation table (TAB)
-- BLAST table, scores to search (TAB)

if there IS a bidirectional hit (does not need to be the best, but any hit)
Sometimes there is not, ususally if we have not run the comparison, but then it becomes difficult to reduce the table by joining with itself
If no bidirectional hit - use one directional, but how do I do this in a query???? 

if gene_name:


# Method of search

if extractOnly:
	execute query for extract

if 	
"""



"""
-- input : blast qseq id
--- how much will I have to make this adjustable? How many %'s and so on
--- will python-mysql allow me to create tables? I want these to be temporary
create table rep as
select a.* from blast a
left join blast b
on b.blast_qseq_id=a.blast_sseq_id
and b.blast_sseq_id=a.blast_qseq_id
and b.blast_qseq_id<>b.blast_sseq_id
where (a.blast_qseq_jg3 like "%[insert pattern here]%" or a.blast_sseq_jg3 like "%[insert pattern here]%") 
and (b.blast_qseq_id is null or b.blast_qseq_id >= a.blast_qseq_id);

-- output : table with the maximum scoring hit for each organism to search.
--- this is NOT the best reciprocal hit of the traget in the organism 
--- ACTUALLY I am not sure, could this actually be the reciprocal? 
--- what would the table look like if geneA <- geneB = 97% and geneA -> geneD = 99%
--- maybe it is the largest between the two!
--- it does not EXCLUDE the match if it is not reciprocal 

create table rep2 as 
select * from rep r
where blast_pident=(select max(cast( blast_pident as DECIMAL(5,2))) from rep s
where s.blast_qseq_jg1=r.blast_qseq_jg1 and s.blast_sseq_jg1=r.blast_sseq_jg1  
group by s.blast_sseq_jg1,s.blast_qseq_jg1);

-- get list of sequence identifyers from this
 select blast_sseq_id from rep2 r union select blast_qseq_id from rep2 s;

-- get sequences for these
--- this returns 30 out of the 47 sequences
select * from proteins where prot_seqname in (select blast_sseq_id from rep2 r union select blast_qseq_id from rep2 s)

--- the difference 
select blast_sseq_id from (select blast_sseq_id from rep2 r union select blast_qseq_id from rep2 s) as a where blast_sseq_id not in ( select prot_seqname from proteins where prot_seqname in (select blast_sseq_id from rep2 r union select blast_qseq_id from rep2 s));
+----------------------------------------------------------------------+
| blast_sseq_id                                                        |
+----------------------------------------------------------------------+
| jgi|Eurhe1|381848|gm1.8563_g                                         |
| jgi|Neofi1|9470|7000001157007457                                     |
| jgi|PenbrAgRF18_1|322598|fgenesh1_kg.3_#_1062_#_Locus255v1rpkm594.95 |
| jgi|Pench1|38193|fgenesh1_pm.6_#_197                                 |
| jgi|Pendi1|3889|PDIG_46270m.01                                       |
| NRRL3_03617                                                          |
| Aspve1_0086035                                                       |
| Afu6g04740                                                           |
| ATET_06973                                                           |
| Aspgl1_0027831                                                       |
| Aspzo1_0019096                                                       |
| Aspsy1_0051743                                                       |
| Aacu16872_046190                                                     |
| AO090701000065                                                       |
| Aspbr1_0046489                                                       |
| Aspfo1_0033769                                                       |
| Asptu1_0048720                                                       |
+----------------------------------------------------------------------+
"""
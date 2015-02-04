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
								description="Get best hits from specific gene name or ID", 
								 usage='%(prog)s -dbname [database name] -out [output filename] -gene [JGI gene name or ID] -strain [JGI organism keys] \n'
								 "Example: python %(prog)s -dbname aspminedb -out hitlist.csv -gene An15g00560 -strain Aspni_DSM_1")
parser.add_argument("-out", "-o", required=False, default = "hitlist.csv", help="Name of output file, default hitlist.csv")
parser.add_argument("-tab", "-t", required=False, action='store_true', help="Create table new or use existing tmp")
parser.add_argument("-dbname", "-d", required=False, default = "aspminedb", help="Database name")
parser.add_argument("-strain", "-s", required=True, default=[None], action='store', help="JGI organism keys")
parser.add_argument("-gene", "-g",  required=True, default=[None], action='store', help="List of gene ID or name")

args 	= parser.parse_args()
gene 	= args.gene
strain = args.strain
outfile = args.out
table 	= args.tab

'''------------------------------------------------------------------
# Print argument values to screen
------------------------------------------------------------------'''
print "#--------------------------------------------------------------\n\
# ARGUMENTS\n\
# Description\t: Script return list of best bidirectional hit for a gene in a specific strain\n\
# Database\t: %s\n\
# Outfile\t: %s\n\
# Gene\t\t: %s\n\
# Strain\t: %s\n\
# Table\t\t: %s\n\
#--------------------------------------------------------------" % (args.dbname, outfile, gene, strain, table)

'''------------------------------------------------------------------
# Connect to specific DB
------------------------------------------------------------------'''
try:
    db = mdb.connect("localhost","asp","1234",str(args.dbname))
except mdb.Error, e:
    sys.exit("# ERROR %d: %s" % (e.args[0],e.args[1]))

try:
	cursor = db.cursor()
except mdb.Error, e: 
	sys.exit("# ERROR %d: %s" % (e.args[0],e.args[1]))

##################################################################
# FUNCTIONS
##################################################################

##################################################################
def make_database(strain_string, gene_string):
##################################################################
	print "# INFO: Creating subset table"

	strain_string	= strain_string.replace("blast", "a.blast")
	gene_string		= gene_string.replace("blast", "a.blast")
	
	'''------------------------------------------------------------------
	# This table includes two way hits, as it gets both sseq and qseq data
	# These hits could be the same but in rare cases they might differ
	# Common case: A -> B AND B -> A # not nessesaraly with the same score
	# Rare case: A -> B BUT B ! -> ! A
	------------------------------------------------------------------'''

	try:
		cursor.execute("DROP table IF EXISTS tmp ;")
		cursor.execute("CREATE table tmp as \
						SELECT a.* FROM blast a \
						LEFT JOIN blast b \
						ON b.blast_qseq_id=a.blast_sseq_id \
						AND b.blast_sseq_id=a.blast_qseq_id \
						AND b.blast_qseq_id!=b.blast_sseq_id \
						WHERE ( " + gene_string + ") \
						AND (" + strain_string + ")\
						AND (b.blast_qseq_id IS NULL OR b.blast_qseq_id > a.blast_qseq_id);")

	except mdb.Error, e:
		sys.exit("# ERROR: executing following query failed:\n%s\n # ERROR %d: %s" % (cursor._last_executed, e.args[0],e.args[1]))

	cursor.execute("SELECT count(*) FROM tmp")
	result = cursor.fetchall()

	if not type(result) is tuple:
		sys.exit("# ERROR : result of query is not a tuple but a %s" % type(result))

	print "# INFO: Created table tmp using gene string:\n# %s\n# and strain string:\n# %s" % (gene_string, strain_string)
	"""
	SELECT *  FROM `tmp` WHERE `blast_qseq_jg1` LIKE 'Aspacu1' OR `blast_sseq_jg1` LIKE 'Aspacu1'
	ORDER BY `tmp`.`blast_pident` ASC

	SELECT *  FROM `blast` WHERE 
	(`blast_qseq_jg1` LIKE 'Aspacu1' OR `blast_sseq_jg1` LIKE 'Aspacu1') and 
	(`blast_qseq_jg1` LIKE 'Aspni7' OR `blast_sseq_jg1` LIKE 'Aspni7') and 
	(`blast_qseq_jg2`="1123159" or `blast_sseq_jg2`="1123159")
	ORDER BY `blast_pident` ASC

	CREATE table tmp as 
	SELECT a.blast_sseq_id, a.blast_qseq_id, a.blast_pident FROM blast as a 
	LEFT JOIN blast as b 
	ON b.blast_qseq_id=a.blast_sseq_id AND 
	b.blast_sseq_id=a.blast_qseq_id AND 
	b.blast_qseq_id<>b.blast_sseq_id WHERE 
	(a.blast_qseq_id like '%1123159%' or a.blast_sseq_id like '%1123159%') AND   
	(a.blast_qseq_jg1 = 'Aspni7' or a.blast_sseq_jg1 = 'Aspni7') AND
	(a.blast_qseq_jg1 LIKE 'Aspacu1' OR a.blast_sseq_jg1 LIKE 'Aspacu1') AND
	(b.blast_qseq_id is null OR b.blast_qseq_id >= a.blast_qseq_id);
	"""		
##################################################################
def get_data ():
##################################################################
	result = ""

	try: 
		cursor.execute("DROP table IF EXISTS tmp_rep ;")
		cursor.execute("CREATE TABLE tmp_rep as \
			SELECT blast_sseq_id,blast_qseq_id,blast_pident, blast_qseq_jg1, blast_sseq_jg1 ,(blast_qend-blast_qstart)/blast_qlen*100 as qcov, (blast_send-blast_sstart)/blast_slen*100 as scov \
			from tmp t  where blast_qseq_jg1 != blast_sseq_jg1 and  \
			blast_pident in ( select greatest(maxa,maxb) as pairmax from  \
			(select blast_qseq_jg1, blast_sseq_jg1, max(blast_pident) as maxa from tmp as a group by blast_qseq_jg1,blast_sseq_jg1) as a \
			join  \
			(select blast_qseq_jg1, blast_sseq_jg1, max(blast_pident) as maxb from tmp as a group by blast_sseq_jg1,blast_qseq_jg1) as b \
			on b.blast_qseq_jg1=a.blast_sseq_jg1 and b.blast_sseq_jg1=a.blast_qseq_jg1 \
			WHERE (t.blast_qseq_jg1 = a.blast_qseq_jg1 or  t.blast_qseq_jg1 = a.blast_sseq_jg1) and \
			(t.blast_sseq_jg1 = a.blast_sseq_jg1 or t.blast_sseq_jg1 = a.blast_qseq_jg1));")

		cursor.execute("SELECT * from tmp_rep")

	except mdb.Error, e:
		sys.exit("# ERROR: executing following query failed:\n%s\n # ERROR %d: %s" % (cursor._last_executed, e.args[0],e.args[1]))
	
	result = cursor.fetchall()

	if not type(result) is tuple:
		sys.exit("# ERROR : result of query is not a tuple but a %s" % type(result))

	if len(result) < 1:
		sys.exit("# ERROR: query returned no rows, verify strains and gene names")	
	#print "# INFO: get_data returns %s rows from tmp" % len(result)	

	return result

##################################################################
def test_table(strain_string, gene_string):
##################################################################
	'''------------------------------------------------------------------
	# Check if table tmp already contains the needed genesa and strains 
	# The user might think that the table is already there.
	# For the case where that is wrong, error and tell to run new table
	------------------------------------------------------------------'''
	try: 
		cursor.execute("SELECT blast_sseq_id,blast_qseq_id FROM tmp WHERE ( " + gene_string + ") AND (" + strain_string + ") ;")
	except mdb.Error, e:
			sys.exit("# ERROR: executing following query failed:\n# %s\n# ERROR %d: %s" % (cursor._last_executed, e.args[0],e.args[1]))

	result = cursor.fetchall()

	if len(result) < 1:
		if not table:
			sys.exit("# ERROR: table tmp does not contain the needed genes/strains\n# Run with -t option to redo the table")
		else:
			print "# INFO: table tmp does not contain the needed genes/strains"	
	else:
		print "# INFO: table tmp was found to contain entries for strain and genes"		
			
##################################################################
# MAIN
##################################################################

'''------------------------------------------------------------------
# Create strain and gene query strings
------------------------------------------------------------------'''
#(strain_string, gene_string) = ("","")

strain_string = "blast_qseq_jg1 = \'" + strain + "\' or blast_sseq_jg1 = \'" + strain + "\'"
gene_string = "blast_qseq_id like \'%" + gene + "%\' or blast_sseq_id like \'%" + gene + "%\'"


'''------------------------------------------------------------------
# TEST content of existing tmp table
# create table if tmp does not contain the right genes and strains
------------------------------------------------------------------'''
if table:
	#print "# INFO: Creating table tmp"
	make_database(strain_string, gene_string)

if cursor.execute("SHOW TABLES LIKE 'tmp' ;"):
	test_table(strain_string, gene_string)
else:
	print "# INFO: Table tmp does not exist, will create table"
	make_database(strain_string, gene_string)

'''------------------------------------------------------------------
# Get data from tmp table
------------------------------------------------------------------'''
result = get_data()

f = open(outfile,'wb')
writer = csv.writer(f, dialect = 'excel', delimiter=";")
columns = map(lambda x:x[0], cursor.description) 	
writer.writerow(columns)
writer.writerows(result)
f.close()

'''------------------------------------------------------------------
# Get seqs
------------------------------------------------------------------'''
try:
	cursor.execute("select distinct prot_seqname, prot_seq from proteins join (\
			select \
			IF (LOCATE('Aacu16872_0', seqid)>0 and LOCATE('|', seqid)<1, SUBSTRING(seqid, LOCATE('Aacu16872_0', seqid)+11), seqid) AS seqid, \
			IF (orgid = '', 'NRRL3', REPLACE(orgid, 'NRRL3', 'Aspni_NRRL3_1')) AS orgid \
			from (select distinct *\
			from (select blast_sseq_id as seqid, blast_sseq_jg1 as orgid from tmp_rep \
				union all \
				select blast_qseq_id as seqid,blast_qseq_jg1 as orgid from tmp_rep) as a) as b) as d\
			on( prot_seqname like concat( '%', seqid, '%') and prot_orgkey like concat(orgid,'%')\
			) order by prot_seqname;")

except mdb.Error, e:
	sys.exit("# ERROR: executing following query failed:\n%s\n # ERROR %d: %s" % (cursor._last_executed, e.args[0],e.args[1]))

seqdata = cursor.fetchall()
f = open("hitseqs.csv",'wb')
writer = csv.writer(f, dialect = 'excel', delimiter=";")
columns = map(lambda x:x[0], cursor.description) 	
writer.writerow(columns)
writer.writerows(seqdata)
f.close()


os.system("sed 's/^/>/' hitseqs.csv | sed 's/;/\\n/' | grep -v prot > hitseqs.fsa")
os.system("rm hitseqs.csv")

'''------------------------------------------------------------------
# Final information to user
------------------------------------------------------------------'''
nr_rows = len(result)

percent = [float(row[2]) for row in result ]
mean_id = float(sum(percent))/len(percent) if len(percent) > 0 else float('nan')

min_id	= min(percent)
max_id	= max(percent)

min_string = [row[0:2] for row in result if float(row[2]) == min_id]
max_string = [row[0:2] for row in result if float(row[2]) == max_id]

min_string = ["%s -- VS -- %s" % tup for tup in min_string]
max_string = ["%s -- VS -- %s" % tup for tup in max_string]

print "#--------------------------------------------------------------\n\
# ARGUMENTS\n\
# Runtime\t\t: %s\n\
# Nr. rows in table\t: %s\n\
# Average percent id\t: %s\n\
# Maximum percent id\t: %s \n\
# Comparisons with max\t: %s\n\
# Minimum percent id\t: %s \n\
# Comparisons with min\t: %s\n\
#--------------------------------------------------------------" % \
( (datetime.now()-startTime), nr_rows, mean_id, max_id, '\n#\t\t\t: '.join(map(str, max_string)), min_id, '\n#\t\t\t: '.join(map(str, min_string)))




'''------------------------------------------------------------------
# NOTES
------------------------------------------------------------------'''

"""

select IF (seqid = '', '_0', REPLACE(seqid, '_0', '%')) AS seqid, orgid  from (select distinct * from (select blast_sseq_id as seqid, blast_sseq_jg1 as orgid from tmp_rep union all select blast_qseq_id as seqid,blast_qseq_jg1 as orgid from tmp_rep) as a) as b;
select 
if(LOCATE('_0', seqid)>0 and LOCATE('|', seqid)<1, SUBSTRING(seqid, LOCATE('_0', seqid)+2), seqid) AS seqid, 
orgid  from (
	select distinct * from (
		select blast_sseq_id as seqid, blast_sseq_jg1 as orgid from tmp_rep 
		union all 
		select blast_qseq_id as seqid,blast_qseq_jg1 as orgid from tmp_rep) 
	as a) 
as b

IF (seqid = '', '_0', REPLACE(seqid, '_', '%')) AS seqid )

select distinct(prot_seqname) from proteins join (
select 
IF (LOCATE('Aacu16872_0', seqid)>0 and LOCATE('|', seqid)<1, SUBSTRING(seqid, LOCATE('Aacu16872_0', seqid)+11), seqid) AS seqid, 
IF (orgid = '', 'NRRL3', REPLACE(orgid, 'NRRL3', 'Aspni_NRRL3_1')) AS orgid from (select distinct * from (select blast_sseq_id as seqid, blast_sseq_jg1 as orgid from tmp_rep union all select blast_qseq_id as seqid,blast_qseq_jg1 as orgid from tmp_rep) as a) as b) as d
on( prot_seqname like concat( "%", seqid, "%") and prot_orgkey like concat(orgid,"%")
) order by prot_seqname;




or (prot_orgkey like concat(seqid,"%")
and ( seqid=prot_seqname or seqid=prot_tailkey 	or prot_seqname like concat( "%", seqid, "%") ))

SELECT blast_sseq_id,blast_qseq_id,blast_pident, blast_qseq_jg1, blast_sseq_jg1 ,
(blast_qend-blast_qstart)/blast_qlen*100 as qcov, (blast_send-blast_sstart)/blast_slen*100 as scov
FROM `tmp` WHERE `blast_qseq_jg1` LIKE 'Aspnov1' OR `blast_sseq_jg1` LIKE 'Aspnov1'
ORDER BY `tmp`.`blast_pident` ASC

for s in range(0,len(strains)):
	strain_string += "blast_qseq_jg1 = \'" + strains[s] + "\' or blast_sseq_jg1 = \'" + strains[s] + "\'"
	if s != len(strains)-1:
		strain_string += " or "


for i in range(0,len(genes)):
	gene_string += "blast_qseq_id like \'%" + genes[i] + "%\' or blast_sseq_id like \'%" + genes[i] + "%\'"
	if i != len(genes)-1:
		gene_string += " or "


SELECT r.blast_sseq_id,r.blast_qseq_id,r.blast_pident, r.blast_qseq_jg1, r.blast_sseq_jg1, 
(r.blast_qend-r.blast_qstart)/r.blast_qlen*100 as qcov, (r.blast_send-r.blast_sstart)/r.blast_slen*100 as scov
from tmp r LEFT JOIN tmp s ON 
s.blast_qseq_id=r.blast_sseq_id and 
s.blast_sseq_id=r.blast_qseq_id and 
s.blast_sseq_id<>s.blast_qseq_id 
WHERE (s.blast_sseq_id IS NULL or s.blast_sseq_id > r.blast_sseq_id) 
and r.blast_pident = (
select max(blast_pident) from tmp a 
where a.blast_sseq_id=r.blast_sseq_id and 
a.blast_qseq_id=r.blast_qseq_id 
group by a.blast_qseq_jg1, a.blast_sseq_jg1)
GROUP BY r.blast_sseq_jg1,r.blast_qseq_jg1;




SELECT r.blast_sseq_id,r.blast_qseq_id,r.blast_pident, r.blast_qseq_jg1, r.blast_sseq_jg1, \
(r.blast_qend-r.blast_qstart)/r.blast_qlen*100 as qcov, (r.blast_send-r.blast_sstart)/r.blast_slen*100 as scov\
from tmp r LEFT JOIN tmp s ON \
s.blast_qseq_id=r.blast_sseq_id and \
s.blast_sseq_id=r.blast_qseq_id and \
s.blast_sseq_id<>s.blast_qseq_id \
WHERE (s.blast_sseq_id IS NULL or s.blast_sseq_id > r.blast_sseq_id) \
and r.blast_pident = (\
select max(blast_pident) from tmp a where group by a.blast_qseq_jg1, a.blast_sseq_jg1)\
GROUP BY r.blast_sseq_jg1,r.blast_qseq_jg1;


SELECT r.blast_sseq_jg1,r.blast_qseq_jg1,max(r.blast_pident), s.blast_sseq_jg1,s.blast_qseq_jg1,max(s.blast_pident)
from tmp r LEFT JOIN tmp s ON \
s.blast_qseq_id=r.blast_sseq_id and \
s.blast_sseq_id=r.blast_qseq_id and \
s.blast_sseq_id<>s.blast_qseq_id \
#WHERE (s.blast_sseq_id IS NULL or s.blast_sseq_id > r.blast_sseq_id)
GROUP BY r.blast_sseq_jg1,r.blast_qseq_jg1;


SELECT r.blast_sseq_id,r.blast_qseq_id, s.blast_sseq_id,s.blast_qseq_id
from tmp r JOIN tmp s ON \
s.blast_qseq_id=r.blast_sseq_id and \
s.blast_sseq_id=r.blast_qseq_id and \
WHERE (s.blast_sseq_id IS NOT NULL);


SELECT blast_sseq_id,blast_qseq_id,blast_pident, blast_qseq_jg1, blast_sseq_jg1 from
tmp a where blast_qseq_jg1 != blast_sseq_jg1 AND 
blast_pident = (SELECT max(b.blast_pident) from tmp b WHERE
b.blast_qseq_jg1=a.blast_qseq_jg1 AND b.blast_sseq_jg1=a.blast_sseq_jg1 AND 
b.blast_sseq_jg1<> b.blast_qseq_jg1)





Get maximum value for each organism pair:
select a.blast_qseq_jg1,a.blast_sseq_jg1,greatest(maxa,maxb) as pairmax from 
(select blast_qseq_jg1,blast_sseq_jg1,max(blast_pident) as maxa from tmp as a group by blast_qseq_jg1,blast_sseq_jg1) as a 
join 
(select blast_qseq_jg1,blast_sseq_jg1,max(blast_pident) as maxb from tmp as a group by blast_sseq_jg1,blast_qseq_jg1) as b 
on (b.blast_qseq_jg1=a.blast_sseq_jg1 and b.blast_sseq_jg1=a.blast_qseq_jg1);




SELECT blast_sseq_id,blast_qseq_id,blast_pident, blast_qseq_jg1, blast_sseq_jg1 ,(blast_qend-blast_qstart)/blast_qlen*100 as qcov, (blast_send-blast_sstart)/blast_slen*100 as scov
from tmp t 
where blast_qseq_jg1 != blast_sseq_jg1 and 
blast_pident in ( select greatest(maxa,maxb) as pairmax from 
(select blast_qseq_jg1,blast_sseq_jg1,max(blast_pident) as maxa from tmp as a group by blast_qseq_jg1,blast_sseq_jg1) as a 
join 
(select blast_qseq_jg1,blast_sseq_jg1,max(blast_pident) as maxb from tmp as a group by blast_sseq_jg1,blast_qseq_jg1) as b 
on b.blast_qseq_jg1=a.blast_sseq_jg1 and b.blast_sseq_jg1=a.blast_qseq_jg1
WHERE (t.blast_qseq_jg1 = a.blast_qseq_jg1 or  t.blast_qseq_jg1 = a.blast_sseq_jg1) and 
(t.blast_sseq_jg1 = a.blast_sseq_jg1 or t.blast_sseq_jg1 = a.blast_qseq_jg1))


select a.blast_qseq_jg1 as qorg,a.blast_sseq_jg1 as sorg,greatest(maxa,maxb) as pairmax from 
(select blast_qseq_jg1,blast_sseq_jg1,max(blast_pident) as maxa from tmp as a group by blast_qseq_jg1,blast_sseq_jg1) as a 
join 
(select blast_qseq_jg1,blast_sseq_jg1,max(blast_pident) as maxb from tmp as a group by blast_sseq_jg1,blast_qseq_jg1) as b 
on b.blast_qseq_jg1=a.blast_sseq_jg1 and b.blast_sseq_jg1=a.blast_qseq_jg1
"""
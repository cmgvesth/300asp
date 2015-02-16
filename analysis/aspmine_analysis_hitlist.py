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
								 usage='%(prog)s -dbname [database name] -out [output filename] -gene [JGI gene name or ID] -strain [JGI organism key] \n'
								 "Example: python %(prog)s -dbname aspminedb -out hitlist.csv -gene An15g00560 -strain Aspni_DSM_1")
parser.add_argument("-out", "-o", required=False, default = "hitlist.csv", help="Name of output file, default hitlist.csv")
parser.add_argument("-tab", "-t", required=False, action='store_true', help="Create table new or use existing tmp")
parser.add_argument("-dbname", "-d", required=False, default = "aspminedb", help="Database name")
parser.add_argument("-strain", "-s", required=True, default=[None], action='store', help="JGI organism key")
parser.add_argument("-gene", "-g",  required=True, default=[None], action='store', help="Gene ID or name")

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
						ON b.q_seqid=a.h_seqid \
						AND b.h_seqid=a.q_seqid \
						AND b.q_seqid!=b.h_seqid \
						WHERE ( " + gene_string + ") \
						AND (" + strain_string + ")\
						AND (b.q_seqid IS NULL OR b.q_seqid > a.q_seqid);")

	except mdb.Error, e:
		sys.exit("# ERROR: executing following query failed:\n%s\n # ERROR %d: %s" % (cursor._last_executed, e.args[0],e.args[1]))

	cursor.execute("SELECT count(*) FROM tmp")
	result = cursor.fetchall()

	if not type(result) is tuple:
		sys.exit("# ERROR : result of query is not a tuple but a %s" % type(result))

	print "# INFO: Created table tmp using gene string:\n# %s\n# and strain string:\n# %s" % (gene_string, strain_string)
	"""
	SELECT *  FROM `tmp` WHERE `q_org` LIKE 'Aspacu1' OR `h_org` LIKE 'Aspacu1'
	ORDER BY `tmp`.`pident` ASC

	SELECT *  FROM `blast` WHERE 
	(`q_org` LIKE 'Aspacu1' OR `h_org` LIKE 'Aspacu1') and 
	(`q_org` LIKE 'Aspni7' OR `h_org` LIKE 'Aspni7') and 
	(`q_seqkey`="1123159" or `h_seqkey`="1123159")
	ORDER BY `pident` ASC

	CREATE table tmp as 
	SELECT a.h_seqid, a.q_seqid, a.pident FROM blast as a 
	LEFT JOIN blast as b 
	ON b.q_seqid=a.h_seqid AND 
	b.h_seqid=a.q_seqid AND 
	b.q_seqid<>b.h_seqid WHERE 
	(a.q_seqid like '%1123159%' or a.h_seqid like '%1123159%') AND   
	(a.q_org = 'Aspni7' or a.h_org = 'Aspni7') AND
	(a.q_org LIKE 'Aspacu1' OR a.h_org LIKE 'Aspacu1') AND
	(b.q_seqid is null OR b.q_seqid >= a.q_seqid);
	"""		
##################################################################
def get_data ():
##################################################################
	result = ""

	try: 
		cursor.execute("DROP table IF EXISTS tmp_rep ;")
		cursor.execute("CREATE TABLE tmp_rep as \
			SELECT h_seqid,q_seqid,pident, q_org, h_org ,(q_end-q_start)/q_len*100 as qcov, (h_end-h_start)/h_len*100 as scov \
			from tmp t  where q_org != h_org and  \
			pident in ( select greatest(maxa,maxb) as pairmax from  \
			(select q_org, h_org, max(pident) as maxa from tmp as a group by q_org,h_org) as a \
			join  \
			(select q_org, h_org, max(pident) as maxb from tmp as a group by h_org,q_org) as b \
			on b.q_org=a.h_org and b.h_org=a.q_org \
			WHERE (t.q_org = a.q_org or  t.q_org = a.h_org) and \
			(t.h_org = a.h_org or t.h_org = a.q_org));")

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
		cursor.execute("SELECT h_seqid,q_seqid FROM tmp WHERE ( " + gene_string + ") AND (" + strain_string + ") ;")
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

strain_string = "q_org = \'" + strain + "\' or h_org = \'" + strain + "\'"
gene_string = "q_seqid like \'%" + gene + "%\' or h_seqid like \'%" + gene + "%\'"


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
			from (select h_seqid as seqid, h_org as orgid from tmp_rep \
				union all \
				select q_seqid as seqid,q_org as orgid from tmp_rep) as a) as b) as d\
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

select IF (seqid = '', '_0', REPLACE(seqid, '_0', '%')) AS seqid, orgid  from (select distinct * from (select h_seqid as seqid, h_org as orgid from tmp_rep union all select q_seqid as seqid,q_org as orgid from tmp_rep) as a) as b;
select 
if(LOCATE('_0', seqid)>0 and LOCATE('|', seqid)<1, SUBSTRING(seqid, LOCATE('_0', seqid)+2), seqid) AS seqid, 
orgid  from (
	select distinct * from (
		select h_seqid as seqid, h_org as orgid from tmp_rep 
		union all 
		select q_seqid as seqid,q_org as orgid from tmp_rep) 
	as a) 
as b

IF (seqid = '', '_0', REPLACE(seqid, '_', '%')) AS seqid )

select distinct(prot_seqname) from proteins join (
select 
IF (LOCATE('Aacu16872_0', seqid)>0 and LOCATE('|', seqid)<1, SUBSTRING(seqid, LOCATE('Aacu16872_0', seqid)+11), seqid) AS seqid, 
IF (orgid = '', 'NRRL3', REPLACE(orgid, 'NRRL3', 'Aspni_NRRL3_1')) AS orgid from (select distinct * from (select h_seqid as seqid, h_org as orgid from tmp_rep union all select q_seqid as seqid,q_org as orgid from tmp_rep) as a) as b) as d
on( prot_seqname like concat( "%", seqid, "%") and prot_orgkey like concat(orgid,"%")
) order by prot_seqname;




or (prot_orgkey like concat(seqid,"%")
and ( seqid=prot_seqname or seqid=prot_tailkey 	or prot_seqname like concat( "%", seqid, "%") ))

SELECT h_seqid,q_seqid,pident, q_org, h_org ,
(q_end-q_start)/q_len*100 as qcov, (h_end-h_start)/h_len*100 as scov
FROM `tmp` WHERE `q_org` LIKE 'Aspnov1' OR `h_org` LIKE 'Aspnov1'
ORDER BY `tmp`.`pident` ASC

for s in range(0,len(strains)):
	strain_string += "q_org = \'" + strains[s] + "\' or h_org = \'" + strains[s] + "\'"
	if s != len(strains)-1:
		strain_string += " or "


for i in range(0,len(genes)):
	gene_string += "q_seqid like \'%" + genes[i] + "%\' or h_seqid like \'%" + genes[i] + "%\'"
	if i != len(genes)-1:
		gene_string += " or "


SELECT r.h_seqid,r.q_seqid,r.pident, r.q_org, r.h_org, 
(r.q_end-r.q_start)/r.q_len*100 as qcov, (r.h_end-r.h_start)/r.h_len*100 as scov
from tmp r LEFT JOIN tmp s ON 
s.q_seqid=r.h_seqid and 
s.h_seqid=r.q_seqid and 
s.h_seqid<>s.q_seqid 
WHERE (s.h_seqid IS NULL or s.h_seqid > r.h_seqid) 
and r.pident = (
select max(pident) from tmp a 
where a.h_seqid=r.h_seqid and 
a.q_seqid=r.q_seqid 
group by a.q_org, a.h_org)
GROUP BY r.h_org,r.q_org;




SELECT r.h_seqid,r.q_seqid,r.pident, r.q_org, r.h_org, \
(r.q_end-r.q_start)/r.q_len*100 as qcov, (r.h_end-r.h_start)/r.h_len*100 as scov\
from tmp r LEFT JOIN tmp s ON \
s.q_seqid=r.h_seqid and \
s.h_seqid=r.q_seqid and \
s.h_seqid<>s.q_seqid \
WHERE (s.h_seqid IS NULL or s.h_seqid > r.h_seqid) \
and r.pident = (\
select max(pident) from tmp a where group by a.q_org, a.h_org)\
GROUP BY r.h_org,r.q_org;


SELECT r.h_org,r.q_org,max(r.pident), s.h_org,s.q_org,max(s.pident)
from tmp r LEFT JOIN tmp s ON \
s.q_seqid=r.h_seqid and \
s.h_seqid=r.q_seqid and \
s.h_seqid<>s.q_seqid \
#WHERE (s.h_seqid IS NULL or s.h_seqid > r.h_seqid)
GROUP BY r.h_org,r.q_org;


SELECT r.h_seqid,r.q_seqid, s.h_seqid,s.q_seqid
from tmp r JOIN tmp s ON \
s.q_seqid=r.h_seqid and \
s.h_seqid=r.q_seqid and \
WHERE (s.h_seqid IS NOT NULL);


SELECT h_seqid,q_seqid,pident, q_org, h_org from
tmp a where q_org != h_org AND 
pident = (SELECT max(b.pident) from tmp b WHERE
b.q_org=a.q_org AND b.h_org=a.h_org AND 
b.h_org<> b.q_org)





Get maximum value for each organism pair:
select a.q_org,a.h_org,greatest(maxa,maxb) as pairmax from 
(select q_org,h_org,max(pident) as maxa from tmp as a group by q_org,h_org) as a 
join 
(select q_org,h_org,max(pident) as maxb from tmp as a group by h_org,q_org) as b 
on (b.q_org=a.h_org and b.h_org=a.q_org);




SELECT h_seqid,q_seqid,pident, q_org, h_org ,(q_end-q_start)/q_len*100 as qcov, (h_end-h_start)/h_len*100 as scov
from tmp t 
where q_org != h_org and 
pident in ( select greatest(maxa,maxb) as pairmax from 
(select q_org,h_org,max(pident) as maxa from tmp as a group by q_org,h_org) as a 
join 
(select q_org,h_org,max(pident) as maxb from tmp as a group by h_org,q_org) as b 
on b.q_org=a.h_org and b.h_org=a.q_org
WHERE (t.q_org = a.q_org or  t.q_org = a.h_org) and 
(t.h_org = a.h_org or t.h_org = a.q_org))


select a.q_org as qorg,a.h_org as sorg,greatest(maxa,maxb) as pairmax from 
(select q_org,h_org,max(pident) as maxa from tmp as a group by q_org,h_org) as a 
join 
(select q_org,h_org,max(pident) as maxb from tmp as a group by h_org,q_org) as b 
on b.q_org=a.h_org and b.h_org=a.q_org
"""
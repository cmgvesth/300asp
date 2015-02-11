import sys,os
sys.path.append(os.path.join(os.path.dirname(__file__), "../../utils/"))
from aspmine_imports import *

#------------------------------------------------------------------
''' Get command line arguments '''
#------------------------------------------------------------------
parser = argparse.ArgumentParser(description="Preprocessing databases for gene clustering", usage="%(prog)s -out filename")
parser.add_argument("-test","-t", required=False, action="store_true", help="Enables testing mode")
parser.add_argument("-out", "-o", required=False, default = "initial.csv", help="Name of output file")
parser.add_argument("-limit", "-lim", required=False, help="Limit of entries to process in Mysql")
parser.add_argument("-user", "-u", type = str, choices=["setd", "jlnr"], required=True, help="Specify user")



args = parser.parse_args()
test = args.test
outfile = args.out
user = args.user



def make_table(outfile):
	outfile = 'newone.csv'
	
#------------------------------------------------------------------
# Connect to specific DB
#------------------------------------------------------------------
	try:
		db = mdb.connect(host="192.38.13.9", user ,passwd="1234",db="testasp")
	except mdb.Error, e:
		sys.exit("Cannot connect to database")#"# ERROR %d: %s" % (e.args[0],e.args[1]))
	try:
		cursor = db.cursor()
	except mdb.Error, e:
		sys.exit("# ERROR %d: %s" % (e.args[0],e.args[1]))
#------------------------------------------------------------------
# Create tables
#------------------------------------------------------------------
	try:
		query = """create table testNidCov as (
				SELECT `blast_qseq_jg3`, `blast_sseq_jg3`,
				ROUND(((((LEAST((`blast_send`-`blast_sstart`), (`blast_qend`-`blast_qstart`)))+1)/GREATEST(`blast_qlen`,`blast_slen`))*100),2) AS 'pct_longCov',
				ROUND(((((LEAST((`blast_send`-`blast_sstart`), (`blast_qend`-`blast_qstart`)))+1)/LEAST(`blast_qlen`,`blast_slen`))*100),2) AS 'pct_shortCov',
				`blast_pident`, `blast_evalue`, `blast_qstart`, `blast_qend`, `blast_qlen`, `blast_sstart`, `blast_send`, `blast_slen` from `blast`
				where `blast_sseq_jg1` like "Aspnid1" and `blast_qseq_jg1` like "Aspnid1");"""
		cursor.execute("drop table IF EXISTS testNidCov ;")
		cursor.execute(query)
	except mdb.Error, e:
		sys.exit('Setting up initial table failed')

	try:
		query="SELECT COUNT(*) FROM testNidCov;"
		cursor.execute(query)
		counter = cursor.fetchone()
		cursor.fetchall()
		print "Created table consists of %r rows" % counter
	except mdb.Error, e:
			sys.exit('Fetching row count failed')

# This query needs to be rewritten, when mysql is running again

	try:
		query = """CREATE TABLE testNidBiHits AS (
				SELECT * FROM (SELECT a.blast_qseq_jg3 AS aq, a.blast_sseq_jg3 AS ass, b.blast_qseq_jg3 AS bq, b.blast_sseq_jg3 AS bs, a.pct_longCov, a.pct_shortCov
				FROM testNidCov AS a LEFT JOIN testNidCov AS b
				ON (a.blast_qseq_jg3=b.blast_sseq_jg3 AND a.blast_sseq_jg3=b.blast_qseq_jg3) 
				WHERE a.blast_qseq_jg3 != a.blast_sseq_jg3 
				) AS c  where (ass > bs or bq is NOT NULL));""" 
		cursor.execute("drop table IF EXISTS testNidBiHits ;")
		cursor.execute(query)
	except mdb.Error, e:
		sys.exit('Setting up testNidBiHits table failed')

	try:
		query = """CREATE TABLE counting AS
			SELECT (SELECT COUNT(*) FROM testNidBiHits WHERE pct_longCov > 80 AND pct_shortCov >80) AS double80,
			(SELECT COUNT(*) FROM testNidBiHits AS longCov80 WHERE pct_longCov > 80) as longCov80 ,
			(SELECT COUNT(*) FROM testNidBiHits WHERE pct_shortCov > 80) AS shortCov80;"""
		cursor.execute("drop table IF EXISTS counting ;")
		cursor.execute(query)
	except mdb.Error, e:
		sys.exit('Counting of events failed')

	result="Sorry no have"

	print "Evaluating differences in coverage"
	try:
		query = "CREATE TABLE testNidCovDif AS (SELECT pct_shortCov-pct_longCov as CovDif FROM testNidBiHits)"
		#cursor.execute("drop table IF EXISTS testNidCovDif ;")
		cursor.execute(query)
	
	except mdb.Error, e:
		sys.exit('Failed to create table with Coverage differences')

	try:
		query = "SELECT * FROM testNidCovDif"
		cursor.execute(query)
		result = cursor.fetchall()
	except:
		print "Cannot fetch data from testNidCovDif table"
	
# SELECT (SELECT COUNT(*) FROM testNidBiHits WHERE pct_longCov > 80) as new WORKING!!!
#	db_cur.execute("SELECT * FROM counting")
#	pt = from_db_cursor(db_cur)

#	print pt
	print result[:10]

	f = open(outfile,'wb')
	writer = csv.writer(f, dialect = 'excel')
	writer.writerows(result)
	f.close()

#def make_plot(outfile):
#	print "# INFO: running Rscript"
#	os.system("R CMD BATCH '--args %s' aspmine_analysis_lengths.R test.out " % outfile)

'''	if test:
		try:
			query="SELECT * FROM initial_screening;"
			cursor.execute(query)
			print "Structure of initial_screening table \n %s" %cursor.fetchone()
		except mdb.Error, e:
			sys.exit('Fetching first table entries failed')''' # Implement this in later versions

make_table()			


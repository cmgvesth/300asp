#!/usr/bin/python

import MySQLdb as mdb
import sys, getopt, argparse
import warnings
warnings.filterwarnings('ignore')

#------------------------------------------------------------------
# Get commandline arguments
#------------------------------------------------------------------
#print "#============================================"
parser = argparse.ArgumentParser(description='# Create MySQL database and tables for the aspmine')

parser.add_argument('-clean', "--clean", action="store_true", required = False, help='Forces remove tables before recreate')
parser.add_argument('-cleandb', "--cleandb", action="store_true", required = False, help='Forces remove database before recreate')

parser.add_argument('-tab',required = True, help='Table to create')
parser.add_argument("-dbname", required=False, default = "aspminedb", help="R|Database name")
args = vars(parser.parse_args())

dbname = args['dbname']

#------------------------------------------------------------------
# Open mysql connection or exit
#------------------------------------------------------------------
try:
    db = mdb.connect("localhost","asp","1234", dbname)

except mdb.Error, e:
    print "Error %d: %s" % (e.args[0],e.args[1])
    sys.exit(1)

cur = db.cursor()

if args['cleandb'] or args["clean"]:
	cur.execute("SET foreign_key_checks=0")

#------------------------------------------------------------------
# Create database unless exists
#------------------------------------------------------------------
if args['cleandb']:
	cur.execute("DROP DATABASE IF EXISTS " + dbname +" ;")
	
	
	cur.execute("CREATE DATABASE IF NOT EXISTS " + dbname +" ;")
	print "INFO: database was recreated. Run again without -cleandb to create tables."
	sys.exit()

#------------------------------------------------------------------
# Create table error handling function
# Generic module for creating a table
# Do not alter the engine as the tables require foreign key functionalities 
#------------------------------------------------------------------
def create_tab(tab_create, tab_cols, tab_keys):
 
    try:
        cursor = db.cursor()
        cursor.execute(tab_create + " ( " + tab_cols + ", " + tab_keys + " ) ENGINE=INNODB;")

    except mdb.Error, e:
        print "Error %d: %s" % (e.args[0],e.args[1])
        sys.exit(1)

    except mdb.IntegrityError, e:
        print "IntegrityError %d: %s" % (e.args[0],e.args[1])
        sys.exit(1)

    except mdb.Warning, e:
        print "Warning %d: %s" % (e.args[0],e.args[1])
        sys.exit(1)

#------------------------------------------------------------------
# Create organism table unless exists 
#------------------------------------------------------------------
if args['tab'] == "organism" or args['tab'] == "all":
	if args["clean"]:
		cur.execute("DROP TABLE IF EXISTS organism ;")
		
	tab_create = "CREATE TABLE IF NOT EXISTS organism"

	tab_cols =	" `org_id` int(11) NOT NULL AUTO_INCREMENT, \
				`name` varchar(100) CHARACTER SET latin1 NOT NULL"

	tab_keys = "PRIMARY KEY (`org_id`),  KEY `iname` (`name`)"

	create_tab(tab_create, tab_cols, tab_keys)

#------------------------------------------------------------------
# Create species table unless exists 
#------------------------------------------------------------------
if args['tab'] == "species" or args['tab'] == "all":
	if args["clean"]:
		cur.execute("DROP TABLE IF EXISTS species ;")

	tab_create = "CREATE TABLE IF NOT EXISTS species"

	tab_cols = "species_id INT NOT NULL AUTO_INCREMENT, \
				species_section INT NOT NULL, \
				species_clade INT NOT NULL, \
				species_desc VARCHAR(100) NOT NULL"

	tab_keys = "PRIMARY KEY (species_id)"

	create_tab(tab_create, tab_cols, tab_keys)

#------------------------------------------------------------------
# Create org_is_species table unless exists 
#------------------------------------------------------------------
if args['tab'] == "org_is_species" or args['tab'] == "all":
	if args["clean"]:
		cur.execute("DROP TABLE IF EXISTS org_is_species ;")
		
	tab_create = "CREATE TABLE IF NOT EXISTS org_is_species"

	tab_cols = "species_id INT NOT NULL, \
				org_id INT NOT NULL "

	tab_keys = "PRIMARY KEY (`species_id`,`org_id`), \
				KEY `org_id` (`org_id`),\
				CONSTRAINT `org_is_species_ibfk_1` FOREIGN KEY (`org_id`) REFERENCES `organism` (`org_id`) ON UPDATE CASCADE,\
				CONSTRAINT `org_is_species_ibfk_2` FOREIGN KEY (`species_id`) REFERENCES `species` (`species_id`) ON UPDATE CASCADE"

	create_tab(tab_create, tab_cols, tab_keys)

#------------------------------------------------------------------
# Create environment table unless exists 
#------------------------------------------------------------------
if args['tab'] == "environment" or args['tab'] == "all":
	if args["clean"]:
		cur.execute("DROP TABLE IF EXISTS environment ;")

	tab_create = "CREATE TABLE IF NOT EXISTS environment"

	tab_cols = "env_id INT NOT NULL AUTO_INCREMENT, \
				env_desc VARCHAR(100) NOT NULL, \
				env_class VARCHAR(100) NOT NULL"

	tab_keys = "PRIMARY KEY (env_id)"

	create_tab(tab_create, tab_cols, tab_keys)

#------------------------------------------------------------------
# Create org_has_environment table unless exists 
#------------------------------------------------------------------
if args['tab'] == "org_has_environment" or args['tab'] == "all":
	if args["clean"]:
		cur.execute("DROP TABLE IF EXISTS org_has_environment ;")

	tab_create = "CREATE TABLE IF NOT EXISTS org_has_environment"

	tab_cols = "env_id INT NOT NULL, \
				org_id INT NOT NULL "

	tab_keys = "PRIMARY KEY (env_id, org_id), \
				FOREIGN KEY (org_id) REFERENCES organism(org_id) ON UPDATE CASCADE, \
				FOREIGN KEY (env_id) REFERENCES environment(env_id) ON UPDATE CASCADE "

	create_tab(tab_create, tab_cols, tab_keys)

#------------------------------------------------------------------
# Create assembly table unless exists - scaffolds unmasked
#------------------------------------------------------------------
if args['tab'] == "assembly" or args['tab'] == "all":
	if args["clean"]:
		cur.execute("DROP TABLE IF EXISTS assembly ;")

	tab_create = "CREATE TABLE IF NOT EXISTS assembly"

	tab_cols = "assembly_jgi1 VARCHAR(100) NOT NULL, \
				assembly_jgi2 VARCHAR(100) NOT NULL, \
				assembly_jgi3 VARCHAR(100) NOT NULL, \
				org_id INT NOT NULL, \
				assembly_seq_name VARCHAR(100) NOT NULL, \
				assembly_datatype VARCHAR(100), \
				assembly_seq LONGBLOB"

	tab_keys = "PRIMARY KEY (`assembly_seq_name`,`org_id`),\
				KEY `org_id` (`org_id`),\
				CONSTRAINT `assembly_ibfk_1` FOREIGN KEY (`org_id`) REFERENCES `organism` (`org_id`) ON DELETE CASCADE"

	create_tab(tab_create, tab_cols, tab_keys)

#------------------------------------------------------------------
# Create cds table unless exists
#------------------------------------------------------------------
if args['tab'] == "cds" or args['tab'] == "all":
	if args["clean"]:
		cur.execute("DROP TABLE IF EXISTS cds ;")
		
	tab_create = "CREATE TABLE IF NOT EXISTS cds"

	tab_cols = "`cds_jgi1` VARCHAR(100) NOT NULL, \
				`cds_jgi2` varchar(100) NOT NULL,\
				`cds_jgi3` varchar(100) NOT NULL,\
				`org_id` INT NOT NULL, \
				`cds_seq_name` VARCHAR(100) NOT NULL, \
				`cds_datatype` VARCHAR(100), \
				`cds_transcript_id` VARCHAR(100), \
				`cds_seq` LONGBLOB"

	tab_keys = "PRIMARY KEY (`cds_seq_name`,`org_id`),\
				KEY `org_id` (`org_id`),\
				CONSTRAINT `cds_ibfk_1` FOREIGN KEY (`org_id`) REFERENCES `organism` (`org_id`) ON DELETE CASCADE"

	create_tab(tab_create, tab_cols, tab_keys)

	tab_create = "CREATE TABLE IF NOT EXISTS cds"

#------------------------------------------------------------------
# Create proteins table unless exists
#------------------------------------------------------------------
if args['tab'] == "proteins" or args['tab'] == "all":
	if args["clean"]:
		cur.execute("DROP TABLE IF EXISTS proteins ;")

	tab_create = "CREATE TABLE IF NOT EXISTS proteins"

	tab_cols = "`prot_jgi1` varchar(100) NOT NULL,\
				`prot_jgi2` varchar(100) NOT NULL,\
				`prot_jgi3` varchar(100) NOT NULL,\
				`org_id` int(11) NOT NULL,\
				`prot_seq_name` varchar(100) NOT NULL,\
				`prot_datatype` varchar(100) DEFAULT NULL,\
				`prot_protein_id` varchar(100) DEFAULT NULL,\
				`prot_seq` longblob"

	tab_keys = "PRIMARY KEY (`prot_seq_name`,`org_id`),\
				KEY `org_id` (`org_id`),\
				KEY `iprot_seq_name` (`prot_seq_name`), \
				CONSTRAINT `proteins_ibfk_1` FOREIGN KEY (`org_id`) REFERENCES `organism` (`org_id`) ON DELETE CASCADE"

	create_tab(tab_create, tab_cols, tab_keys)

	tab_create = "CREATE TABLE IF NOT EXISTS proteins"

#------------------------------------------------------------------
# Create transcripts table unless exists
#------------------------------------------------------------------
if args['tab'] == "transcripts" or args['tab'] == "all":
	if args["clean"]:
		cur.execute("DROP TABLE IF EXISTS transcripts ;")

	tab_create = "CREATE TABLE IF NOT EXISTS transcripts"

	tab_cols = "`trans_jgi1` varchar(100) NOT NULL,\
				`trans_jgi2` varchar(100) NOT NULL,\
				`trans_jgi3` varchar(100) NOT NULL,\
				`org_id` int(11) NOT NULL DEFAULT '0',\
				`trans_seq_name` varchar(100) NOT NULL,\
				`trans_datatype` varchar(100) DEFAULT NULL,\
				`trans_transcript_id` varchar(100) DEFAULT NULL,\
				`trans_seq` longblob"
				

	tab_keys = "PRIMARY KEY (`trans_seq_name`,`org_id`),\
				KEY `org_id` (`org_id`),\
				KEY `itrans_transcript_id` (`trans_transcript_id`), \
				CONSTRAINT `transcripts_ibfk_1` FOREIGN KEY (`org_id`) REFERENCES `organism` (`org_id`) ON DELETE CASCADE"

	create_tab(tab_create, tab_cols, tab_keys)

	tab_create = "CREATE TABLE IF NOT EXISTS transcripts"

#------------------------------------------------------------------
# Create annotation tables unless exists
#------------------------------------------------------------------
if args['tab'] == "gff" or args['tab'] == "all":
	if args["clean"]:
		cur.execute("DROP TABLE IF EXISTS annotation ;")
	
	tab_create = "CREATE TABLE IF NOT EXISTS gff"

	# "seq_id", "source", "type", "start", "end", "score", "strand", "phase", "attributes"

	tab_cols = "gff_name VARCHAR(50), \
				gff_seqorigin VARCHAR(100) NOT NULL, \
				org_id INT, \
				gff_type VARCHAR(100) NOT NULL, \
				gff_protein_id VARCHAR(100), \
				gff_transcript_id VARCHAR(100), \
				gff_start INT, \
				gff_end INT, \
				gff_length INT, \
				gff_score VARCHAR(10), \
				gff_datatype VARCHAR(100), \
				gff_strand VARCHAR(10), \
				gff_phase VARCHAR(10), \
				gff_attributes VARCHAR(500)"

	tab_keys = "PRIMARY KEY (org_id, gff_attributes, gff_start), FOREIGN KEY (org_id) REFERENCES organism(org_id) ON DELETE CASCADE"

	create_tab(tab_create, tab_cols, tab_keys)

#------------------------------------------------------------------
# Create blast tables unless exists
#------------------------------------------------------------------
if args['tab'] == "blast" or args['tab'] == "all":
	tab_create = "CREATE TABLE IF NOT EXISTS blast"

	# qseq_id sseq_id pident qlen qstart qend slen sstart send evalue bitscore
	if args["clean"]:
		cur.execute("DROP TABLE IF EXISTS blast")

	tab_cols = "qseq_id VARCHAR(100), \
				sseq_id VARCHAR(100), \
				pident VARCHAR(50), \
				qlen INT, \
				qstart INT, \
				qend INT, \
				slen INT, \
				sstart INT, \
				send INT, \
				bitscore INT, \
				evalue VARCHAR(10)"
				
	tab_keys = "PRIMARY KEY (bitscore,qseq_id,sseq_id)"

	create_tab(tab_create, tab_cols, tab_keys)
	
#------------------------------------------------------------------
# Create blast tables unless exists
#------------------------------------------------------------------
if args['tab'] == "blast2org" or args['tab'] == "all":
	tab_create = "CREATE TABLE IF NOT EXISTS blast2org"

	# qseq_id sseq_id pident qlen qstart qend slen sstart send evalue bitscore
	if args["clean"]:
		cur.execute("DROP TABLE IF EXISTS blast2org")

	tab_cols = "qseq_id VARCHAR(50), \
				sseq_id VARCHAR(50), \
				qorg_id INT, \
				sorg_id INT, \
				qprot_seq_name VARCHAR(100), \
				sprot_seq_name VARCHAR(100)"

	tab_keys = "PRIMARY KEY (qseq_id,sseq_id), FOREIGN KEY (qorg_id) REFERENCES organism(org_id) ON DELETE CASCADE, FOREIGN KEY (sorg_id) REFERENCES organism(org_id) ON DELETE CASCADE"

	create_tab(tab_create, tab_cols, tab_keys)	

#------------------------------------------------------------------
# Create KEGG
# protein_id, ecNum, definition, catalyticActivity, cofactors, associatedDiseases, pathway, pathway_class, pathway_type, org_id, prot_seq_name
#------------------------------------------------------------------
if args['tab'] == "kegg" or args['tab'] == "all":
	tab_create = "CREATE TABLE IF NOT EXISTS kegg"

	#transcript_id	protein_id	kogid	kogdefline	kogClass	kogGroup

	if args["clean"]:
		cur.execute("DROP TABLE IF EXISTS kegg")

	tab_cols = "kegg_id INT AUTO_INCREMENT, \
				kegg_protein_id VARCHAR(50), \
				kegg_ecNum VARCHAR(50), \
				kegg_definition VARCHAR(500), \
				kegg_catalyticActivity VARCHAR(500), \
				kegg_cofactors VARCHAR(500), \
				kegg_associatedDiseases VARCHAR(500), \
				kegg_pathway VARCHAR(500), \
				kegg_pathway_class VARCHAR(500), \
				kegg_pathway_type VARCHAR(500), \
				org_id INT, \
				kegg_prot_seq_name VARCHAR(100)"

	tab_keys = "KEY (kegg_id), PRIMARY KEY (kegg_protein_id,kegg_id,org_id), FOREIGN KEY (org_id) REFERENCES organism(org_id) ON DELETE CASCADE"

	create_tab(tab_create, tab_cols, tab_keys)
	
#------------------------------------------------------------------
# Create KOG
#------------------------------------------------------------------
if args['tab'] == "kog" or args['tab'] == "all":
	tab_create = "CREATE TABLE IF NOT EXISTS kog"

	#transcript_id	protein_id	kogid	kogdefline	kogClass	kogGroup

	if args["clean"]:
		cur.execute("DROP TABLE IF EXISTS kog")

	tab_cols = "kog_transcript_id VARCHAR(50), \
				kog_protein_id VARCHAR(50), \
				kog_id VARCHAR(50), \
				kog_defline VARCHAR(500), \
				kog_Class VARCHAR(500), \
				kog_Group VARCHAR(500), \
				org_id INT, \
				kog_prot_seq_name VARCHAR(100)"

	tab_keys = "PRIMARY KEY (kog_transcript_id,kog_id,org_id), FOREIGN KEY (org_id) REFERENCES organism(org_id) ON DELETE CASCADE"

	create_tab(tab_create, tab_cols, tab_keys)
					
#------------------------------------------------------------------
# Create SignalP
#------------------------------------------------------------------
if args['tab'] == "sigp" or args['tab'] == "all":
	tab_create = "CREATE TABLE IF NOT EXISTS sigp"

	#protein_id	nn_cutpos	neuro_net_vote	hmm_cutpos	hmm_signalpep_probability

	if args["clean"]:
		cur.execute("DROP TABLE IF EXISTS sigp")

	tab_cols = "sigp_protein_id VARCHAR(50), \
				sigp_nn_cutpos INT, \
				sigp_neuro_net_vote INT, \
				sigp_hmm_cutpos INT, \
				sigp_hmm_signalpep_probability FLOAT,\
				org_id INT, \
				sigp_prot_seq_name VARCHAR(100)"

	tab_keys = "PRIMARY KEY (sigp_protein_id,sigp_nn_cutpos,org_id), FOREIGN KEY (org_id) REFERENCES organism(org_id) ON DELETE CASCADE"

	create_tab(tab_create, tab_cols, tab_keys)
	
#------------------------------------------------------------------
# Create GO
#------------------------------------------------------------------
if args['tab'] == "go" or args['tab'] == "all":
	tab_create = "CREATE TABLE IF NOT EXISTS go"

	#protein_id	nn_cutpos	neuro_net_vote	hmm_cutpos	hmm_signalpep_probability

	if args["clean"]:
		cur.execute("DROP TABLE IF EXISTS go")

	tab_cols = "go_protein_id VARCHAR(100), \
				go_term_id TEXT, \
				go_name TEXT, \
				go_termtype TEXT, \
				go_acc TEXT,\
				org_id INT, \
				go_prot_seq_name VARCHAR(100)"

	tab_keys = "PRIMARY KEY (go_protein_id,go_acc,org_id), FOREIGN KEY (org_id) REFERENCES organism(org_id) ON DELETE CASCADE"

	create_tab(tab_create, tab_cols, tab_keys)
#protein_id, goterm_id, goname, gotermType, goacc, org_id, prot_seq_name


#------------------------------------------------------------------
# Create InterPro
#protein_id, ipr_id, iprDesc, domainDb, domain_id, domainDesc, numHits, domainStarts, domainEnds, score
#------------------------------------------------------------------
if args['tab'] == "ipr" or args['tab'] == "all":
	tab_create = "CREATE TABLE IF NOT EXISTS ipr"

	#protein_id	nn_cutpos	neuro_net_vote	hmm_cutpos	hmm_signalpep_probability

	if args["clean"]:
		cur.execute("DROP TABLE IF EXISTS ipr")

	tab_cols = "ipr_protein_id VARCHAR(50), \
				ipr_id VARCHAR(50), \
				ipr_desc TEXT, \
				ipr_domaindb VARCHAR(50), \
				ipr_domain_id VARCHAR(100),\
				ipr_domaindesc TEXT,\
				ipr_numHits INT,\
				ipr_domain_starts VARCHAR(100),\
				ipr_domain_ends VARCHAR(100),\
				ipr_score TINYTEXT,\
				org_id INT, \
				ipr_prot_seq_name VARCHAR(100)"

	tab_keys = "PRIMARY KEY (ipr_protein_id,ipr_id,org_id,ipr_domain_id,ipr_domain_starts), FOREIGN KEY (org_id) REFERENCES organism(org_id) ON DELETE CASCADE"

	create_tab(tab_create, tab_cols, tab_keys)

#------------------------------------------------------------------
# CLOSE database
#------------------------------------------------------------------
cur.execute("SET foreign_key_checks=1")
db.commit()	
db.close()



#!/usr/bin/python

#------------------------------------------------------------------
# Imports
#------------------------------------------------------------------
import MySQLdb as mdb
import sys, getopt, argparse, re, glob, os, gzip
from datetime import datetime
sys.path.append('/home/tcve/github/db_mysql/')
from utilsFreedb import * # custom functions

##############################################################################################################################################
##############################################################################################################################################
####################################################################### ARGUMENTS ############################################################
##############################################################################################################################################
##############################################################################################################################################

db = connect_testdb()
cursor = db.cursor()

sql = "CREATE TABLE rep AS SELECT a.*, b.* from blast AS a JOIN blast AS b USING (a.blast_qseq_jg1=b.blast_sseq_jg1 and b.blast_qseq_jg1=a.blast_sseq_jg1 and a.blast_qseq_jg1<b.blast_sseq_jg1)"
print type(db)



"""

SELECT 
a.blast_qseq_jg2, 
a.blast_sseq_jg2,
a.blast_pident, 
b.blast_qseq_jg2, 
b.blast_sseq_jg2,
b.blast_pident
from blast AS a JOIN blast AS b ON (a.blast_qseq_jg1=b.blast_qseq_jg1) limit 10




"""
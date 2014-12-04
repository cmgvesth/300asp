#!/usr/bin/python

#------------------------------------------------------------------
# Imports
#------------------------------------------------------------------
import MySQLdb as mdb
import sys, getopt, argparse, re, glob, os, gzip
from datetime import datetime

##############################################################################################################################################
##############################################################################################################################################
####################################################################### ARGUMENTS ############################################################
##############################################################################################################################################
##############################################################################################################################################

def connect_testdb ():
	#startTime = datetime.now() # record runtime

	#------------------------------------------------------------------
	# Connect to specific DB
	#------------------------------------------------------------------
	try:
		host = "sql2.freemysqlhosting.net"	
		user = "sql260120"	
		password = "kE7!uV6*"
		database ="sql260120"
		db = mdb.connect(host, user, password, database)
		print "# INFO: success connectiong to existing database %s" % database

	except mdb.Error, e:
		print "# INFO: could not connect to db %s" % database
		sys.exit("# ERROR %d: %s" % (e.args[0],e.args[1]))


	cursor = db.cursor()
	sql = "SELECT * from blast limit 1"
	cursor.execute(sql)
	results = str(cursor.fetchone())

	test = "('Aspni_DSM_1', 158691L, 'An04g00060m.01', 'Aspni_DSM_1', 164093L, 'An11g08840m.01', Decimal('66.98'), 233L, 23L, 233L, 220L, 18L, 220L, 298L, Decimal('0.0'), 'jgi|Aspni_DSM_1|164093|An11g08840m.01', 'jgi|Aspni_DSM_1|158691|An04g00060m.01', 'Anigercbs51388VsAnigercbs51388Table.txt')"

	"""
	print test
	print results
	print type(test), type(results)
	print cmp(results, test)
	print results == test
	"""

	if results == test:
		print "# INFO testdb function test: query returned correct output"
	if results != test: 
		sys.exit("# ERROR query did not return the expected result")

	return db	
#!/usr/bin/python

#------------------------------------------------------------------
# Imports
#------------------------------------------------------------------
'''##################################################################
# Imports
##################################################################'''
import sys 
from aspmine_imports import *

#######################################################################
# Execute query or error
# Returns both column names and result: 	
# (columns, result)	 = executeQuery(cursor, tquery)
#######################################################################
def executeQuery(cursor, query):
	try:
		cursor.execute(query)
		columns = map(lambda x:x[0], cursor.description) 	
		result = cursor.fetchall()

	except mdb.Error, e:
		sys.exit("# ERROR %d: %s" % (e.args[0],e.args[1]))
	return columns, result	

#######################################################################
# Connect to db or error
# Returns cursor	
#######################################################################
def DBconnect(host, user, passw, dbname):
	try:
		db = mdb.connect( host, user, passw, dbname )
		cursor = db.cursor()
	except mdb.Error, e:
		sys.exit( "# ERROR %d: %s" % ( e.args[0],e.args[1] ) )
	return cursor	

#######################################################################
# Write mysql resutl to CSV
#######################################################################
def cursor2csv(columns, result, filename):
	f = open(filename,'wb')
	writer = csv.writer(f, dialect = 'excel', delimiter=";")
	writer.writerow(columns)
	writer.writerows(result)
	f.close()
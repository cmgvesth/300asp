def asp_con(path, user, pw):
	import MySQLdb as mdb
#------------------------------------------------------------------
# Connect to specific DB
#------------------------------------------------------------------
	try:
		db = mdb.connect(host=path, user=user ,passwd=pw,db= 'aspminedb') #host="localhost", user="root" ,passwd="Gru3n3r+T33",db= database)
	except mdb.Error, e:
		sys.exit("Cannot connect to database")#"# ERROR %d: %s" % (e.args[0],e.args[1]))
	try:
		cursor = db.cursor()
	except mdb.Error, e:
		sys.exit("# ERROR %d: %s" % (e.args[0],e.args[1]))
	return cursor
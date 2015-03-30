import sys,os
sys.path.append("../utils/")
from aspmine_imports import *
import pandas as pd


template=





f = open(outfile,'wb')

writer = csv.writer(f, dialect = 'excel')
writer.writerows(result)
f.close()


def to_txt():
	SELECT orderNumber, status, orderDate, requiredDate, comments
	FROM orders
	WHERE status = 'Cancelled'
	INTO OUTFILE 'C:/tmp/cancelled_orders.csv'
	FIELDS ENCLOSED BY '"' TERMINATED BY ';' ESCAPED BY '"'
	LINES TERMINATED BY '\r\n';
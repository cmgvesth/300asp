#!/usr/bin/python

#------------------------------------------------------------------
# Imports
#------------------------------------------------------------------
import MySQLdb as mdb
import getopt, argparse, re, glob, os, gzip
from datetime import datetime
import time
import xlrd, csv

''' bio python '''
from Bio import SeqIO
from Bio import SearchIO

from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

''' local libs '''
#sys.path.append('../utils/')
from utilsArgparse import * # custom functions
from utilsDataFormats import *
from utilsFileOperations import *
from utilsDBoperations import *
from utilsLogger import *
from utilsMySQLload import *

import warnings
warnings.filterwarnings("ignore", "Unknown table.*")
warnings.filterwarnings("ignore", "Data truncated.*")

sys.path.append(os.path.join(os.path.dirname(__file__), "../utils/util_snip"))
sys.path.append(os.path.join(os.path.dirname(__file__), "../utils/bio_snip"))


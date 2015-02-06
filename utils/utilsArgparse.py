#!/usr/bin/python
#######################################################################
"""
Costum CustomArgumentParser and SmartFormatter

Classes for custom argument parser.
To include these in code import using the following:
from utilsArgparse import * # custom functions

To use the custom parser create argument parser as follows:
parser = CustomArgumentParser(
	formatter_class=SmartFormatter, 
	description="Load data into local mysql database", 
	usage='%(prog)s -filetype filetype -action action -source filename\n')

To make use of the SmartFormatter add a "R|" to the beginning of the help text for each argument:
parser.add_argument('-action' , "-a", required=True, default = "test", choices=['load', 'test'],
					help='R|Test or test and load data\n'
					"-test\t: Read files and verify formats\n"
					"-load\t: Read files and verify formats followed by storing data in MySQL tables\n")

"""
#######################################################################

'''##################################################################
# Imports
##################################################################'''
import getopt, argparse, re, glob, os, gzip, sys

#from aspmine_imports import *

#######################################################################
# Format argument helper
#######################################################################
class SmartFormatter(argparse.HelpFormatter):
	width = 100
	def _split_lines(self, text, width):
		if text.startswith('R|'):
			return text[2:].splitlines()  
		return argparse.HelpFormatter._split_lines(self, text, width)

#######################################################################
# Format argument parser, allows help display when no arguments are given
#######################################################################
class CustomArgumentParser(argparse.ArgumentParser):
    def error(self, message):
		print "#--------------------------------------------------------------"
		print "# HELP:"
		print "#--------------------------------------------------------------"
		self.print_help()
		print "#--------------------------------------------------------------"
		print "# ERROR: %s" % message
		print "#--------------------------------------------------------------"
		sys.exit()
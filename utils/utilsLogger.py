#!/usr/bin/python
#######################################################################
import getopt, argparse, re, glob, os, gzip, sys

class Logger(object):
    def __init__(self, filename):
        self.terminal = sys.stdout
        self.log = open(filename, "w")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)  
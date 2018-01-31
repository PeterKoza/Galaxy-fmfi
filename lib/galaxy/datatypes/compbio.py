import logging, os, sys, time, sets, tempfile, shutil
import data
from galaxy import util
from galaxy.datatypes.sniff import *
from cgi import escape
import urllib
from bx.intervals.io import *
from galaxy.datatypes.metadata import MetadataElement


class Descriptor(data.Text):
    """delimited data in descriptor format"""
    file_ext = "des"

    MetadataElement( name="deskriptor", default=0, desc="RNA motif deskriptor", readonly=True )

    def sniff( self, filename ):
        print "\n" * 10
        print "here"
        print "\n" * 10
        has_header = False
        self.header = "" 
        try:
            fh = open(filename)
            while True:
                line = fh.readline().strip()
                if line[0] == "#":
                    pass
                if has_header == False:
                    self.header = line.split(" ")
                if line[0:2] not in self.header:
                    return False
                if not line:
                    break  # EOF
            fh.close()
        except Exception:
            pass
        return True
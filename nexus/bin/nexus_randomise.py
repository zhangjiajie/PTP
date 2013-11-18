#!/usr/bin/env python
import sys
import os

from nexus import NexusReader, NexusWriter, VERSION
from nexus.tools import shufflenexus

__author__ = 'Simon Greenhill <simon@simon.net.nz>'
__doc__ = """randomise - python-nexus tools v%(version)s
Shuffles the characters between each taxon to create a new nexus
""" % {'version': VERSION,}

if __name__ == '__main__':
    #set up command-line options
    from optparse import OptionParser
    parser = OptionParser(usage="usage: %prog -n 100 fudge.nex output.nex")
    parser.add_option("-n", "--numchars", dest="numchars", 
            action="store", default=False, 
            help="Number of Characters to Generate")
    options, args = parser.parse_args()
    
    try:
        nexusname = args[0]
    except IndexError:
        print __doc__
        print "Author: %s\n" % __author__
        parser.print_help()
        sys.exit()
        
    try:
        newnexus = args[1]
    except IndexError:
        newnexus = None
    
    if options.numchars != False:
        try:
            options.numchars = int(options.numchars)
        except ValueError:
            print "numchars needs to be a number!"
            raise
        
    nexus = NexusReader(nexusname)
    nexus = shufflenexus(nexus, options.numchars)
    if newnexus is not None:
        nexus.write_to_file(newnexus)
        print "New random nexus written to %s" % newnexus
    else:
        print nexus.write()    
        
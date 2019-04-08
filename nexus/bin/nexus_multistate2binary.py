#!/usr/bin/env python
import sys
import os
from nexus import NexusReader, NexusWriter, NexusFormatException, VERSION
from nexus.tools import binarise
__author__ = 'Simon Greenhill <simon@simon.net.nz>'
__doc__ = """nexus_multistate2binary - python-nexus tools v%(version)s

Converts multistate nexuses to binary present/absent form.
""" % {'version': VERSION,}

if __name__ == '__main__':
    #set up command-line options
    from optparse import OptionParser
    parser = OptionParser(usage="usage: %prog old.nex new.nex")
    parser.add_option("-1", "--onefile", dest="onefile", 
            action="store_true", default=False, 
            help="One nexus file for each multistate character")
    options, args = parser.parse_args()
    
    try:
        nexusname = args[0]
        newnexusname = args[1]
    except IndexError:
        print(__doc__)
        print("Author: %s\n" % __author__)
        parser.print_help()
        sys.exit()
        
    nexus = NexusReader(nexusname)
    
    new = binarise(nexus, one_nexus_per_block=options.onefile)
    if isinstance(new, NexusWriter):
        new.write_to_file(newnexusname)
    elif len(new) > 1:
        newnexusname, ext = os.path.splitext(newnexusname)
        for nex in new:
            nex.write_to_file("%s-%s%s" % (newnexusname, nex.clean(nex.characters[0]), ext))
            
    

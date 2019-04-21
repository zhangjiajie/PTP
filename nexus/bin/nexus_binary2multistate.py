#!/usr/bin/env python
import sys
from nexus import NexusReader, VERSION
from nexus.tools import multistatise, combine_nexuses

__author__ = 'Simon Greenhill <simon@simon.net.nz>'
__doc__ = """nexus_binary2multistate - python-nexus tools v%(version)s

Converts binary nexuses to a multistate nexus.
""" % {'version': VERSION,}

if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser(usage="usage: %prog [-o output.nex] nex1.nex nex2.nex ... nexN.nex")
    parser.add_option("-o", "--output", dest="output", 
            action="store", default=None, type="string",
            help="output nexus file")
    options, nexuslist = parser.parse_args()
    
    if len(nexuslist) < 1:
        print(__doc__)
        parser.print_help()
        sys.exit()
        
    if options.output is not None:
        outfile = options.output
    else:
        outfile = 'multistate.nex'
    
    nexuslist2 = []
    for nfile in nexuslist:
        n = NexusReader(nfile)
        n = multistatise(n)
        nexuslist2.append(n)
        
    out = combine_nexuses(nexuslist2)
    
    out.write_to_file(outfile, charblock=True, interleave=False)
    print("Written to %s" % outfile)

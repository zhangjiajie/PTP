#!/usr/bin/env python
import sys
from nexus import NexusReader, VERSION
from nexus.tools import combine_nexuses

__author__ = 'Simon Greenhill <simon@simon.net.nz>'
__doc__ = """combine-nexus - python-nexus tools v%(version)s
combines a series of nexuses into one nexus.
""" % {'version': VERSION,}


if __name__ == '__main__':
    #set up command-line options
    from optparse import OptionParser
    parser = OptionParser(usage="usage: %prog nex1.nex nex2.nex ... nexN.nex")
    options, nexuslist = parser.parse_args()
    
    if len(nexuslist) <= 1:
        print(__doc__)
        parser.print_help()
        sys.exit()
    
    nexuslist = [NexusReader(n) for n in nexuslist]
    out = combine_nexuses(nexuslist)
    out.write_to_file('combined.nex', charblock=False, interleave=False)
    print("Written to combined.nex")

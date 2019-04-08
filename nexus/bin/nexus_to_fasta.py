#!/usr/bin/env python
import sys
from nexus import NexusReader, VERSION
from textwrap import wrap

__author__ = 'Simon Greenhill <simon@simon.net.nz>'
__doc__ = """nexus_to_fasta - python-nexus tools v%(version)s
converts a nexus file to a fasta file.
""" % {'version': VERSION,}

if __name__ == '__main__':
    #set up command-line options
    from optparse import OptionParser
    parser = OptionParser(usage="usage: %prog nexus")
    options, args = parser.parse_args()
    
    try:
        nexusname = args[0]
    except IndexError:
        parser.print_help()
        sys.exit()
        
    n = NexusReader(nexusname)
    for taxon in sorted(n.data.matrix):
        print('>%s' % taxon)
        for line in wrap("".join(n.data.matrix[taxon]), 70):
            print(line)
    

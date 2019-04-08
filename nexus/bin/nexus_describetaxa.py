#!/usr/bin/env python
"""Displays some simple nexus info"""
__author__ = 'Simon Greenhill <simon@simon.net.nz>'

import sys
from nexus import NexusReader
from nexus.tools import check_for_valid_NexusReader


def print_taxa_stats(nexus_obj):
    """
    Prints the taxa state statistics for a given `nexus_obj`
    
    :param nexus_obj: A `NexusReader` instance
    :type nexus_obj: NexusReader 
    
    :raises AssertionError: if nexus_obj is not a nexus
    :raises NexusFormatException: if nexus_obj does not have a `data` block
    """
    check_for_valid_NexusReader(nexus_obj, required_blocks=['data'])
    
    for taxon in sorted(nexus_obj.data.matrix):
        tally = {}
        for site in nexus_obj.data.matrix[taxon]:
            tally[site] = tally.get(site, 0) + 1
            
        tally = ", ".join(['%s x %s' % (k,tally[k]) for k in sorted(tally)])
        print(taxon.ljust(20), tally)
    return

if __name__ == '__main__':
    #set up command-line options
    from optparse import OptionParser
    parser = OptionParser(usage="usage: %prog nexusname")
    options, args = parser.parse_args()
    
    try:
        nexusname = args[0]
    except IndexError:
        parser.print_help()
        sys.exit()
    
    n = NexusReader(nexusname)
    print(n)
    for k, v in n.blocks.items():
        print(' ', k, v)
        
        if k == 'data':
            print_taxa_stats(n)

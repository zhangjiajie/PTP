#!/usr/bin/env python
from __future__ import division
import sys
import os
from nexus import NexusReader, NexusWriter, NexusFormatException, VERSION
from nexus.tools import count_site_values, \
    find_constant_sites, find_unique_sites, new_nexus_without_sites

__author__ = 'Simon Greenhill <simon@simon.net.nz>'
__doc__ = """nexusmanip - python-nexus tools v%(version)s

Performs a number of nexus manipulation methods.
""" % {'version': VERSION,}


def print_site_values(nexus_obj, characters=['-', '?']):
    """
    Prints out counts of the number of sites with state in `characters` in a nexus.

    (Wrapper around `count_site_values`)

    :param nexus_obj: A `NexusReader` instance
    :type nexus_obj: NexusReader 
    """
    count = count_site_values(nexus_obj, characters)
    print ("Number of %s in %s" % (",".join(characters), nexus_obj.filename))
    for taxon in sorted(count):
        print("%s: %d/%d (%0.2f%%)" % \
            (taxon.ljust(20), count[taxon], nexus.data.nchar, (count[taxon]/nexus.data.nchar)*100)
        )
    print('-'*76)
    total_count = sum([x for x in count.values()])
    total_data = nexus.data.nchar * nexus.data.ntaxa
    print('TOTAL: %d/%d (%0.2f%%)' % \
        (total_count, total_data, (total_count/total_data)*100)
    )

def print_character_stats(nexus_obj):
    """
    Prints the number of states and members for each site in `nexus_obj`
    
    :param nexus_obj: A `NexusReader` instance
    :type nexus_obj: NexusReader 
    
    :return: A list of the state distribution
    """
    state_distrib = []
    for i in range(0, nexus_obj.data.nchar):
        tally = {}
        for taxa, characters in nexus_obj.data:
            c = characters[i]
            tally[c] = tally.get(c, 0) + 1
        
        print("%5d" % i)
        for state in tally:
            print("%sx%d" % (state, tally[state]),
            state_distrib.append(tally[state]))
        
    return state_distrib
    



if __name__ == '__main__':
    #set up command-line options
    from optparse import OptionParser
    parser = OptionParser(usage="usage: %prog old.nex [new.nex]")
    parser.add_option("-n", "--number", dest="number", 
            action="store_true", default=False, 
            help="Count the number of characters")
    parser.add_option("-c", "--constant", dest="constant", 
            action="store_true", default=False, 
            help="Remove the constant characters")
    parser.add_option("-u", "--unique", dest="unique", 
            action="store_true", default=False, 
            help="Remove the unique characters")
    parser.add_option("-s", "--stats", dest="stats", 
            action="store_true", default=False, 
            help="Print character-by-character stats")
    options, args = parser.parse_args()
    
    try:
        nexusname = args[0]
    except IndexError:
        print(__doc__)
        print("Author: %s\n" % __author__)
        parser.print_help()
        sys.exit()
        
    try:
        newnexusname = args[1]
    except IndexError:
        newnexusname = None
        
    
    nexus = NexusReader(nexusname)
    newnexus = None
    
    if options.number:
        if newnexusname is not None:
            print_site_values(nexus, newnexusname)
        else:
            print_site_values(nexus)
            
    elif options.constant:
        const = find_constant_sites(nexus)
        newnexus = new_nexus_without_sites(nexus, const)
        print("Constant Sites: %s" % ",".join([str(i) for i in const]))
    elif options.unique:
        unique = find_unique_sites(nexus)
        newnexus = new_nexus_without_sites(nexus, unique)
        print("Unique Sites: %s" % ",".join([str(i) for i in unique]))
    elif options.stats:
        d = print_character_stats(nexus)
        
    else:
        exit()
        
    # check for saving
    if newnexus is not None and newnexusname is not None:
        newnexus.write_to_file(newnexusname)
        print("New nexus written to %s" % newnexusname)
        
    

#!/usr/bin/env python
import os

from textwrap import TextWrapper

from nexus import NexusReader, NexusWriter, NexusFormatException, VERSION
from nexus.tools import tally_by_taxon, tally_by_site

__author__ = 'Simon Greenhill <simon@simon.net.nz>'
__doc__ = """nexusmanip - python-nexus tools v%(version)s

Performs a number of nexus counting/tallying methods.
""" % {'version': VERSION,}

def print_tally(tally): 
    wrapper = TextWrapper(initial_indent=" ", subsequent_indent="\t", width=65)
    for tkey in sorted(tally):
        print(tkey)
        for skey in sorted(tally[tkey]):
            s = " ".join(sorted(tally[tkey][skey]))
            print(" - %s: " % skey) 
            for w in wrapper.wrap(s):
                print(w)
    return
    

if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser(usage="usage: %prog [taxa/sites] nexus.nex")
    options, commands = parser.parse_args()
    
    if len(commands) != 2:
        print(__doc__)
        parser.print_help()
        quit()
    
    command, nex = commands
    
    try:
        nex = NexusReader(nex)
    except IOError:
        raise IOError("Unable to read %s" % nex)
        
    if command in ('taxa', 't'):
        tally = tally_by_taxon(nex)
    elif command in ('site', 's'):
        tally = tally_by_site(nex)
    else:
        quit("Invalid tally command. Only 'taxa' and 'site' are valid.")
    
    print_tally(tally)

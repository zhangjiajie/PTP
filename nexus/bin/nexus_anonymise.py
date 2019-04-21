#!/usr/bin/env python
import sys
import hashlib
from nexus import NexusReader, VERSION

__author__ = 'Simon Greenhill <simon@simon.net.nz>'
__doc__ = """nexus_anonymise - python-nexus tools v%(version)s

Anonymises the taxa in a nexus
""" % {'version': VERSION,}

def anonymise(nexus_obj):
    """Anonymises a nexus object"""
    for block in nexus_obj.blocks:
        if block == 'taxa':
            for idx, t in enumerate(nexus_obj.blocks[block].taxa):
                nexus_obj.blocks[block].taxa[idx] = hash(nexus_obj.filename, t)
        elif block == 'trees':
            if nexus_obj.blocks[block].was_translated:
                for idx in nexus_obj.blocks[block].translators:
                    nexus_obj.blocks[block].translators[idx] = hash(nexus_obj.filename,
                                                            nexus_obj.blocks[block].translators[idx])
            else:
                raise NotImplementedError("Unable to anonymise untranslated trees")
        elif block == 'data':
            newmatrix = {}
            for t in nexus_obj.blocks[block].matrix:
                newmatrix[hash(nexus_obj.filename, t)] = nexus_obj.blocks[block].matrix[t]
            nexus_obj.blocks[block].matrix = newmatrix
            
        else:
            raise NotImplementedError("Unable to anonymise %s blocks" % block)
    return nexus_obj


def hash(salt, taxon):
    return hashlib.md5("%s-%s" % (salt, taxon)).hexdigest()

if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser(usage="usage: %prog fudge.nex output.nex")
    options, nexuslist = parser.parse_args()
    
    try:
        nexusname = args[0]
    except IndexError:
        print(__doc__)
        print("Author: %s\n" % __author__)
        parser.print_help()
        sys.exit()
    
    try:
        newnexus = args[1]
    except IndexError:
        newnexus = None
        
    nexus = NexusReader(nexusname)
    nexus = anonymise(nexus)
    
    if newnexus is not None:
        nexus.write_to_file(newnexus)
        print("New nexus written to %s" % newnexus)
    else:
        print(nexus.write_to_file(hash('filename', filename)))

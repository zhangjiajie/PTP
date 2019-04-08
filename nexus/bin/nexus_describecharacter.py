#!/usr/bin/env python
import sys
from nexus import NexusReader, VERSION

__author__ = 'Simon Greenhill <simon@simon.net.nz>'
__doc__ = """describecharacter - python-nexus tools v%(version)s
Describes the given character.
""" % {'version': VERSION,}
import sys

from textwrap import TextWrapper

from nexus import NexusReader, NexusFormatException
from nexus.tools import check_for_valid_NexusReader

wrapper = TextWrapper(initial_indent="  ", subsequent_indent="  ")

def print_character_stats(nexus_obj, character_index):
    """
    Prints the character/site statistics for a given `nexus_obj` and character index
    
    :param nexus_obj: A `NexusReader` instance
    :type nexus_obj: NexusReader 
    
    :param character_index: The character index of the character to summarise
    :type character_index: Int or String 
    
    :raises AssertionError: if nexus_obj is not a nexus
    :raises IndexError: if character_index is not in nexus data block
    :raises NexusFormatException: if nexus_obj does not have a `data` block
    """
    check_for_valid_NexusReader(nexus_obj, required_blocks=['data'])
    
    index = None
    if character_index in nexus_obj.data.characters:
        index = character_index # string index
    else:
        try:
            character_index = int(character_index)
        except ValueError:
            pass
        
        if character_index in nexus_obj.data.characters:
            index = character_index
        
    if index is None:
        raise IndexError("Character '%s' is not in the nexus" % char)
    
    states = {}
    for taxon, state in nexus_obj.data.characters[index].items():
        states[state] = states.get(state, [])
        states[state].append(taxon)
    
    for state in sorted(states):
        print( 'State: %s (%d / %d = %0.2f)' % (state, 
            len(states[state]), nexus_obj.data.ntaxa, 
            (len(states[state]) / nexus_obj.data.ntaxa * 100)
        ))
        print("\n".join(wrapper.wrap(", ".join(states[state]))))
    
    return
    

if __name__ == '__main__':
    #set up command-line options
    from optparse import OptionParser
    parser = OptionParser(usage="usage: %prog site_index nexusfile.nex")
    options, args = parser.parse_args()
    
    try:
        char = args[0]
        nexusname = args[1]
    except IndexError:
        parser.print_help()
        sys.exit()
    
    print_character_stats(NexusReader(nexusname), char)

"""Contains Nexus Manipulation Tools that operate on Site/Characters"""
from nexus import NexusWriter
from nexus.tools import check_for_valid_NexusReader


def find_constant_sites(nexus_obj):
    """
    Returns a list of the constant sites in a nexus

    :param nexus_obj: A `NexusReader` instance
    :type nexus_obj: NexusReader 

    :return: A list of constant site positions.
    :raises AssertionError: if nexus_obj is not a nexus
    :raises NexusFormatException: if nexus_obj does not have a `data` block
    """
    check_for_valid_NexusReader(nexus_obj, required_blocks=['data'])
    
    const = []
    for i in range(0, nexus_obj.data.nchar):
        states = []
        for taxa, data in nexus_obj.data:
            c = data[i]
            if c in ('?', '-'):
                continue
            elif c not in states:
                states.append(c)

        if len(states) == 1:
            const.append(i)
    return const


def find_unique_sites(nexus_obj):
    """
    Returns a list of the unique sites in a binary nexus
    i.e. sites with only one taxon belonging to them.
        (this only really makes sense if the data is coded as presence/absence)

    :param nexus_obj: A `NexusReader` instance
    :type nexus_obj: NexusReader 

    :return: A list of unique site positions.
    :raises AssertionError: if nexus_obj is not a nexus
    :raises NexusFormatException: if nexus_obj does not have a `data` block
    """
    check_for_valid_NexusReader(nexus_obj, required_blocks=['data'])
    
    unique = []
    for i in range(0, nexus_obj.data.nchar):
        members = {}
        missing = 0
        for taxa, characters in nexus_obj.data:
            c = characters[i]
            if c in (u'?', u'-'):
                missing += 1
            else:
                members[c] = members.get(c, 0) + 1

        # a character is unique if there's only two states
        # AND there's a state with 1 member 
        # AND the state with 1 member is NOT the 0 (absence) state
        if len(members) == 2:
            for state, count in members.items():
                if state != '0' and count == 1:
                    unique.append(i)
    return unique



def count_site_values(nexus_obj, characters=('-', '?')):
    """
    Counts the number of sites with values in `characters` in a nexus

    :param nexus_obj: A `NexusReader` instance
    :type nexus_obj: NexusReader 

    :param characters: An iterable of the characters to count
    :type characters: tuple

    :return: (A dictionary of taxa and missing counts, and a list of log comments)
    :raises AssertionError: if nexus_obj is not a nexus
    :raises NexusFormatException: if nexus_obj does not have a `data` block
    """
    try:
        iter(characters) # just check it's iterable. Don't _change_ it to an
            # iterable or we won't be able to compare the characters properly below.
    except TypeError:
        raise TypeError("characters should be iterable")

    check_for_valid_NexusReader(nexus_obj, required_blocks=['data'])

    tally = {}
    for taxon, sites in nexus_obj.data:
        tally[taxon] = tally.get(taxon, 0)
        for site in sites:
            if site in characters:
                tally[taxon] += 1
    return tally
    
    

def new_nexus_without_sites(nexus_obj, sites_to_remove):
    """
    Returns a new NexusReader instance with the sites in 
    `sites_to_remove` removed.

    :param nexus_obj: A `NexusReader` instance
    :type nexus_obj: NexusReader 

    :param sites_to_remove: A list of site numbers
    :type sites_to_remove: List

    :return: A NexusWriter instance
    :raises AssertionError: if nexus_obj is not a nexus
    :raises NexusFormatException: if nexus_obj does not have a `data` block
    """
    check_for_valid_NexusReader(nexus_obj, required_blocks=['data'])

    # make new nexus
    nexout = NexusWriter()
    nexout.add_comment(
        "Removed %d sites: %s" %
        (len(sites_to_remove), ",".join(["%s" % s for s in sites_to_remove]))
    )
    new_sitepos = 0
    for sitepos in range(nexus_obj.data.nchar):
        if sitepos in sites_to_remove:
            continue # skip!
        for taxon, data in nexus_obj.data:
            nexout.add(taxon, new_sitepos, data[sitepos])
        new_sitepos += 1
    return nexout


def tally_by_site(nexus_obj):
    """
    Counts the number of taxa per state per site (i.e. site 1 has three taxa 
    coded as "A", and 1 taxa coded as "G")

    Returns a dictionary of the cognate sets by members in the nexus

    :param nexus_obj: A `NexusReader` instance
    :type nexus_obj: NexusReader 

    :return: A Dictionary
    :raises AssertionError: if nexus_obj is not a nexus
    :raises NexusFormatException: if nexus_obj does not have a `data` block

    e.g. {
        'site1': {'state1': ['taxon1', 'taxon2'], 'state0': ['taxon3'], } 
        'site2': {'state1': ['taxon2'], 'state0': ['taxon1', 'taxon3'], }
    }
    """
    check_for_valid_NexusReader(nexus_obj, required_blocks=['data'])
    tally = {}
    for site, data in nexus_obj.data.characters.items():
        tally[site] = tally.get(site, {})
        for taxon, state in data.items():
            tally[site][state] = tally[site].get(state, [])
            tally[site][state].append(taxon)
    return tally
    


def tally_by_taxon(nexus_obj):
    """
    Counts the number of states per site that each taxon has (i.e. taxon 1 
    has three sites coded as "A" and 1 coded as "G") 

    Returns a dictionary of the cognate sets by members in the nexus

    :param nexus_obj: A `NexusReader` instance
    :type nexus_obj: NexusReader 

    :return: A Dictionary
    :raises AssertionError: if nexus_obj is not a nexus
    :raises NexusFormatException: if nexus_obj does not have a `data` block

    e.g. {
        'taxon1': {'state1': ['site1', 'site2'], 'state0': ['site3'], } 
        'taxon2': {'state1': ['site2'], 'state0': ['site1', 'site3'], }
    }
    """
    check_for_valid_NexusReader(nexus_obj, required_blocks=['data'])
    tally = {}
    for taxon, characters in nexus_obj.data:
        tally[taxon] = {}
        for pos, char in enumerate(characters):
            label = nexus_obj.data.charlabels.get(pos, pos)
            tally[taxon][char] = tally[taxon].get(char, [])
            tally[taxon][char].append(label)
    return tally
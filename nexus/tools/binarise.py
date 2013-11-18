from nexus import NexusWriter
from nexus.tools import check_for_valid_NexusReader

def _recode_to_binary(char, keep_zero=False):
    """
    Recodes a dictionary to binary data. 
    
    :param char: A dictionary of taxa to state values
    :type char: dict 
    
    :param keep_zero: A boolean flag denoting whether to 
        treat '0' as a missing state or not. The default
        (False) is to ignore '0' as a trait absence.
        
        Setting this to True will treat '0' as a unique
        state.
    :type keep_zero: Boolean

    :return: A dictionary of taxa to recoded values.
    :raises ValueError: if any of the states in the 
        `char` dictionary is not a string (i.e. 
        integer or None values)
    
    >>> recode, table = _recode_to_binary({'Maori': '1', 'Dutch': '2', 'Latin': '1'})
    >>> recode['Maori']
    '10'
    >>> recode['Dutch']
    '01'
    >>> recode['Latin']
    '10'
    """
    
    newdata = {}
    
        # unwanted states
    unwanted_states = ['-', '?']
    if not keep_zero:
        unwanted_states.append('0')
    
    # get unique states
    states = list(set([c for c in char.values() if c not in unwanted_states])) 
    
    if not all(isinstance(s, basestring) for s in states):
        raise ValueError('Data must be strings')
    
    states = sorted(states)
    num_states = len(states)
    for taxon, value in char.items():
        newdata[taxon] = ['0' for x in range(num_states)]
        if value not in unwanted_states: # ignore missing values
            newdata[taxon][states.index(value)] = '1'
        newdata[taxon] = "".join(newdata[taxon])
        assert len(newdata[taxon]) == num_states
    
    return newdata

def binarise(nexus_obj, one_nexus_per_block=False):
    """
    Returns a binary variant of the given `nexus_obj`.
    If `one_nexus_per_block` then we return a list of NexusWriter instances.

    :param nexus_obj: A `NexusReader` instance
    :type nexus_obj: NexusReader 

    :param one_nexus_per_block: Whether to return a single NexusWriter, or a 
                                list of NexusWriter's (one per character)
    :type one_nexus_per_block: Boolean

    :return: A NexusWriter instance or a list of NexusWriter instances.
    :raises AssertionError: if nexus_obj is not a nexus
    :raises NexusFormatException: if nexus_obj does not have a `data` block
    """
    check_for_valid_NexusReader(nexus_obj, required_blocks=['data'])
    nexuslist = []
    n = NexusWriter()
    
    for i in sorted(nexus_obj.data.charlabels):
        label = nexus_obj.data.charlabels[i] # character label
        char = nexus_obj.data.characters[label] # character dict (taxon->state)
        recoding = _recode_to_binary(char) # recode
        
        new_char_length = len(recoding[recoding.keys()[0]])
        
        # loop over recoded data
        for j in range(new_char_length):
            for taxon, state in recoding.items():
                # make new label
                new_label = "%s_%d" % (str(label), j)
                # add to nexus
                n.add(taxon, new_label, state[j])
            
        if one_nexus_per_block:
            nexuslist.append(n)
            n = NexusWriter()
                
    if one_nexus_per_block:
        return nexuslist
    else:
        return n


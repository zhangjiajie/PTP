from nexus import NexusWriter
from nexus.tools import check_for_valid_NexusReader

def multistatise(nexus_obj):
    """
    Returns a multistate variant of the given `nexus_obj`.

    :param nexus_obj: A `NexusReader` instance
    :type nexus_obj: NexusReader 

    :return: A NexusReader instance
    :raises AssertionError: if nexus_obj is not a nexus
    :raises NexusFormatException: if nexus_obj does not have a `data` block
    """
    check_for_valid_NexusReader(nexus_obj, required_blocks=['data'])
        
    site_idx = 0
    nexout = NexusWriter()
    missing = []
    
    charlabel = getattr(nexus_obj, 'short_filename', 1)
    
    for site, data in nexus_obj.data.characters.items():
        multistate_value = chr(65 + site_idx)
        for taxon, value in data.items():
            assert value == str(value)
            if value in ('?', '-'):
                missing.append(taxon)
                
            if value == '1':
                nexout.add(taxon, charlabel, multistate_value)
                if taxon in missing: # remove taxon if we've seen a non-? entry
                    missing.remove(taxon)
        site_idx += 1
        assert site_idx < 26, "Too many characters to handle! - run out of A-Z"
        
    # add missing state for anything that is all missing, and has not been
    # observed anywhere
    for taxon in nexus_obj.data.taxa:
        if taxon not in nexout.data[str(charlabel)]:
            nexout.add(taxon, charlabel, '?')
    return nexout._convert_to_reader()
    

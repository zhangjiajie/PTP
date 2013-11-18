from nexus import NexusReader, NexusFormatException

def check_for_valid_NexusReader(nexus_obj, required_blocks=[]):
    """
    Performs some checking to make sure we received a valid NexusReader
    
    :param nexus_obj: A `NexusReader` instance
    :type nexus_obj: NexusReader 
    
    :param required_blocks: A list of nexus 'blocks' that are required
    :type one_nexus_per_block: List
    
    :return: Boolean True
    :raises TypeError: if nexus_obj is not a nexus
    :raises NexusFormatException: if nexus_obj does not have a `required_block`
    """
    if isinstance(nexus_obj, NexusReader) == False:
        raise TypeError("Nexus_obj should be a NexusReader instance")
    for b in required_blocks:
        if hasattr(nexus_obj, b) == False:
            raise NexusFormatException("Requires a `%s` block, but one was not found in nexus" % b)
    return True


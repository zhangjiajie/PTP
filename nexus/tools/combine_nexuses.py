import os
from nexus import NexusWriter
from nexus.tools import check_for_valid_NexusReader

def combine_nexuses(nexuslist):
    """
    Combines a list of NexusReader instances into a single nexus
    
    :param nexuslist: A list of NexusReader instances
    :type nexuslist: List 
    
    :return: A NexusWriter instance
    
    :raises TypeError: if nexuslist is not a list of NexusReader instances
    :raises IOError: if unable to read an file in nexuslist
    :raises NexusFormatException: if a nexus file does not have a `data` block
    """
    out = NexusWriter()
    charpos = 0
    for nex_id, nex in enumerate(nexuslist, 1):
        check_for_valid_NexusReader(nex, required_blocks=['data'])
        
        if hasattr(nex, 'short_filename'):
            nexus_label = os.path.splitext(nex.short_filename)[0]
        elif hasattr(nex, 'label'):
            nexus_label = nex.label
        else:
            nexus_label = str(nex_id)
        
        out.add_comment("%d - %d: %s" % (charpos, charpos + nex.data.nchar -1, nexus_label))
        for site_idx, site in enumerate(sorted(nex.data.characters), 0):
            data = nex.data.characters.get(site)
            charpos += 1
            # work out character label
            charlabel = nex.data.charlabels.get(site_idx, site_idx + 1)
            label = '%s.%s' % (nexus_label, charlabel)
            
            for taxon, value in data.items():
                out.add(taxon, label, value)
    return out
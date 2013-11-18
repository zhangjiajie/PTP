import os
import unittest

from nexus import NexusReader
from nexus.tools import find_unique_sites

EXAMPLE_DIR = os.path.join(os.path.dirname(__file__), '../../examples')

class Test_FindUniqueSites(unittest.TestCase):
    """Test find_unique_sites"""
    def test_find_unique_sites_1(self):
        nexus = NexusReader(os.path.join(EXAMPLE_DIR, 'example.nex'))
        assert len(find_unique_sites(nexus)) == 0

    def test_find_unique_sites_2(self):
        nexus = NexusReader()
        nexus.read_string("""Begin data;
        Dimensions ntax=4 nchar=7;
        Format datatype=standard symbols="01" gap=-;
        Matrix
        Harry              10000?-
        Simon              1100011
        Betty              1110000
        Louise             1111000
        ;""")
        unique = find_unique_sites(nexus)

        # site 1 should NOT be in the uniques (3x1 and 1x0)
        # - i.e. are we ignoring sites with ONE absent taxon
        assert 1 not in unique
        # these should also NOT be in unique
        assert 0 not in unique
        assert 2 not in unique
        assert 4 not in unique # constant
        # site 3 is a simple unique site - check we found it
        assert 3 in unique
        # sites 5 and 6 should also be unique 
        # - are we handling missing data appropriately?
        assert 5 in unique
        assert 6 in unique
        
if __name__ == '__main__':
    unittest.main()

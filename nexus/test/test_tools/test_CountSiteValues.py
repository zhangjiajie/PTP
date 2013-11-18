import os
import unittest

from nexus import NexusReader
from nexus.tools import count_site_values

EXAMPLE_DIR = os.path.join(os.path.dirname(__file__), '../../examples')

class Test_CountSiteValues(unittest.TestCase):
    def test_failure_on_notiterableargument(self):
        self.assertRaises(TypeError, count_site_values, 1)

    def test_count_missing_one(self):
        nexus = NexusReader(os.path.join(EXAMPLE_DIR, 'example.nex'))
        missing = count_site_values(nexus)
        for taxon in missing:
            assert missing[taxon] == 0

    def test_count_missing_two(self):
        expected = {'Harry': 0, 'Simon': 1, 'Peter': 1, 'Betty': 2, 'Louise': 3}
        nexus = NexusReader()
        nexus.read_string("""#NEXUS 
        Begin data;
        Dimensions ntax=5 nchar=3;
        Format datatype=standard symbols="01" gap=-;
        Matrix
        Harry              010  [No missing]
        Simon              0?0  [one missing]
        Peter              0-0  [one gap]
        Betty              ?-1  [one gap and one missing = 2 missing]
        Louise             ???  [three missing]
            ;
        End;
        """)
        missing = count_site_values(nexus)
        for taxon in missing:
            assert missing[taxon] == expected[taxon]

    def test_count_other_values_one(self):
        expected = {'Harry': 1, 'Simon': 1, 'Peter': 0, 'Betty': 0, 'Louise': 0}
        nexus = NexusReader()
        nexus.read_string("""#NEXUS 
        Begin data;
        Dimensions ntax=5 nchar=3;
        Format datatype=standard symbols="01" gap=-;
        Matrix
        Harry              0A0  [No missing]
        Simon              0A0  [one missing]
        Peter              0-0  [one gap]
        Betty              ?-1  [one gap and one missing = 2 missing]
        Louise             ???  [three missing]
            ;
        End;
        """)
        count = count_site_values(nexus, 'A')
        for taxon in count:
            assert count[taxon] == expected[taxon]

    def test_count_other_values_two(self):
        expected = {'Harry': 1, 'Simon': 2, 'Peter': 1, 'Betty': 0, 'Louise': 0}
        nexus = NexusReader()
        nexus.read_string("""#NEXUS 
        Begin data;
        Dimensions ntax=5 nchar=3;
        Format datatype=standard symbols="01" gap=-;
        Matrix
        Harry              0A0  [No missing]
        Simon              0AB  [one missing]
        Peter              0-B  [one gap]
        Betty              ?-1  [one gap and one missing = 2 missing]
        Louise             ???  [three missing]
            ;
        End;
        """)
        count = count_site_values(nexus, ['A', 'B'])
        for taxon in count:
            assert count[taxon] == expected[taxon]


if __name__ == '__main__':
    unittest.main()

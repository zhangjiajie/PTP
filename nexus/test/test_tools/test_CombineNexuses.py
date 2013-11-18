import os
import re
import unittest

from nexus import NexusReader
from nexus.tools import combine_nexuses

EXAMPLE_DIR = os.path.join(os.path.dirname(__file__), '../../examples')

class Test_CombineNexuses(unittest.TestCase):
    """Test combine_nexuses"""
    def setUp(self):
        self.nex1 = NexusReader()
        self.nex1.read_string(
            """Begin data;
            Dimensions ntax=2 nchar=1;
            Format datatype=standard symbols="12" gap=-;
            Matrix
            Harry              1
            Simon              2
            ;"""
        )
        self.nex2 = NexusReader()
        self.nex2.read_string(
            """Begin data;
            Dimensions ntax=2 nchar=1;
            Format datatype=standard symbols="34" gap=-;
            Matrix
            Harry              3
            Simon              4
            ;"""
        )
        self.nex3 = NexusReader()
        self.nex3.read_string(
            """Begin data;
            Dimensions ntax=3 nchar=1;
            Format datatype=standard symbols="345" gap=-;
            Matrix
            Betty              3
            Boris              4
            Simon              5
            ;"""
        )
    
    def test_failure_on_nonlist_1(self):
        self.assertRaises(TypeError, combine_nexuses, "I am not a list")
        
    def test_failure_on_nonlist_2(self):
        self.assertRaises(TypeError, combine_nexuses, ["hello",]) # should be NexusReader instances
        
    def test_combine_simple(self):
        newnex = combine_nexuses([self.nex1, self.nex2])
        assert newnex.data['1.1']['Harry'] == '1'
        assert newnex.data['1.1']['Simon'] == '2'
        assert newnex.data['2.1']['Harry'] == '3'
        assert newnex.data['2.1']['Simon'] == '4'
    
    def test_combine_simple_generated_matrix(self):
        newnex = combine_nexuses([self.nex1, self.nex2])
        assert re.search(r"""\bSimon\s+24\b""", newnex.write())
        assert re.search(r"""\bHarry\s+13\b""", newnex.write())
    
    def test_combine_simple_generated_formatline(self):
        newnex = combine_nexuses([self.nex1, self.nex2])
        assert re.search(r"""\bNTAX=2\b""", newnex.write())
        assert re.search(r"""\bNCHAR=2\b""", newnex.write())
        assert re.search(r'\sSYMBOLS="1234"[\s;]', newnex.write())
        
    def test_combine_missing(self):
        newnex = combine_nexuses([self.nex1, self.nex3])
        assert newnex.data['1.1']['Harry'] == '1'
        assert newnex.data['1.1']['Simon'] == '2'
        assert newnex.data['2.1']['Betty'] == '3'
        assert newnex.data['2.1']['Boris'] == '4'
        
    def test_combine_missing_generated_matrix(self):
        newnex = combine_nexuses([self.nex1, self.nex3])
        assert re.search(r"""\bSimon\s+25\b""", newnex.write())
        assert re.search(r"""\bHarry\s+1\\?\b""", newnex.write())
        assert re.search(r"""\bBetty\s+\?3\b""", newnex.write())
        assert re.search(r"""\bBoris\s+\?4\b""", newnex.write())
        
    def test_combine_missing_generated_formatline(self):
        newnex = combine_nexuses([self.nex1, self.nex3])
        assert re.search(r"""\bNTAX=4\b""", newnex.write())
        assert re.search(r"""\bNCHAR=2\b""", newnex.write())
        assert re.search(r'\sSYMBOLS="12345"[\s;]', newnex.write())

    def test_combine_with_character_labels(self):
        n1 = NexusReader()
        n1.read_string(
            """
            BEGIN DATA;
                DIMENSIONS NTAX=3 NCHAR=3;
                FORMAT DATATYPE=STANDARD MISSING=0 GAP=-  SYMBOLS="123";
                CHARSTATELABELS
            		1 char1,
            		2 char2,
            		3 char3
            ;
            MATRIX
            Tax1         123
            Tax2         123
            Tax3         123
            ;
            """
        )
        n2 = NexusReader()
        n2.read_string(
            """
            BEGIN DATA;
                DIMENSIONS NTAX=3 NCHAR=3;
                FORMAT DATATYPE=STANDARD MISSING=0 GAP=-  SYMBOLS="456";
                CHARSTATELABELS
            		1 char1,
            		2 char2,
            		3 char3
            ;
            MATRIX
            Tax1         456
            Tax2         456
            Tax3         456
            ;
            """
        )
        newnex = combine_nexuses([n1, n2])
        assert re.search(r"""\bNTAX=3\b""", newnex.write())
        assert re.search(r"""\bNCHAR=6\b""", newnex.write())
        assert re.search(r'\sSYMBOLS="123456"[\s;]', newnex.write())
        
        for tax in [1,2,3]:
            assert re.search(r"""\bTax%d\s+123456\b""" % tax, newnex.write())
        
        counter = 1
        for nex_id in [1,2]:
            for char_id in [1,2,3]:
                assert re.search(
                    r"""\b%d\s+%d.char%d\b""" % (counter, nex_id, char_id), 
                    newnex.write(charblock=True)
                )
                counter += 1


if __name__ == '__main__':
    unittest.main()

import os
import unittest

from nexus import NexusReader
from nexus.tools import multistatise

EXAMPLE_DIR = os.path.join(os.path.dirname(__file__), '../../examples')

class Test_Multistatise(unittest.TestCase):
    """Test multistatise"""
    def setUp(self):
        self.nex = NexusReader()
        self.nex.read_string(
        """Begin data;
        Dimensions ntax=4 nchar=4;
        Format datatype=standard symbols="01" gap=-;
        Matrix
        Harry              1000
        Simon              0100
        Betty              0010
        Louise             0001
        ;""")
        self.nex = multistatise(self.nex)
        
    def test_nexusreader_transformation(self):
        assert isinstance(self.nex, NexusReader), "Nexus_obj should be a NexusReader instance"

    def test_block_find(self):
        assert 'data' in self.nex.blocks
    
    def test_ntaxa_recovery(self):
        assert self.nex.data.ntaxa == 4
        
    def test_nchar_recovery(self):
        assert self.nex.data.nchar == 1
        
    def test_matrix(self):
        assert self.nex.data.matrix['Harry'][0] == 'A'
        assert self.nex.data.matrix['Simon'][0] == 'B'
        assert self.nex.data.matrix['Betty'][0] == 'C'
        assert self.nex.data.matrix['Louise'][0] == 'D'
    
    
    def test_regression_include_invisible_taxa(self):
        """Include taxa that have no entries"""
        data = """
        #NEXUS
        
        BEGIN DATA;
            DIMENSIONS  NTAX=15 NCHAR=7;
            FORMAT DATATYPE=STANDARD MISSING=? GAP=- INTERLEAVE=YES;
        MATRIX
        
        Gertrude                0000001
        Debbie                  0001000
        Zarathrustra            0000000
        Christie                0010000
        Benny                   0100000
        Bertha                  0100000
        Craig                   0010000
        Fannie-May              0000010
        Charles                 0010000
        Annik                   1000000
        Frank                   0000010
        Amber                   1000000
        Andreea                 1000000
        Edward                  0000100
        Donald                  0001000
        ;
        END;
        """
        
        nex = NexusReader()
        nex.read_string(data)
        msnex = multistatise(nex)
        
        for taxon,sites in msnex.data.matrix.items():
            if taxon[0] == 'Z':
                continue # will check later
            
            # first letter of taxa name is the expected character state
            assert taxon[0] == sites[0], "%s should be %s not %s" % (taxon, taxon[0], sites[0])
        # deal with completely missing taxa
        assert 'Zarathrustra' in msnex.data.matrix
        assert msnex.data.matrix['Zarathrustra'][0] == '?'


if __name__ == '__main__':
    unittest.main()

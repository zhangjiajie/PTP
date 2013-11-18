import re
import unittest

from nexus import NexusReader
from nexus.tools import binarise
from nexus.tools.binarise import _recode_to_binary

class Test_Recode_To_Binary(unittest.TestCase):
    
    def test_error_on_integer(self):
        orig = {'Maori': 1, 'Dutch': '1', 'Latin': '1'}
        with self.assertRaises(ValueError):
            recoded = _recode_to_binary(orig)
            
    def test_error_on_badvalue(self):
        orig = {'Maori': None, 'Dutch': '1', 'Latin': '1'}
        with self.assertRaises(ValueError):
            recoded = _recode_to_binary(orig)
    
    def test_one(self):
        orig = {'Maori': '1', 'Dutch': '1', 'Latin': '1'}
        recoded = _recode_to_binary(orig)
        assert recoded['Maori'] == '1'
        assert recoded['Dutch'] == '1'
        assert recoded['Latin'] == '1'
        
    def test_two_state(self):
        orig = {'Maori': '1', 'Dutch': '2', 'Latin': '1'}
        recoded = _recode_to_binary(orig)
        assert recoded['Maori'] == '10'
        assert recoded['Dutch'] == '01'
        assert recoded['Latin'] == '10'
        
    def test_three_state(self):
        orig = {'Maori': '1', 'Dutch': '2', 'Latin': '3'}
        recoded = _recode_to_binary(orig)
        assert recoded['Maori'] == '100'
        assert recoded['Dutch'] == '010'
        assert recoded['Latin'] == '001'
        
    def test_absent_state(self):
        orig = {'Maori': '0', 'Dutch': '2', 'Latin': '1'}
        recoded = _recode_to_binary(orig)
        assert recoded['Maori'] == '00', recoded
        assert recoded['Dutch'] == '01', recoded
        assert recoded['Latin'] == '10', recoded
        
    def test_absent_state_but_keep_zero(self):
        orig = {'Maori': '0', 'Dutch': '2', 'Latin': '1'}
        recoded = _recode_to_binary(orig, keep_zero=True)
        assert recoded['Maori'] == '100', recoded
        assert recoded['Dutch'] == '001', recoded
        assert recoded['Latin'] == '010', recoded
        
    def test_missing_state(self):
        orig = {'Maori': '?', 'Dutch': '2', 'Latin': '1'}
        recoded = _recode_to_binary(orig)
        assert recoded['Maori'] == '00'
        assert recoded['Dutch'] == '01'
        assert recoded['Latin'] == '10'
        
    def test_gap_state(self):
        orig = {'Maori': '-', 'Dutch': '2', 'Latin': '1'}
        recoded = _recode_to_binary(orig)
        assert recoded['Maori'] == '00'
        assert recoded['Dutch'] == '01'
        assert recoded['Latin'] == '10'
        
    def test_noncontiguous_states(self):
        orig = {'Maori': '1', 'Dutch': '5', 'Latin': '9'}
        recoded = _recode_to_binary(orig)
        assert recoded['Maori'] == '100'
        assert recoded['Dutch'] == '010'
        assert recoded['Latin'] == '001'
    
    @unittest.expectedFailure
    def test_polymorphic_states(self):
        orig = {'Maori': '1,3', 'Dutch': '2', 'Latin': '3'}
        recoded = _recode_to_binary(orig)
        assert recoded['Maori'] == '101', recoded
        assert recoded['Dutch'] == '010', recoded
        assert recoded['Latin'] == '001', recoded



class Test_Binarise(unittest.TestCase):
    def setUp(self):
        self.nex = NexusReader()
        self.nex.read_string(
        """Begin data;
        Dimensions ntax=3 nchar=2;
        Format datatype=standard symbols="01" gap=-;
        Charstatelabels
            1 char1, 2 char2;
        Matrix
        Maori               14
        Dutch               25
        Latin               36
        ;""")
        self.nex = binarise(self.nex)
    
    def test_to_binary(self):
        """Test Nexus -> Binary: Two Character"""
        expected = {
            'char1_0': {"Maori": '1', "Dutch": "0", "Latin": "0"},
            'char1_1': {"Maori": '0', "Dutch": "1", "Latin": "0"},
            'char1_2': {"Maori": '0', "Dutch": "0", "Latin": "1"},
            'char2_0': {"Maori": '1', "Dutch": "0", "Latin": "0"},
            'char2_1': {"Maori": '0', "Dutch": "1", "Latin": "0"},
            'char2_2': {"Maori": '0', "Dutch": "0", "Latin": "1"},
        }
        
        for char, data in expected.items():
            for taxon, exp_value in data.items():
                assert self.nex.data[char][taxon] == exp_value
    
    def test_to_binary_nchar(self):
        """Test Nexus -> Binary: Number of Characters"""
        assert len(self.nex.characters) == 6
        
    def test_to_binary_symbollist(self):
        """Test Nexus -> Binary: Update Symbol List"""
        # check symbol list was updated
        assert len(self.nex.symbols) == 2
        assert '1' in self.nex.symbols
        assert '0' in self.nex.symbols
        
    def test_to_binary_nexus(self):
        """Test Nexus -> Binary: Nexus"""
        nexus = self.nex.make_nexus(interleave=False)
        assert re.search("Dutch\s+010010", nexus)
        assert re.search("Maori\s+100100", nexus)
        assert re.search("Latin\s+001001", nexus)
        


if __name__ == '__main__':
    unittest.main()

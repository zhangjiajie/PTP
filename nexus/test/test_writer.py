import re
import copy
import unittest
from nexus import NexusWriter

data = {
    'char1': {'French': 1, 'English': 2, 'Latin': 3},
    'char2': {'French': 4, 'English': 5, 'Latin': 6},
}

class Test_NexusWriter_1(unittest.TestCase):
    def setUp(self):
        self.nex = NexusWriter()
        
    def test_char_adding1(self):
        """Test Character Addition 1"""
        for tx, value in data['char1'].items():
            self.nex.add(tx, 'char1', value)
        assert self.nex.data['char1']['French'] == '1'
        assert self.nex.data['char1']['English'] == '2'
        assert self.nex.data['char1']['Latin'] == '3'
        
    def test_char_adding2(self):
        """Test Character Addition 2"""
        for tx, value in data['char2'].items():
            self.nex.add(tx, 'char2', value)
        assert self.nex.data['char2']['French'] == '4'
        assert self.nex.data['char2']['English'] == '5'
        assert self.nex.data['char2']['Latin'] == '6'
        

class Test_NexusWriter_2(unittest.TestCase):
    def setUp(self):
        self.nex = NexusWriter()
        for char, b in data.items():
            for taxon, value in b.items():
                self.nex.add(taxon, char, value)

    def test_nexus_noninterleave(self):
        """Test Nexus Generation - Non-Interleaved"""
        n = self.nex.make_nexus(interleave=False)
        assert re.search("#NEXUS", n)
        assert re.search("BEGIN DATA;", n)
        assert re.search("DIMENSIONS NTAX=3 NCHAR=2;", n)
        assert re.search("MATRIX", n)
        assert re.search("Latin\s+36", n)
        assert re.search("French\s+14", n)
        assert re.search("English\s+25", n)
        assert re.search("FORMAT.*MISSING\=(.+?)", n).groups()[0] == '?'
        assert re.search("FORMAT.*DATATYPE\=(\w+)\s", n).groups()[0] == 'STANDARD'
        assert re.search('FORMAT.*SYMBOLS\="(\d+)";', n).groups()[0] == '123456'
        
    def test_nexus_charblock(self):
        """Test Nexus Generation - with characters block"""
        n = self.nex.make_nexus(charblock=True)
        assert re.search("#NEXUS", n)
        assert re.search("BEGIN DATA;", n)
        assert re.search("DIMENSIONS NTAX=3 NCHAR=2;", n)
        assert re.search("CHARSTATELABELS", n)
        assert re.search("1 char1,", n)
        assert re.search("2 char2", n)
        assert re.search("MATRIX", n)
        assert re.search("Latin\s+36", n)
        assert re.search("French\s+14", n)
        assert re.search("English\s+25", n)
        assert re.search("FORMAT.*MISSING\=(.+?)", n).groups()[0] == '?'
        assert re.search("FORMAT.*DATATYPE\=(\w+)\s", n).groups()[0] == 'STANDARD'
        assert re.search('FORMAT.*SYMBOLS\="(\d+)";', n).groups()[0] == '123456'
    
    def test_nexus_interleave(self):
        """Test Nexus Generation - Interleaved"""
        n = self.nex.make_nexus(interleave=True)
        assert re.search("#NEXUS", n)
        assert re.search("BEGIN DATA;", n)
        assert re.search("DIMENSIONS NTAX=3 NCHAR=2;", n)
        assert re.search("MATRIX", n)
        # char1
        assert re.search("Latin\s+3", n)
        assert re.search("French\s+1", n)
        assert re.search("English\s+2", n)
        # char2
        assert re.search("Latin\s+6", n)
        assert re.search("French\s+4", n)
        assert re.search("English\s+5", n)
        
        assert re.search("FORMAT.*MISSING\=(.+?)", n).groups()[0] == '?'
        assert re.search("FORMAT.*DATATYPE\=(\w+)\s", n).groups()[0] == 'STANDARD'
        assert re.search("FORMAT.*(INTERLEAVE)", n).groups()[0] == 'INTERLEAVE'
        assert re.search('FORMAT.*SYMBOLS\="(\d+)";', n).groups()[0] == '123456'


class RegressionTests(unittest.TestCase):
    def test_regression_format_string_has_datatype_first(self):
        """Regression: Format string should contain 'datatype' as the first element"""
        # SplitsTree complains otherwise.
        nex = NexusWriter()
        for char, b in data.items():
            for taxon, value in b.items():
                nex.add(taxon, char, value)
        out = nex.make_nexus()
        assert "FORMAT DATATYPE=STANDARD" in out
    
    def test_regression_format_string_has_quoted_symbols(self):
        """Regression: Symbols in the format string should be quoted"""
        nex = NexusWriter()
        for char, b in data.items():
            for taxon, value in b.items():
                nex.add(taxon, char, value)
        out = nex.make_nexus()
        assert 'SYMBOLS="123456"' in out

    
if __name__ == '__main__':
    unittest.main()

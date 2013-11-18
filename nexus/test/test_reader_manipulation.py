"""Tests for nexus reading manipulation"""
import os
import re
import unittest
from nexus import NexusReader
from nexus.reader import GenericHandler, DataHandler, TreeHandler

EXAMPLE_DIR = os.path.join(os.path.dirname(__file__), '../examples')

class Test_Manipulation_Data(unittest.TestCase):
    """Test the manipulation of data in the NexusReader"""
    
    def setUp(self):
        self.nex = NexusReader(os.path.join(EXAMPLE_DIR, 'example.nex'))
    
    def test_add_taxa(self):
        assert self.nex.data.ntaxa == 4
        self.nex.data.add_taxon('Elvis', ['1', '2'])
        assert self.nex.data.ntaxa == 5
        assert self.nex.data.matrix['Elvis'] == ['1', '2']
        assert 'Elvis' in self.nex.data.taxa
        assert 'Elvis' in self.nex.data.matrix
        expected_patterns = [
            '^begin data;$',
            '^\s+dimensions ntax=5 nchar=2;$',
            '^\s+format datatype=standard symbols="01" gap=-;$',
            '^matrix$',
            '^Simon\s+01$',
            '^Louise\s+11$',
            '^Betty\s+10$',
            '^Harry\s+00$',
            '^Elvis\s+12$',
            '^\s+;$',
            '^end;$',
        ]
        written = self.nex.write()
        for expected in expected_patterns:
            assert re.search(expected, written, re.MULTILINE), 'Expected "%s"' % expected
    
    def test_delete_taxa(self):
        assert self.nex.data.ntaxa == 4
        self.nex.data.del_taxon('Simon')
        assert self.nex.data.ntaxa == 3
        
        assert 'Simon' not in self.nex.data.taxa
        assert 'Simon' not in self.nex.data.matrix
        
        expected_patterns = [
            '^begin data;$',
            '^\s+dimensions ntax=3 nchar=2;$',
            '^\s+format datatype=standard symbols="01" gap=-;$',
            '^matrix$',
            '^Louise\s+11$',
            '^Betty\s+10$',
            '^Harry\s+00$',
            '^\s+;$',
            '^end;$',
        ]
        written = self.nex.write()
        for expected in expected_patterns:
            assert re.search(expected, written, re.MULTILINE), 'Expected "%s"' % expected
        
        # should NOT be here
        assert re.search('^Simon\s+01$', written, re.MULTILINE) == None, \
            'Expected Taxon "Simon" to be Deleted'
        
    def test_add_character(self):
        pass
        
    def test_delete_character(self):
        pass

    def test_edit_charlabels(self):
        pass
        

    
# TreeHandler
# self.translators = {}
# self.attributes = []
# self.taxa = []
# self.trees = []
# ntrees
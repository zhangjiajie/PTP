"""Tests for utils in bin directory"""
import os
import unittest

from nexus import NexusReader
from nexus.bin.nexus_anonymise import anonymise, hash

EXAMPLE_DIR = os.path.join(os.path.dirname(__file__), '../../examples')

class Test_Anonymise(unittest.TestCase):
    
    def test_anonymise_data(self):
        nex = NexusReader(os.path.join(EXAMPLE_DIR, 'example.nex'))
        nex = anonymise(nex)
        for old_taxon in ['Harry', 'Simon', 'Betty', 'Louise']:
            assert old_taxon not in nex.data.matrix, '%s should have been anonymised' % old_taxon
        
        assert nex.data.matrix['894a76c65225a9812d31ff75edf38feb'] == ['1', '0']
        assert nex.data.matrix['a0434190848c0d64332dce12a8a27961'] == ['0', '0']
        assert nex.data.matrix['bbf0da40d536d862e184a6eccb433a73'] == ['0', '1']
        assert nex.data.matrix['d24eb4091c14b87b6cd0bd94fd0704be'] == ['1', '1']

    def test_anonymise_data_with_labels(self):
        nex = NexusReader(os.path.join(EXAMPLE_DIR, 'example2.nex'))
        nex = anonymise(nex)
        for old_taxon in ['John', 'Paul', 'George', 'Ringo']:
            assert old_taxon not in nex.data.matrix, '%s should have been anonymised' % old_taxon
            
        # check data block
        assert nex.data.matrix['e5b17eef97cc46782008af612450c650'] == ['a', 'c', 't', 'g']
        assert nex.data.matrix['fc38854f2edead0c7a799a35de58bd60'] == ['a', 'c', 't', 'g']
        assert nex.data.matrix['c12b8557f6be8a49aef5fcee99c40724'] == ['a', 'c', 't', 'g']
        assert nex.data.matrix['de57ceacb7eb62511210794be2cff5ab'] == ['a', 'c', 't', 'g']
        
        # check taxa block
        assert 'e5b17eef97cc46782008af612450c650' in nex.taxa.taxa
        assert 'fc38854f2edead0c7a799a35de58bd60' in nex.taxa.taxa
        assert 'c12b8557f6be8a49aef5fcee99c40724' in nex.taxa.taxa
        assert 'de57ceacb7eb62511210794be2cff5ab' in nex.taxa.taxa
        
    def test_anonymise_data_with_interleave(self):
        nex = NexusReader(os.path.join(EXAMPLE_DIR, 'example3.nex'))
        nex = anonymise(nex)
        for old_taxon in ['Harry', 'Simon']:
            assert old_taxon not in nex.data.matrix, '%s should have been anonymised' % old_taxon
        assert nex.data.matrix['6eb7148a2d4155085e517979410b9f23'] == ['0', '1', '2', '3', '4', '5']
        assert nex.data.matrix['698de77f637e7fae18ead22f2172102a'] == ['0', '1', '2', '3', '4', '5']
    
    def test_anonymise_untranslated_trees(self):
        # NOT IMPLEMENTED
        nex = NexusReader(os.path.join(EXAMPLE_DIR, 'example.trees'))
        self.assertRaises(NotImplementedError, anonymise, nex)
    
    def test_anonymise_translated_trees(self):
        nex = NexusReader(os.path.join(EXAMPLE_DIR, 'example-translated.trees'))
        nex = anonymise(nex)
        expected = ['Chris', 'Bruce', 'Tom', 'Henry', 'Timothy', 'Mark', 'Simon', 'Fred', 'Kevin', 'Roger', 'Michael', 'Andrew', 'David']
        assert len(nex.trees.taxa) == len(expected)
        for taxon in expected:
            hashtaxon = hash(os.path.join(EXAMPLE_DIR, 'example-translated.trees'), taxon)
            assert hashtaxon in nex.trees.taxa
    
    def test_anonymise_beast_treefile(self):
        nex = NexusReader(os.path.join(EXAMPLE_DIR, 'example-beast.trees'))
        nex = anonymise(nex)
        expected = [
            "R1", "B2", "S3", "T4", "A5", "E6", "U7", "T8", "T9", "F10", "U11", 
            "T12", "N13", "F14", "K15", "N16", "I17", "L18", "S19", "T20", "V21", 
            "R22", "M23", "H24", "M25", "M26", "M27", "R28", "T29", "M30", "P31", 
            "T32", "R33", "P34", "R35", "W36", "F37", "F38"
        ]
        # check taxa block
        for taxon in expected:
            assert taxon not in nex.taxa.taxa, '%s should have been anonymised' % taxon
            hashtaxon = hash(os.path.join(EXAMPLE_DIR, 'example-beast.trees'), taxon)
            assert hashtaxon in nex.taxa.taxa
        
        # check trees block
        for taxon in expected:
            assert taxon not in nex.trees.taxa, '%s should have been anonymised' % taxon
            hashtaxon = hash(os.path.join(EXAMPLE_DIR, 'example-beast.trees'), taxon)
            assert hashtaxon in nex.trees.taxa
        
    def test_anonymise_characters(self):
        nex = NexusReader(os.path.join(EXAMPLE_DIR, 'example-characters.nex'))
        nex = anonymise(nex)
        
        expected_taxa = ['A', 'B', 'C', 'D', 'E']
        # check taxa block
        for taxon in expected_taxa:
            assert taxon not in nex.taxa.taxa, '%s should have been anonymised' % taxon
            hashtaxon = hash(os.path.join(EXAMPLE_DIR, 'example-characters.nex'), taxon)
            assert hashtaxon in nex.taxa.taxa
        
        # check characters block
        for taxon in expected_taxa:
            assert taxon not in nex.data.taxa, '%s should have been anonymised' % taxon
            hashtaxon = hash(os.path.join(EXAMPLE_DIR, 'example-characters.nex'), taxon)
            assert hashtaxon in nex.data.taxa
    
    
    

if __name__ == '__main__':
    unittest.main()

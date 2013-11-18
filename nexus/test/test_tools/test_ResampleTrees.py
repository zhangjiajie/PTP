"""Tests for utils in bin directory"""
import os
import unittest

from nexus import NexusReader
from nexus.bin.nexus_treemanip import run_resample

EXAMPLE_DIR = os.path.join(os.path.dirname(__file__), '../../examples')

class Test_ResampleTrees(unittest.TestCase):
    """Test nexus_treemanip.run_resample"""
    def setUp(self):
        self.nexus = NexusReader(os.path.join(EXAMPLE_DIR, 'example.trees'))
        
    def test_resample(self):
        newnex = run_resample(2, self.nexus)
        assert len(newnex.trees.trees) == 1
        
    def test_resample_one(self):
        newnex = run_resample(1, self.nexus)
        assert len(newnex.trees.trees) == 3
    

if __name__ == '__main__':
    unittest.main()

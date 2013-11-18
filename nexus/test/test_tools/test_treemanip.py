"""Tests for utils in bin directory"""
import os
import unittest

from nexus import NexusReader
from nexus.bin.nexus_treemanip import TreeListException
from nexus.bin.nexus_treemanip import parse_deltree, run_deltree
from nexus.bin.nexus_treemanip import run_resample
from nexus.bin.nexus_treemanip import run_removecomments
from nexus.bin.nexus_treemanip import run_detranslate
from nexus.bin.nexus_treemanip import run_random

EXAMPLE_DIR = os.path.join(os.path.dirname(__file__), '../../examples')

class Test_TreeManip_run_deltree(unittest.TestCase):
    
    def test_run_deltree(self):
        nex = NexusReader(os.path.join(EXAMPLE_DIR, 'example.trees'))
        new_nex = run_deltree('2', nex, do_print=False)
        assert len(new_nex.trees.trees) == 2
        assert new_nex.trees.ntrees == 2
        assert new_nex.trees[0].startswith('tree tree.0.1065.603220')
        assert new_nex.trees[1].startswith('tree tree.20000.883.396049')
    

class Test_TreeManip_run_resample(unittest.TestCase):

    def test_run_resample_1(self):
        # shouldn't resample anything..
        nex = NexusReader(os.path.join(EXAMPLE_DIR, 'example.trees'))
        new_nex = run_resample('1', nex, do_print=False)
        assert len(new_nex.trees.trees) == 3
        assert new_nex.trees.ntrees == 3
        assert new_nex.trees[0].startswith('tree tree.0.1065.603220')
        assert new_nex.trees[1].startswith('tree tree.10000.874.808756')
        assert new_nex.trees[2].startswith('tree tree.20000.883.396049')


class Test_TreeManip_run_removecomments(unittest.TestCase):
    
    def test_run_removecomments(self):
        nex = NexusReader(os.path.join(EXAMPLE_DIR, 'example-beast.trees'))
        new_nex = run_removecomments(nex, do_print=False)
        assert '[&lnP=-15795.47019648783]' not in new_nex.trees[0]

        
class Test_TreeManip_run_randomise(unittest.TestCase):
    
    def test_failure_on_nonint(self):
        nex = NexusReader(os.path.join(EXAMPLE_DIR, 'example-translated.trees'))
        self.assertRaises(ValueError, run_random, 'fudge', nex)
        
    def test_run_randomise_sample1(self):
        nex = NexusReader(os.path.join(EXAMPLE_DIR, 'example-translated.trees'))
        new_nex = run_random(1, nex)
        assert new_nex.trees.ntrees == len(new_nex.trees.trees) == 1
        
    def test_run_randomise_sample2(self):
        nex = NexusReader(os.path.join(EXAMPLE_DIR, 'example-translated.trees'))
        new_nex = run_random(2, nex)
        assert new_nex.trees.ntrees == len(new_nex.trees.trees) == 2
    
    def test_run_randomise_sample_toobig(self):
        # raises ValueError, sample size too big (only 3 trees in this file)
        nex = NexusReader(os.path.join(EXAMPLE_DIR, 'example-translated.trees'))
        self.assertRaises(ValueError, run_random, 10, nex)
        
        
if __name__ == '__main__':
    unittest.main()

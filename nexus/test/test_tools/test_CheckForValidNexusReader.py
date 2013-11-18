import os
import unittest

from nexus import NexusReader, NexusFormatException
from nexus.tools import check_for_valid_NexusReader

EXAMPLE_DIR = os.path.join(os.path.dirname(__file__), '../../examples')

class Test_CheckForValidNexusReader(unittest.TestCase):
    """Test check_for_valid_NexusReader"""
    
    def test_valid_NexusReader(self):
        check_for_valid_NexusReader(NexusReader())
    
    def test_failure_on_nonnexus_1(self):
        self.assertRaises(TypeError, check_for_valid_NexusReader, "I am not a nexus")
    
    def test_failure_on_nonnexus_2(self):
        self.assertRaises(TypeError, check_for_valid_NexusReader, 1)
        
    def test_failure_on_nonnexus_3(self):
        self.assertRaises(TypeError, check_for_valid_NexusReader, [1, 2 ,3])
    
    def test_failure_on_required_block_one(self):
        nexus_obj = NexusReader(os.path.join(EXAMPLE_DIR, 'example.nex'))
        with self.assertRaises(NexusFormatException):
            check_for_valid_NexusReader(nexus_obj, ['trees'])
    
    def test_failure_on_required_block_two(self):
        nexus_obj = NexusReader(os.path.join(EXAMPLE_DIR, 'example2.nex'))
        with self.assertRaises(NexusFormatException):
            check_for_valid_NexusReader(nexus_obj, ['r8s'])
    
    def test_valid_with_required_block_one(self):
        nexus_obj = NexusReader(os.path.join(EXAMPLE_DIR, 'example.nex'))
        check_for_valid_NexusReader(nexus_obj, ['data'])
        
    def test_valid_with_required_block_two(self):
        nexus_obj = NexusReader(os.path.join(EXAMPLE_DIR, 'example2.nex'))
        check_for_valid_NexusReader(nexus_obj, ['data', 'taxa'])

if __name__ == '__main__':
    unittest.main()

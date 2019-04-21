#!/usr/bin/env python
import sys
import os
from random import sample

from nexus import NexusReader, VERSION, NexusFormatException
from nexus.tools import check_for_valid_NexusReader


__author__ = 'Simon Greenhill <simon@simon.net.nz>'
__doc__ = """treemanip.py - python-nexus tools v%(version)s
Performs some functions on trees
""" % {'version': VERSION,}

__usage__ = """
Deleting trees:
    nexus_treemanip.py -d 1 old.trees new.trees   - delete tree #1
    nexus_treemanip.py -d 1-5 old.trees new.trees - delete trees #1-5
    nexus_treemanip.py -d 1,5 old.trees new.trees - delete trees #1 and #5
    nexus_treemanip.py -d 1,4,20-30 old.trees new.trees - delete trees #1, #4, #20-30
    
Resampling trees:
    nexus_treemanip.py -r 10 old.trees new.trees   - resample every 10th tree

Remove translate block:
    nexus_treemanip.py -t old.trees new.trees
    
Remove comments:
    nexus_treemanip.py -c old.trees new.trees
"""

class TreeListException(Exception):
    """Generic Exception for Tree Lists"""
    def __init__(self, arg):
        Exception.__init__(self, arg)
        self.arg = arg

def parse_deltree(dstring):
    """
    Returns a list of trees to be deleted
    
    :param dstring: A string
    :type dstring: string 
    
    :return: A list of trees to be deleted.
    :raises TreeListException: if dstring is invalid.
    
    >>> parse_deltree('1')
    [1]
    >>> parse_deltree('1,2,3')
    [1, 2, 3]
    >>> parse_deltree('1,3,5')
    [1, 3, 5]
    >>> parse_deltree('1-10')
    [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    >>> parse_deltree('1,3,4-6')
    [1, 3, 4, 5, 6]
    >>> parse_deltree('1,3,4-6,8,9-10')
    [1, 3, 4, 5, 6, 8, 9, 10]

    # alternate syntax
    >>> parse_deltree('1:10')
    [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    >>> parse_deltree('1,3,4:6')
    [1, 3, 4, 5, 6]
    >>> parse_deltree('1,3,4:6,8,9:10')
    [1, 3, 4, 5, 6, 8, 9, 10]
    """
    out = []
    for token in dstring.split(','):
        token = token.replace(':', '-')
        if '-' in token:
            try:
                start, stop = token.split("-")
                out.extend([x for x in range(int(start), int(stop)+1)])
            except (ValueError, IndexError):
                raise TreeListException("'%s' is not a valid token for a tree list" % token)
        else:
            try:
                out.append(int(token))
            except (ValueError):
                raise TreeListException("'%s' is not a valid token for a tree list" % token)
    return sorted(out)


def run_deltree(deltree, nexus_obj, do_print=False):
    """
    Returns a list of trees to be deleted
    
    :param deltree: A string of trees to be deleted.
    :type deltree: String 
    
    :param nexus_obj: A `NexusReader` instance
    :type nexus_obj: NexusReader 
    
    :param do_print: flag to print() logging information or not
    :type do_print: Boolean
    
    :return: A NexusReader instance with the given trees removed.
    
    :raises AssertionError: if nexus_obj is not a nexus
    :raises NexusFormatException: if nexus_obj does not have a `trees` block
    """
    check_for_valid_NexusReader(nexus_obj, required_blocks=['trees'])
        
    new = []
    delitems = parse_deltree(deltree)
    
    if do_print: 
        print('Deleting: %d trees' % len(delitems))
        
    for index, tree in enumerate(nexus_obj.trees, 1):
        if index in delitems:
            if do_print: print('Deleting tree %d' % index)
        else:
            new.append(tree)
    nexus_obj.trees.trees = new
    return nexus_obj


def run_resample(resample, nexus_obj, do_print=False):
    """
    Resamples the trees in a nexus
    
    :param resample: Resample every `resample` trees
    :type resample: Integer 
    
    :param nexus_obj: A `NexusReader` instance
    :type nexus_obj: NexusReader 
    
    :param do_print: flag to print() logging information or not
    :type do_print: Boolean
    
    :return: A NexusReader instance with the given trees removed.
    
    :raises AssertionError: if nexus_obj is not a nexus
    :raises NexusFormatException: if nexus_obj does not have a `trees` block
    """
    check_for_valid_NexusReader(nexus_obj, required_blocks=['trees'])
    
    new = []
    try:
        every = int(resample)
    except ValueError:
        sys.exit("Invalid resample option %s - should be an integer" % resample)
    
    if do_print: 
        print('Resampling ever %d trees' % every)
    
    ignore_count = 0 
    for index, tree in enumerate(nexus_obj.trees, 1):
        if index % every == 0:
            new.append(tree)
        else:
            ignore_count += 1
    
    if do_print: 
        print("Ignored %d trees" % ignore_count)
    
    nexus_obj.trees.trees = new
    return nexus_obj
    
def run_removecomments(nexus_obj, do_print=False):
    """
    Removes comments from the trees in a nexus
    
    :param nexus_obj: A `NexusReader` instance
    :type nexus_obj: NexusReader 
    
    :param do_print: flag to print() logging information or not
    :type do_print: Boolean
    
    :return: A NexusReader instance with the comments removed.
    
    :raises AssertionError: if nexus_obj is not a nexus
    :raises NexusFormatException: if nexus_obj does not have a `trees` block
    """
    check_for_valid_NexusReader(nexus_obj, required_blocks=['trees'])
    
    new = []
    for index, tree in enumerate(nexus_obj.trees, 1):
        new.append(nexus_obj.trees.remove_comments(tree))
    
    if do_print:
        print("Removed comments")
    
    nexus_obj.trees.trees = new
    return nexus_obj

def run_detranslate(nexus_obj, do_print=False):
    """
    Removes comments from the trees in a nexus
    
    :param nexus_obj: A `NexusReader` instance
    :type nexus_obj: NexusReader 
    
    :param do_print: flag to print() logging information or not
    :type do_print: Boolean
    
    :return: A NexusReader instance with the comments removed.
    
    :raises AssertionError: if nexus_obj is not a nexus
    :raises NexusFormatException: if nexus_obj does not have a `trees` block
    """
    check_for_valid_NexusReader(nexus_obj, required_blocks=['trees'])
    nexus_obj.trees.detranslate()
    return nexus_obj

def run_random(num_trees, nexus_obj, do_print=False):
    """
    Returns a specified number (`num_trees`) of random trees from the nexus file.
    
    :param num_trees: The number of trees to resample
    :type num_trees: Integer
    
    :param nexus_obj: A `NexusReader` instance
    :type nexus_obj: NexusReader 
    
    :return: A NexusReader instance.
    
    :raises AssertionError: if nexus_obj is not a nexus
    :raises NexusFormatException: if nexus_obj does not have a `trees` block
    :raises ValueError: if num_trees is not an integer
    :raises ValueError: if num_trees is larger than population
    """
    check_for_valid_NexusReader(nexus_obj, required_blocks=['trees'])
    
    try:
        num_trees = int(num_trees)
    except ValueError:
        raise ValueError("num_trees should be an integer")
    
    if num_trees > nexus_obj.trees.ntrees:
        raise ValueError("Treefile only has %d trees in it." % nexus_obj.trees.ntrees)
    elif num_trees == nexus_obj.trees.ntrees:
        return nexus_obj # um. ok.
    else:
        if do_print:
            print("%d trees read... sampling %d" % (nexus_obj.trees.ntrees, num_trees))
        nexus_obj.trees.trees = sample(nexus_obj.trees.trees, num_trees)
    return nexus_obj
    
    
    

if __name__ == '__main__':
    #set up command-line options
    from optparse import OptionParser
    parser = OptionParser(usage="usage: %prog old.trees new.trees")
    parser.add_option("-d", "--deltree", dest="deltree", 
            action="store", default=False, 
            help="Remove the listed trees")
    parser.add_option("-r", "--resample", dest="resample", 
            action="store", default=False, 
            help="Resample the trees every Nth tree")
    parser.add_option("-n", "--random", dest="random", 
            action="store", default=False, 
            help="Randomly sample N trees from the treefile")
    parser.add_option("-c", "--removecomments", dest="removecomments", 
            action="store_true", default=False, 
            help="Remove comments from the trees")
    parser.add_option("-t", "--detranslate", dest="detranslate", 
            action="store_true", default=False, 
            help="Remove taxa translation block from the trees")
    parser.add_option("-q", "--quiet", dest="quiet", 
            action="store_true", default=False, 
            help="Be quiet (no logging information displayed)")
    options, args = parser.parse_args()
    
    try:
        nexusname = args[0]
    except IndexError:
        print(__doc__)
        print(__usage__)
        print("Author: %s\n" % __author__)
        parser.print_help()
        sys.exit()
        
    try:
        newnexus = args[1]
    except IndexError:
        newnexus = None
    
    nexus = NexusReader(nexusname)
    if 'trees' not in nexus.blocks:
        sys.exit("No trees found in file %s!" % nexusname)
    if nexus.trees.ntrees == 0:
        sys.exit("No trees found in found %s!" % nexusname)
    if options.quiet is False:
        print( "%d trees found with %d translated taxa" % \
            (nexus.trees.ntrees, len(nexus.trees.translators)))
    
    # Delete trees
    if options.deltree:
        nexus = run_deltree(options.deltree, nexus, options.quiet)
    
    # Resample trees
    if options.resample:
        nexus = run_resample(options.resample, nexus, options.quiet)
    
    # Randomly sample trees
    if options.random:
        nexus = run_random(options.random, nexus, options.quiet)
    
    # remove comments
    if options.removecomments:
        nexus = run_removecomments(nexus, options.quiet)
    
    # detranslate
    if options.detranslate:
        nexus = run_detranslate(nexus, options.quiet)
    
    if newnexus is not None:
        nexus.write_to_file(newnexus)
        if options.quiet is False:
            print("New nexus with %d trees written to %s" % (nexus.trees.ntrees, newnexus))
    else:
        print(nexus.write())
        

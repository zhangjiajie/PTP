python-nexus - Generic nexus (.nex, .trees) reader/writer for python
================================================================================


Copyright (c) 2009-2010, Simon J. Greenhill <simon@simon.net.nz>
 * http://simon.net.nz/python-nexus/
 * http://bitbucket.org/simongreenhill/python-nexus/


Description
================================================================================


python-nexus provides simple nexus file-format reading/writing tools, and a small
collection of nexus manipulation scripts.



Usage
================================================================================


Reading a Nexus:
    
>>> from nexus import NexusReader
>>> n = NexusReader()
>>> n.read_file('examples/example.nex')

or, more simply: 

>>> n = NexusReader('examples/example.nex')

You can also load from a string:

>>> n = NexusReader()
>>> n.read_string('#NEXUS\n\nbegin foo; ... end;')


NexusReader will load each of the nexus `blocks` it identifies using specific 
`handlers`. 

>>> n.blocks
{   'data': <NexusDataBlock: 2 characters from 4 taxa>, 
    'trees': <NexusTreeBlock: 2 trees>
}

A dictionary mapping blocks to handlers is available at .handlers:

>>> print n.handlers = {
    'data': DataHandler,
    'characters': DataHandler,
    'trees': TreeHandler,
    'taxa': TaxaHandler,
}

Any blocks that aren't in this dictionary will be parsed using GenericHandler.

NexusReader can then write the nexus to a string using .write() or to another 
file using .write_to_file(filename):

>>> nexus = n.write()
>>> n.write_to_file("mynewnexus.nex")

NOTE: if you want more fine-grained control over generating nexus files, then try
NexusWriter discussed below.


Block Handlers:
================================================================================


There are specific "Handlers" to parse certain known nexus blocks, including the
common 'data', 'trees', and 'taxa' blocks. Any blocks that are unknown will be 
parsed with GenericHandler.

ALL handlers extend the `GenericHandler` class and have the following methods.

* `parse(self, data)`
    parse is called by NexusReader to parse the contents of the block (in `data`)
    appropriately.

* `write(self)`
    write is called by NexusReader to write the contents of a block to a string 
    (i.e. for regenerating the nexus format for saving a file to disk)


All blocks have access to the following:

The raw block content (as a list of lines)
>>> n.blockname.block
[....]

A helper function to remove all the comments in a nexus file.
>>> n.block.remove_comments("hello [comment!] world")
"hello world"

To find out what file the nexus was loaded from:

>>> n.filename
>>> n.short_filename
'example.nex'


`generic` block handler
-----------------------

The generic block handler simply stores each line of the block in `.block`:

>>> n.blockname.block
['line1', 'line2', ... ]


`data` block handler
--------------------

These are the main blocks encountered in nexus files - and contain the data matrix.

So, given the following nexus file with a data block:

    #NEXUS 
    
    Begin data;
    Dimensions ntax=4 nchar=2;
    Format datatype=standard symbols="01" gap=-;
        Matrix
    Harry              00
    Simon              01
    Betty              10
    Louise             11
        ;
    End;
    
    begin trees;
        tree A = ((Harry:0.1,Simon:0.2):0.1,Betty:0.2):Louise:0.1);
        tree B = ((Simon:0.1,Harry:0.2):0.1,Betty:0.2):Louise:0.1);
    end;



You can do the following:

Find out how many characters:
>>> n.data.nchar
2

Ask about how many taxa:
>>> n.data.ntaxa
4

Get the taxa names:
>>> n.data.taxa
['Harry', 'Simon', 'Betty', 'Louise']


Get the `format` info:
>>> n.data.format
{'datatype': 'standard', 'symbols': '01', 'gap': '-'}


The actual data matrix is a dictionary, which you can get to in `.matrix`:
>>> n.data.matrix
{'Simon': ['0', '1'], 'Louise': ['1', '1'], 'Betty': ['1', '0'], 'Harry': ['0', '0']}

Or, you could access the data matrix via taxon:
>>> n.data.matrix['Simon']
['0', '1']

Or even loop over it like this:
>>> for taxon, characters in n.data:
>>>     print taxon, characters

You can also iterate over the sites (rather than the taxa):

>>> for site, data in n.data.characters.items():
>>>     print site, data
0 {'Simon': '0', 'Louise': '1', 'Betty': '1', 'Harry': '0'}
1 {'Simon': '1', 'Louise': '1', 'Betty': '0', 'Harry': '0'}

..or you can access the characters matrix directly:

>>> n.data.characters[0]
{'Simon': '0', 'Louise': '1', 'Betty': '1', 'Harry': '0'}

NOTE: that sites are zero-indexed!



`trees` block handler
---------------------

If there's a `trees` block, then you can do the following

You can get the number of trees:
>>> n.trees.ntrees
2

You can access the trees via the `.trees` dictionary:
>>> n.trees.trees[0]
'tree A = ((Harry:0.1,Simon:0.2):0.1,Betty:0.2):Louise:0.1);'

Or loop over them:
>>> for tree in n.trees:
>>>     print tree

>>>n.blocks['trees'].detranslate()

`taxa` block handler
--------------------

These are the alternate nexus file format found in programs like SplitsTree:

    BEGIN Taxa;
    DIMENSIONS ntax=4;
    TAXLABELS
    [1] 'John'
    [2] 'Paul'
    [3] 'George'
    [4] 'Ringo'
    ;
    END; [Taxa]


In a taxa block you can get the number of taxa and the taxa list:

>>> n.taxa.ntaxa
4
>>> n.taxa.taxa
['John', 'Paul', 'George', 'Ringo']

NOTE: with this alternate nexus format the Characters blocks *should* be parsed by
DataHandler.


Writing a Nexus File using NexusWriter
================================================================================


NexusWriter provides more fine-grained control over writing nexus files, and 
is useful if you're programmatically generating a nexus file rather than loading
a pre-existing one.

>>> from nexus import NexusWriter
>>> n = NexusWriter()

Add a comment to appear in the header of the file
>>> n.add_comment("I am a comment")

Data are added by using the "add" function - which takes 3 arguments, a taxon, 
a character name, and a value.

>>> n.add('taxon1', 'Character1', 'A')
>>> n.data
{'Character1': {'taxon1': 'A'}}

>>> n.add('taxon2', 'Character1', 'C')
>>> n.add('taxon3', 'Character1', 'A')

Characters and values can be strings or integers
>>> n.add('taxon1', 2, 1)
>>> n.add('taxon2', 2, 2)
>>> n.add('taxon3', 2, 3)

NexusWriter will interpolate missing entries (i.e. taxon2 in this case)
>>> n.add('taxon1', "Char3", '4')
>>> n.add('taxon3', "Char3", '4')

... when you're ready, you can generate the nexus using `make_nexus` or `write_to_file`:
>>> data = n.make_nexus(interleave=True, charblock=True)

>>> n.write_to_file(filename="output.nex", interleave=True, charblock=True)

... you can make an interleaved nexus by setting `interleave` to True, and you can
include a character block in the nexus (if you have character labels for example) 
by setting charblock to True.


Nexus manipulation scripts included
================================================================================

nexus_combine_nexus.py
----------------------

Combines a series of nexuses into one nexus ('combined.nex').

Usage: python nexus_combine_nexus.py filename1.nex filename2.nex ... filenameN.nex

nexus_deinterleave.py
---------------------

Converts an interleaved nexus to a simple nexus.

Usage: python nexus_deinterleave.py oldnexus.nex newnexus.nex

nexus_detranslate.py
--------------------

Converts a nexus tree file with 'translated' taxa labels to a detranslated form.

Usage: python nexus_detranslate.py oldnexus.trees newnexus.trees

nexus_randomise.py
------------------

Randomly shuffles the character states between taxa in a nexus file.

Usage: python nexus_randomise.py oldnexus.nex randomised.nex


nexus_treemanip.py
------------------

Provides a number of functions for manipulating trees.

Usage: python nexus_treemanip.py [option] oldnexus.trees newnexus.trees

Deleting trees:
    nexus_treemanip.py -d 1 oldnexus.trees newnexus.trees   - delete tree #1
    nexus_treemanip.py -d 1-5 oldnexus.trees newnexus.trees - delete trees #1-5
    nexus_treemanip.py -d 1,5 oldnexus.trees newnexus.trees - delete trees #1 and #5
    nexus_treemanip.py -d 1,4,20-30 oldnexus.trees newnexus.trees - delete trees #1, #4, #20-30
    
Resampling trees:
    nexus_treemanip.py -r 10 oldnexus.trees newnexus.trees   - resample every 10th tree
    
Remove comments:
    nexus_treemanip.py -c oldnexus.trees newnexus.trees

Sampling N Random trees:
    nexus_treemanip.py -n 100 oldnexus.trees newnexus.trees  - randomly sample 100 trees.


nexus_nexusmanip.py
-------------------

Provides a number of functions for manipulating nexus character files.

Usage: python nexus_nexusmanip.py [option] oldnexus.nex [newnexus.nex]

Count missing characters:
    Usage: python nexus_nexusmanip.py -m oldnexus.nex

Remove constant characters:
    Usage: python nexus_nexusmanip.py -c oldnexus.nex newnexus.nex
    
Remove the unique characters:
    Usage: python nexus_nexusmanip.py -u oldnexus.nex newnexus.nex
    

nexinfo.py
----------

Displays quick information about the nexus file loaded.



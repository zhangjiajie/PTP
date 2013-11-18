"""
python-nexus - Generic nexus (.nex, .trees) reader for python
=============================================================

Reading a Nexus
---------------

>>> import os
>>> EXAMPLE_DIR = os.path.join(os.path.dirname(__file__), 'examples')
>>>
>>> from nexus import NexusReader
>>> n = NexusReader()
>>> n.read_file(os.path.join(EXAMPLE_DIR, 'example.nex'))
...
>>> n = NexusReader(os.path.join(EXAMPLE_DIR, 'example.nex'))
...

# display blocks found in data file
>>> n.blocks
{'data': <NexusDataBlock: 2 characters from 4 taxa>}

`data` blocks
-------------

>>> n.data.nchar
2

>>> n.data.ntaxa
4

>>> n.data.format
{'datatype': 'standard', 'symbols': '01', 'gap': '-'}

>>> n.data.matrix
{'Simon': ['0', '1'], 'Louise': ['1', '1'], 'Betty': ['1', '0'], 'Harry': ['0', '0']}

>>> n.data.matrix['Simon']
['0', '1']

>>> sorted(n.data.taxa)
['Betty', 'Harry', 'Louise', 'Simon']

>>> sorted(n.data.matrix.keys())
['Betty', 'Harry', 'Louise', 'Simon']

>>> for taxon, characters in n.data: #doctest: +SKIP


`tree` blocks
-------------

>>> n = NexusReader(os.path.join(EXAMPLE_DIR, 'example.trees'))
>>> n.trees.ntrees
3
>>> n.trees.trees[0]
'tree tree.0.1065.603220 = (((((((Chris:0.0668822155,Bruce:0.0173144449):0.0062091603,Tom:0.0523825242):0.0206190840,(Henry:0.0482653647,Timothy:0.0744964092):0.0183093750):0.0401805957,(Mark:0.0066961591,Simon:0.0755275882):0.0264078188):0.0536464636,((Fred:0.0428499135,Kevin:0.0734738565):0.0937536292,Roger:0.0538708492):0.0438297939):0.0453008384,(Michael:0.0953237112,Andrew:0.0654710419):0.0803079594):0.0630363263,David:0.0855948485);'

>>> for tree in n.trees: #doctest: +SKIP



Writing a Nexus File
--------------------
>>> from nexus import NexusWriter
>>> n = NexusWriter()

Add a comment to appear in the header of the file
>>> n.add_comment("I am a comment")

data are added by using the "add" function - 
which takes 3 arguments, a taxon, a character name, and a value

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
n.make_nexus(interleave=True, charblock=True)
n.write_to_file(filename="output.nex", interleave=True, charblock=True)

"""
__author__ = 'Simon Greenhill <simon@simon.net.nz>'

from reader import *
from writer import NexusWriter

__version__ = "0.87"
PACKAGE_NAME = "python-nexus"
PACKAGE_VERSION = __version__
VERSION = __version__
PACKAGE_AUTHOR = "Simon J. Greenhill"
PACKAGE_COPYRIGHT = "Copyright 2009-2011 Simon J. Greenhill"
PACKAGE_LICENSE = """
Copyright (c) 2009-2011, Simon J. Greenhill
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

    1. Redistributions of source code must retain the above copyright notice, 
       this list of conditions and the following disclaimer.
    
    2. Redistributions in binary form must reproduce the above copyright 
       notice, this list of conditions and the following disclaimer in the
       documentation and/or other materials provided with the distribution.

    3. Neither the name of python-nexus nor the names of its contributors may be
       used to endorse or promote products derived from this software without
       specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""


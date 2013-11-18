"""
Tools for writing a nexus file
"""

TEMPLATE = """
#NEXUS 
%(comments)s
BEGIN DATA;
    DIMENSIONS NTAX=%(ntax)d NCHAR=%(nchar)d;
    FORMAT DATATYPE=%(datatype)s MISSING=%(missing)s GAP=%(gap)s %(interleave)s SYMBOLS="%(symbols)s";
    %(charblock)s
MATRIX
%(matrix)s
;
END;
"""
from reader import NexusReader

class NexusWriter:
    
    MISSING = '?'
    GAP = '-'
    DATATYPE = 'STANDARD'
    
    def __init__(self):
        self.taxalist = []
        self.comments = []
        self.characters = []
        self.clean_characters = {}
        self.symbols = set()
        self.data = {}
        self.is_binary = False
        
    def clean(self, s):
        """Removes unsafe characters"""
        replacements = {' ': '', '\\': '', '(':'_', ')':'', ':': '', '/':'', '?': '', '-': ''}
        for f,t in replacements.items():
            s = s.replace(f, t)
        return s
        
    def _add_char(self, charlabel):
        """Adds a character"""
        charlabel = str(charlabel)
        if charlabel not in self.characters:
            self.characters.append(charlabel)
            self.data[charlabel] = {}
            self.clean_characters[charlabel] = self.clean(charlabel)
        return charlabel
        
    def _add_taxa(self, taxon):
        """Adds a taxa"""
        if taxon not in self.taxalist:
            self.taxalist.append(taxon)
        return taxon
        
    def _make_charlabel_block(self):
        """Generates a character label block"""
        out = ["CHARSTATELABELS"]
        for i, c in enumerate(self.characters, 1):
            out.append("\t\t%d %s," % (i, self.clean_characters[c]))
        out[-1] = out[-1].strip(',') # remove trailing comma
        out.append(";")
        return "\n".join(out)
        
    def _make_matrix_block(self, interleave):
        """Generates a matrix block"""
        out = []
        if interleave:
            for c in self.characters:
                for t in self.taxalist:
                    out.append("%s %s" % (t.ljust(25), self.data[c].get(t, self.MISSING)))
                out.append("")
        else:
            for t in sorted(self.taxalist):
                s = []
                for c in self.characters:
                    value = self.data[c].get(t, self.MISSING)
                    if len(value) > 1: # wrap equivocal states in ()'s
                        value = "(%s)" % value
                    s.append(value)
                out.append("%s %s" % (t.ljust(25), ''.join(s)))
        return "\n".join(out)
        
    def _make_comments(self):
        """Generates a comments block"""
        return "\n".join(["[%s]" % c.ljust(70) for c in self.comments])
        
    def add_comment(self, comment):
        """Adds a `comment` into the nexus file"""
        self.comments.append(comment)
        
    def add(self, taxon, character, value):
        """Adds a `character` for the given `taxon` and sets it to `value`"""
        assert self.is_binary == False, "Unable to add data to a binarised nexus form"
        character = self._add_char(character)
        taxon = self._add_taxa(taxon)
        value = str(value)
        # have multiple entries
        #assert (taxon in self.data[character]) == False, "Duplicate entry for %s-%s" % (taxon, character)
        if taxon in self.data[character]:
            self.data[character][taxon] += value
        else:
            self.data[character][taxon] = value
        
        # add to symbols
        if value not in ['?', '-']:
            [self.symbols.add(v) for v in value]
        
    def recode_to_binary(self):
        """Recodes the matrix to binary form"""
        
    def write(self, interleave=False, charblock=False):
        """
        Generates a string representation of the nexus 
        (basically a wrapper around make_nexus)
        
        :param interleave: Generate interleaved matrix or not
        :type interleave: Boolean 
        :param charblock: Include a characters block or not 
        :type charblock: Boolean 
        
        :return: String
        """
        return self.make_nexus(interleave, charblock)
        
    def make_nexus(self, interleave=False, charblock=False):
        """
        Generates a string representation of the nexus
        
        :param interleave: Generate interleaved matrix or not
        :type interleave: Boolean 
        :param charblock: Include a characters block or not 
        :type charblock: Boolean 
        
        :return: String
        """
        return TEMPLATE.strip() % {
            'ntax': len(self.taxalist),
            'nchar': len(self.characters),
            'charblock': self._make_charlabel_block() if charblock else '',
            'matrix': self._make_matrix_block(interleave=interleave),
            'interleave': 'INTERLEAVE' if interleave else '',
            'comments': self._make_comments(),
            'symbols': ''.join(sorted(self.symbols)),
            'missing': self.MISSING,
            'gap': self.GAP,
            'datatype': self.DATATYPE,
        }
    
    def write_to_file(self, filename="output.nex", interleave=False, charblock=False):
        """
        Generates a string representation of the nexus
        
        :param filename: Filename to store nexus as
        :type filename: String 
        :param interleave: Generate interleaved matrix or not
        :type interleave: Boolean 
        :param charblock: Include a characters block or not 
        :type charblock: Boolean 
        
        :return: None
        """
        
        handle = open(filename, 'w+')
        handle.write(self.make_nexus(interleave, charblock))
        handle.close()
    
    def write_as_table(self):
        """
        Generates a simple table of the nexus
        """
        out = []
        for t in sorted(self.taxalist):
            s = []
            for c in self.characters:
                value = self.data[c].get(t, self.MISSING)
                if len(value) > 1: # wrap equivocal states in ()'s
                    value = "(%s)" % value
                s.append(value)
            out.append("%s %s" % (t.ljust(25), ''.join(s)))
        return "\n".join(out)
    
    def _convert_to_reader(self):
        """
        This is a hack to convert a NexusWriter object to a NexusReader instance
        
        One day I'll refactor this all so Reader and Writer subclass something,
        which will make this unnecessary.
        """
        n = NexusReader()
        out = self.make_nexus(interleave=False, charblock=True)
        n.read_string(out)
        return n
        
        

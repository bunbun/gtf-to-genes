#!/usr/bin/env python
"""

    Minimal class for parsing single entries (from a single line of a GTF file)
    
    Original code from Andreas Heger
    
    Tested with GTF version 2.2. data from Ensembl

    From http://mblab.wustl.edu/GTF22.html:

        GTF stands for Gene transfer format. It borrows from GFF, but has additional 
        structure that warrants a separate definition and format name.
        Structure is as GFF, so the fields are:
        <seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments] 

"""

################################################################################
#
#   minimal_gtf_iterator
#
#
#   Copyright (c) 3/11/2010 Leo Goodstadt
#   
#   Permission is hereby granted, free of charge, to any person obtaining a copy
#   of this software and associated documentation files (the "Software"), to deal
#   in the Software without restriction, including without limitation the rights
#   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#   copies of the Software, and to permit persons to whom the Software is
#   furnished to do so, subject to the following conditions:
#   
#   The above copyright notice and this permission notice shall be included in
#   all copies or substantial portions of the Software.
#   
#   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
#   THE SOFTWARE.
#################################################################################

import sys, os
import os.path


# add self to search path for testing
if __name__ == '__main__':
    exe_path = os.path.split(os.path.abspath(sys.argv[0]))[0]
    # Use import path from <<../python_modules>>
    sys.path.append(os.path.abspath(os.path.join(exe_path,"..", "python_modules")))
    module_name = os.path.split(sys.argv[0])[1]
    module_name = os.path.splitext(module_name)[0];
else:
    module_name = __name__



#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#   options        


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888


if __name__ == '__main__':
    from optparse import OptionParser
    import StringIO
    
    parser = OptionParser(version="%prog 1.0", usage = "\n\n    %progs [options]")
    parser.add_option("-i", "--input_file", dest="input_file",
                      metavar="FILE", 
                      type="string",
                      help="Name and path of input file. "
                          "Defaults to reading from STDIN.")
    
    #
    #   general options: verbosity / logging
    # 
    parser.add_option("-v", "--verbose", dest = "verbose",
                      action="count", default=0,
                      help="Print more verbose messages for each additional verbose level.")
    parser.add_option("-L", "--log_file", dest="log_file",
                      metavar="FILE", 
                      type="string",
                      help="Name and path of log file")
    parser.add_option("--skip_parameter_logging", dest="skip_parameter_logging",
                        action="store_true", default=False,
                        help="Do not print program parameters to log.")
    parser.add_option("--debug", dest="debug",
                        action="count", default=0,
                        help="Set default program parameters in debugging mode.")
    
    
    
    
    # get help string
    f =StringIO.StringIO()
    parser.print_help(f)
    helpstr = f.getvalue()
    
    arguments = " ".join(sys.argv)
    
    (options, remaining_args) = parser.parse_args()
    
    
    #vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
    #                                             #
    #   Debug: Change these                       #
    #                                             #
    #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    if options.debug:
        options.log_file                = os.path.join("minimal_gtf_iterator.log")
        options.verbose                 = 5
        options.log_parameters          = True
    #vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
    #                                             #
    #   Debug: Change these                       #
    #                                             #
    #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    
    
    
#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#   imports        


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888


#from json import dumps
#from collections import defaultdict



#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#   classes


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
#_________________________________________________________________________________________

#   Error

#_________________________________________________________________________________________
class Error(Exception):
    """Base class for exceptions in this module."""
    def __str__(self):
        return str(self.message)
    def _get_message(self, message): return self._message
    def _set_message(self, message): self._message = message
    message = property(_get_message, _set_message)

#_________________________________________________________________________________________

#   ParsingError

#_________________________________________________________________________________________
class ParsingError(Error):
    """Exception raised for errors in the input.

    Attributes:
        message -- explanation of the error
    """

    def __init__(self, message):
        self.message = message

#_________________________________________________________________________________________

#   gtf_entry

#_________________________________________________________________________________________
class gtf_entry:
    """read gtf formatted entry representing one line in a GTF file

    The coordinates are kept internally in python coordinates (0-based, open-closed), but are
    output as inclusive 1-based coordinates according to

    http://www.sanger.ac.uk/Software/formats/GFF/
    """

    #_____________________________________________________________________________________

    #   init

    #_____________________________________________________________________________________
    def __init__(self):
        self.mContig = "."
        self.mSource = "."
        self.mFeature = "."
        self.mFrame = "."
        self.mStart = 0
        self.mEnd = 0
        self.mScore = "."
        self.mStrand = "."
        self.mGeneId = None
        self.mTranscriptId =  None
        self.mAttributes = {}

    #_____________________________________________________________________________________

    #   read

    #_____________________________________________________________________________________
    def read( self, line ):
        """read gff entry from line.
        
        <seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments]
        """
        
        data = line[:-1].split("\t")

        try:
            (self.mContig, self.mSource, self.mFeature,
             self.mStart, self.mEnd, self.mScore, self.mStrand,
             self.mFrame ) = data[:8]
        except ValueError:
            raise ValueError( "parsing error in line `%s`" % line )

        ## note: frame might be .
        (self.mStart, self.mEnd) = map(int, (self.mStart, self.mEnd))
        self.mStart -= 1

        self.parseInfo( data[8], line )

    #_____________________________________________________________________________________

    #   parseInfo

    #_____________________________________________________________________________________
    def parseInfo( self, attributes, line ):
        """parse attributes.
        """
        # remove comments
        attributes = attributes.split( "#" )[0]
        # separate into fields
        fields = map( lambda x: x.strip(), attributes.split(";")[:-1])
        self.mAttributes = {}
        
        for f in fields:
            
            d = map( lambda x: x.strip(), f.split(" "))
            
            n,v = d[0], d[1]
            if len(d) > 2: v = d[1:]

            if v[0] == '"' and v[-1] == '"':
                v = v[1:-1]
            else:
                ## try to convert to a value
                try:
                    v = float( v )
                    v = int( v )
                except ValueError:
                    pass
                except TypeError:
                    pass
                
            if n == "gene_id": 
                self.mGeneId = v     
            elif n == "transcript_id": 
                self.mTranscriptId = v
            else: 
                self.mAttributes[n] = v

        if not self.mGeneId:
            raise ParsingError( "missing attribute 'gene_id' in line %s" % line)
        if not self.mTranscriptId:
            raise ParsingError( "missing attribute 'transcript_id' in line %s" % line)





#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#   Functions


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#_________________________________________________________________________________________

#   iterator

#_________________________________________________________________________________________
def iterator( infile ):
    """return a simple iterator over all entries in a file."""
    while 1:
        line = infile.readline()
        if not line: raise StopIteration
        if line.startswith("#"): continue
        gtf = gtf_entry()
        gtf.read( line )
        yield gtf

#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#   Logger


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

if __name__ == '__main__':
    #
    #   set up log
    # 
    import logging
    logger = logging.getLogger(module_name)
    from default_logger import  setup_default_log, MESSAGE
    setup_default_log(logger, options.log_file, options.verbose)
    

    #
    #   log programme parameters
    # 
    if not options.skip_parameter_logging:
        logger.info("%s" % (arguments))

#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#   Main logic


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
if __name__ == '__main__':
    from itertools import izip
    from json import dumps
    
    #
    #   assumes debug mode for unit testing
    # 
    if options.debug:
        import unittest
        class Test_minimal_gtf_iterator(unittest.TestCase):

            #       self.assertEqual(self.seq, range(10))
            #       self.assert_(element in self.seq)
            #       self.assertRaises(ValueError, random.sample, self.seq, 20)



            def test_function(self):
                """
                    test 
                """
                test_results = []
                for i, gtf_entry in izip(range(2),  iterator(open(os.path.join(exe_path, "test_data", "test.shortish")))):
                    test_results.append(dumps(gtf_entry.__dict__))
                self.assertEqual(test_results, [
                    '{"mSource": "pseudogene", "mFeature": "exon", "mGeneId": "ENSG00000224777", "mScore": ".", "mTranscriptId": "ENST00000424047", "mFrame": ".", "mAttributes": {"exon_number": "1", "transcript_name": "OR4F2P-001", "gene_name": "OR4F2P"}, "mEnd": 87586, "mStart": 86648, "mStrand": "-", "mContig": "11"}',
                    '{"mSource": "protein_coding", "mFeature": "exon", "mGeneId": "ENSG00000230724", "mScore": ".", "mTranscriptId": "ENST00000382784", "mFrame": ".", "mAttributes": {"exon_number": "1", "transcript_name": "AC069287.3-201", "gene_name": "AC069287.3"}, "mEnd": 129388, "mStart": 129059, "mStrand": "-", "mContig": "11"}'])

        #
        #   call unit test without parameters
        #     

        if sys.argv.count("--debug"):
            sys.argv.remove("--debug")
        unittest.main()










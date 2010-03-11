#!/usr/bin/env python
"""

    dump_object.py
    [--log_file PATH]
    [--verbose]

"""

################################################################################
#
#   dump_object
#
#
#   Copyright (c) 3/2/2010 Leo Goodstadt
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
from collections import defaultdict


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#   Functions        


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
#_________________________________________________________________________________________    

#   __str__

#_________________________________________________________________________________________    
import pprint
import re
object_regex = re.compile("<.+ object at 0x.*>")

def dump_object (o, linker_str = "\n"):
    """
    dumps object to string
    """
        
    linker_str += "    "
    attribute_names = [attr for attr in dir(o) \
                       if not callable(getattr(o,attr)) and 
                        not attr.startswith("__")]
    if not len(attribute_names):
        return ""
    max_len = max(len(s) for s in attribute_names)
    format_str = "%-" + str(max_len) + "s = %s"
    data_strs = []
    
    pprinter = pprint.PrettyPrinter(indent= 4)

    for attr_name in sorted(attribute_names):
        if len(attr_name) >= 2 and attr_name[0:2] == "__":
            continue
        d = getattr(o, attr_name)
        indent_str = " " * (max_len + 3)

        #
        #   pprinter first
        # 
        object_d_str = pprinter.pformat(d)
        if not object_regex.match(object_d_str):
            d_str = object_d_str.replace("\n", linker_str + indent_str)
        else:
            try:
                d_str = str(d)
                d_str = d_str.replace("\n", linker_str)
            except:
                try:
                    d_str = linker_str + indent_str + "\n" + dump_object(d, linker_str )[4:]
                except:
                    d_str = object_d_str
        
        data_strs.append(format_str % (attr_name, d_str))
        
    return (linker_str).join(data_strs)
    

#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#   Testing


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
import unittest
class CC(object):
    def __init__(self):
        """
        """
        self.CC_e = 3
        self.CC_f = [2,3,7, "Very long string", "Very long string", ["Very long string", "Very long string", "Very long string", "Very long string"], "Very long string", "Very long string", ]

    def __repr__(self):
        return dump_object(self, "\n")
        
class BB(object):
    def __init__(self):
        """
        """
        self.BB_c = 3
        self.BB_d = [2,3,5, CC()]
        self.BB_CC = CC()

    def __str__(self):
        """
        """
        #data_str, format_str = dump_object(self, "\n", set(["BB_d"]))
        data_str = dump_object(self, "\n")
        
        return data_str
        
        
class AA(object):
    def __init__(self):
        """
        """
        self.AA_a = 3
        self.AA_b = [2,3,6]
        self.AA_bb = BB()

    def __str__(self):
        """
        """
        return "    " + dump_object(self,  "\n")
class Test_dump_object(unittest.TestCase):

    #       self.assertEqual(self.seq, range(10))
    #       self.assert_(element in self.seq)
    #       self.assertRaises(ValueError, random.sample, self.seq, 20)


    
    def test_function(self):
        """
            test 
        """
        a = dict()
        a[1] = 2
        a[45]= 3
        print str(AA())
        import pprint
        pp = pprint.PrettyPrinter(indent=4)
        
        

#
#   debug code not run if called as a module
#     
if __name__ == '__main__':
    if sys.argv.count("--debug"):
        sys.argv.remove("--debug")
    unittest.main()



    

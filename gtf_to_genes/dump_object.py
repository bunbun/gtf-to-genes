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
import array
import textwrap

#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#   Functions


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
__No    = 0
__TRY   = 1
__Yes   = 2
__LIST  = 0
__TUPLE = 1
__SET   = 2
#_________________________________________________________________________________________

#   __str__

#_________________________________________________________________________________________
import pprint
pprinter4 = pprint.PrettyPrinter(indent= 4)
import re
object_regex = re.compile("<.+ object at 0x.*>")

#_________________________________________________________________________________________

#   find_recursion

#_________________________________________________________________________________________
def is_recursive (item, context):
    """

    """
    if id(item) in context:
        return True

    context.add(id(item))

    if item is None:
        return False

    if isinstance(item, dict):
        for k, v in item.iteritems():
            if (is_recursive(k, context) or
                is_recursive(v, context)):
                return True


    elif isinstance(item, (list, set, tuple)):
        for v in item:
            if is_recursive(v, context):
                return True

    elif isinstance(item, (basestring, int, float, long, complex)):
        return False

    else:
        for attr in dir(item):
            v= getattr(item, attr)
            if not callable(v) and is_recursive(v, context):
                return True

    return False



#_________________________________________________________________________________________

#   do_get_item_str

#_________________________________________________________________________________________
def do_get_item_str (item, context, use_new_line, line_length, max_levels, sorted_attribute_names = [], ignored_attribute_names =set(), prevent_recursive_str_call = False):

    nested_levels = 0
    if item is None:
        d_str = repr(item)
    elif isinstance(item, dict):
        #print >>sys.stderr, " ^^ dump dict"
        d_str, nested_levels = do_dump_object_or_dict(item, context, use_new_line, line_length, max_levels, sorted_attribute_names, ignored_attribute_names)
    elif isinstance(item, list):
        #print >>sys.stderr, " ^^ dump list"
        d_str, nested_levels = do_dump_list_or_tuple(item, context, use_new_line, __LIST, line_length, max_levels)
    elif isinstance(item, set):
        #print >>sys.stderr, " ^^ dump set"
        d_str, nested_levels = do_dump_list_or_tuple(sorted(item), context, use_new_line, __SET, line_length, max_levels)
    elif isinstance(item, tuple):
        #print >>sys.stderr, " ^^ dump tuple"
        d_str, nested_levels = do_dump_list_or_tuple(item, context, use_new_line, __TUPLE, line_length, max_levels)
    elif isinstance(item, float):
        d_str = str(item)
    elif isinstance(item, array.array):
        #d_str, nested_levels = do_dump_list_or_tuple(list(item), context, __No, __LIST, line_length, max_levels)
        d_str = textwrap.fill(repr(item),  line_length)
    elif isinstance(item, (basestring, int, float, long, complex)):
        d_str = textwrap.fill(repr(item),  line_length)
    elif prevent_recursive_str_call and is_recursive (item, set()):
        #print >> sys.stderr, "prevent_recursive_str_call", item.__class__
        d_str, nested_levels = do_dump_object_or_dict(item, context, use_new_line, line_length, max_levels, sorted_attribute_names, ignored_attribute_names)
    else:
        try:
            r = getattr(type(item), "__str__", None)
            # if you have not redefined __str__, I will go down the rabbit hole
            if isinstance(item, object) and r is object.__str__:
                #print >> sys.stderr, "Original __str__", item.__class__
                d_str, nested_levels = do_dump_object_or_dict(item, context, use_new_line, line_length, max_levels, sorted_attribute_names, ignored_attribute_names)
            else:
                #print >> sys.stderr, "Redefined __str__", item.__class__
                d_str = str(item)
                if object_regex.match(d_str):
                    #print >> sys.stderr, "Rubbish __str__", item.__class__
                    d_str, nested_levels = do_dump_object_or_dict(item, context, use_new_line, line_length, max_levels, sorted_attribute_names, ignored_attribute_names)
        except:
            raise

    if use_new_line == __No and len(d_str) > line_length:
        if "\n" in d_str:
            # We have use new line already: ignore
            # something via the string or pformat ended up using new lines
            return d_str, nested_levels

        # OK try the same but with new lines trying to shorten the line
        #print >>sys.stderr, " !! %5s >>%s<< length=%d, line=%d" % (use_new_line, d_str, len(d_str), line_length)
        return do_get_item_str(item, context, __TRY, line_length, max_levels, sorted_attribute_names, ignored_attribute_names, prevent_recursive_str_call)



    #print >>sys.stderr, " ?? %5s >>%s<< length=%d, line=%d" % (use_new_line, d_str, len(d_str), line_length)
    return d_str, nested_levels


#_________________________________________________________________________________________

#   do_dump_list_or_tuple

#_________________________________________________________________________________________
def do_dump_list_or_tuple (o, context, use_new_line, islist, line_length, max_levels):

    if id(o) in context:
        #return "'Back reference to <%s> object'" % (o.__class__.__name__), use_new_line, 0
        return "'Back reference to <%s> object'" % (o.__class__.__name__), 0
    context.add(id(o))


    if islist == __LIST:
        brackets1 = '['
        brackets2 = ']'
        end_comma = ''
    elif islist == __TUPLE:
        brackets1 = '('
        brackets2 = ')'
        end_comma = ','
    elif islist == __SET:
        brackets1 = 'set(['
        brackets2 = '])'
        end_comma = ''
    else:
        raise Exception("islist (%s) not understood" % (islist))

    data_strs = []

    indent4 = " " * 4
    line_length -= 4
    nested_levels = 0

    #print >>sys.stderr, " --", "%5s" % use_new_line, ">>%s<<" % str(o), use_new_line, 0

    if use_new_line == __TRY:
        sub_items_use_new_line = __No
    else:
        sub_items_use_new_line = use_new_line
    max_nested_levels = 0
    any_use_new_line = use_new_line
    data_strs = []
    for item in o:
        item_str, nested_levels = do_get_item_str (item, context, sub_items_use_new_line, line_length, max_levels)
        item_use_new_line = "\n" in item_str
        max_nested_levels = max(max_nested_levels, nested_levels)
        any_use_new_line = max(any_use_new_line, item_use_new_line)
        data_strs.append(item_str)

    if max_nested_levels >= max_levels:
        any_use_new_line = __Yes

    if any_use_new_line in (__Yes, __TRY):
        for i, s in enumerate(data_strs):
            data_strs[i] = indent4 + s.replace("\n", "\n" + indent4)

    #print >>sys.stderr, " ++", "%5s" % use_new_line, ">>%s<<" % str(o), any_use_new_line, max_nested_levels

    if len(data_strs) >= 2:
        end_comma = ''

    context.remove(id(o))
    if any_use_new_line in (__Yes, __TRY):
        ret_str = (brackets1 + "\n" +  ",\n".join(data_strs) + end_comma + "\n" + brackets2)
        return ret_str, max_nested_levels + 1
    else:
        ret_str = brackets1 +  ", ".join(data_strs) + end_comma + brackets2
        if len(ret_str) > line_length + 4:
            return do_dump_list_or_tuple (o, context, __TRY, islist, line_length + 4, max_levels)
        return ret_str, max_nested_levels + 1


#_________________________________________________________________________________________

#   get_keys
#   get_item
#
#       Helper functions so that we can pretend dicts are objects and vice versa

#_________________________________________________________________________________________
def get_keys (o):
    """
    Helper function to get list of keys for dicts or attribute names for objects
    """
    if isinstance(o, dict):
        return o.keys()

    return [attr for attr in dir(o)
                if hasattr(o, attr) and not callable(getattr(o,attr)) and not attr.startswith("__")]

def get_slots (o):
    """
    Helper function to get list of keys for dicts or attribute names for objects
    """
    if isinstance(o, dict):
        return []

    if "__slots__" in dir(o) and o.__slots__ != None:
        return o.__slots__

    return []

def get_item (o, key):
    """
    Helper function to get item given a key (for dict) or attribute name (for objects)
    """
    if isinstance(o, dict):
        return o[key], key

    attr = getattr(o, key)

    if callable(attr):
        if key.startswith("get_"):
            key = key[4:]
        return attr(), key
    else:
        return attr, key

#_________________________________________________________________________________________

#   do_dump_object

#_________________________________________________________________________________________
def do_dump_object_or_dict (o, context, use_new_line, line_length, max_levels, sorted_attribute_names, ignored_attribute_names):

    if id(o) in context:
        #return "'Back reference to <%s> object'" % (o.__class__.__name__), use_new_line, 0
        return "'Back reference to <%s> object'" % (o.__class__.__name__), 0
    context.add(id(o))

    if not len(sorted_attribute_names):
        sorted_attribute_names = get_slots(o)
    #
    #   list any sorted keys before adding the rest in sorted order
    #
    keys = list(sorted_attribute_names)
    attribute_names = get_keys (o)
    if not len(attribute_names):
        return "{}", 0

    #print >>sys.stderr,  "ignored_attribute_names=", ignored_attribute_names
    for key in sorted(attribute_names):
        if key not in keys and key not in ignored_attribute_names:
            keys.append(key)


    #print >>sys.stderr, "!!!keys = ", keys

    data_strs = []

    max_len = max(len(repr(s)) for s in keys)
    format_str = "%-" + str(max_len) + "s : %s"

    line_length -= 4 + max_len + 3
    indent = " " * (4 + max_len + 3)
    indent4 = " " * 4
    nested_levels = 0

    #print >>sys.stderr, " --", "%5s" % use_new_line, ">>%s<<" % str(o), use_new_line, 0

    if use_new_line == __TRY:
        sub_items_use_new_line = __No
    else:
        sub_items_use_new_line = use_new_line
    max_nested_levels = 0
    any_use_new_line = use_new_line
    data_strs = []


    for key in keys:
        item, key_name = get_item(o, key)
        item_str, nested_levels = do_get_item_str (item, context, sub_items_use_new_line, line_length, max_levels)
        item_use_new_line = __Yes if "\n" in item_str else 0
        any_use_new_line = max(any_use_new_line, item_use_new_line)
        max_nested_levels = max(max_nested_levels, nested_levels)
        data_strs.append((repr(key_name), item_str))

    if max_nested_levels >= max_levels:
        any_use_new_line = __Yes

    if any_use_new_line in (__Yes, __TRY):
        for i in range(len(data_strs)):
            data_strs[i] = indent4 + (format_str % (data_strs[i])).replace("\n", "\n" + indent)
    else:
        for i in range(len(data_strs)):
            data_strs[i] = "%s: %s" % (data_strs[i])


    context.remove(id(o))

    if any_use_new_line in (__Yes, __TRY):
        ret_str = "{\n" +  ",\n".join(data_strs) + "\n}"
        return ret_str, max_nested_levels + 1
    else:
        ret_str = '{' +  ", ".join(data_strs) + '}'
        if len(ret_str) > line_length + 4 + max_len + 5:
            return do_dump_object_or_dict (o, context, __TRY, line_length + 4 + max_len + 5, max_levels, sorted_attribute_names, ignored_attribute_names)
        return ret_str, max_nested_levels + 1


#_________________________________________________________________________________________

#   dump_object

#_________________________________________________________________________________________
def dump_object (item, line_length = 80, max_levels = 4, sorted_attribute_names = [], ignored_attribute_names =set()):
    """
    Dump anything to string recursively
    """
    use_new_line = __No
    if max_levels == 0:
        use_new_line = __Yes
    context = set()
    return do_get_item_str (item, context, use_new_line, line_length,
                            max_levels, sorted_attribute_names, ignored_attribute_names, prevent_recursive_str_call = True)[0]
#_________________________________________________________________________________________

#   dump_list
#   dump_tuple

#_________________________________________________________________________________________
#def dump_list (o, line_length = 45, max_levels = 4):
#    """
#    dumps list to string
#    """
#    use_new_line = __No
#    if max_levels == 0:
#        use_new_line = __Yes
#
#    context = set()
#    return do_dump_list_or_tuple (o, context, use_new_line, True, line_length, max_levels - 1)[0]
#
#def dump_tuple (o, line_length = 45, max_levels = 4):
#    """
#    dumps list to string
#    """
#    use_new_line = __No
#    if max_levels == 0:
#        use_new_line = __Yes
#    context = set()
#    return do_dump_list_or_tuple (o, context, use_new_line, False, line_length, max_levels - 1)[0]
#
#def dump_dict (o, line_length = 45, max_levels = 4, sorted_attribute_names = []):
#    """
#    dumps list to string
#    """
#    use_new_line = __No
#    if max_levels == 0:
#        use_new_line = __Yes
#    context = set()
#    return do_dump_object_or_dict (o, context, use_new_line, line_length, max_levels - 1, sorted_attribute_names)[0]
#
#def dump_object (o, line_length = 80, max_levels = 4, sorted_attribute_names = []):
#    """
#    dumps list to string
#    """
#    use_new_line = __No
#    if max_levels == 0:
#        use_new_line = __Yes
#    context = set()
#    return do_dump_object_or_dict (o, context, use_new_line, line_length, max_levels - 1, sorted_attribute_names)[0]




#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#   Testing


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
import unittest
class CC(object):
    def __init__(self):
        """
        """
        self.CC_e = 3
        self.CC_f = [2,3,7, "Very long string", "Very long string", range(10), ["Very long string", "Very long string", "Very long string", "Very long string"], "Very long string", "Very long string", ]

    def __repr__(self):
        return dump_object(self)

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
        data_str = dump_object(self)

        return data_str


class AA(object):
    def __init__(self):
        """
        """
        self.AA_a = 3
        self.AA_b = [2,3,6]
        self.AA_bb = BB()
        self.s = self

    def __str__(self):
        """
        """
        return dump_object(self)
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
        print >>sys.stderr, str(BB())
        #print >>sys.stderr, dump_list ([1,2,[3,[5,[7,range(8), 8], 4]],3, range(12)])
        #print >>sys.stderr, dump_tuple ((1,2,(3,(5,(7,tuple(range(8)), 8), 4)),3, tuple(range(12))))
        #print >>sys.stderr, dump_dict (dict(zip(range(26), "abcdefghijklmnop")))
        print >>sys.stderr, dump_object (dict(zip("abcdefghijklmnopqrstuvwxyz",range(5))))
        d = {
                ('a', 0)  : 4,
                ('b', 1)  : 5,
                ('c', 2)  : 6,
                ('d', 3)  : range(8),
                ('e', 4)  : 8,
                ('f', 5)  : ([1,2,[3,[5,[7,range(5), 8], 4]],3, range(12)]),
                ('g', 6)  : 10,
                ('h', 7)  : 11,
                ('i', 8)  : 12,
                ('j', 9)  : 13,
                ('k', 10) : 14}

        #print >>sys.stderr, dump_dict (d, 60, 0)
        print >>sys.stderr, str(AA())



#
#   debug code not run if called as a module
#
if __name__ == '__main__':
    if sys.argv.count("--debug"):
        sys.argv.remove("--debug")
    unittest.main()



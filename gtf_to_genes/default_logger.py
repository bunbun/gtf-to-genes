#!/usr/bin/env python
"""
Standard way of logging

    logger.debug(...) goes to stderr if verbose is set
    logger.info(...)  goes to the log file if one is set up
    logger.log(MESSAGE,...) goes to both
"""
################################################################################
#
#   default_logger
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
import sys

import logging
import logging.handlers

MESSAGE = 15
logging.addLevelName(MESSAGE, "MESSAGE")

def setup_default_log (logger, log_file, verbose):
    """
    set up log
    """
    class debug_filter(logging.Filter):
        """
        Ignore INFO messages
        """
        def filter(self, record):
            return logging.INFO != record.levelno

    class NullHandler(logging.Handler):
        """
        for when there is no logging    
        """
        def emit(self, record):
            pass

    # We are interesting in all messages
    logger.setLevel(logging.DEBUG)
    has_handler = False

    # log to file if that is specified
    if log_file:
        handler = logging.FileHandler(log_file, delay=False)
        handler.setFormatter(logging.Formatter("%(asctime)s - %(name)s - %(levelname)6s - %(message)s"))
        handler.setLevel(MESSAGE)
        logger.addHandler(handler)
        has_handler = True

    # log to stderr if verbose
    if verbose:
        stderrhandler = logging.StreamHandler(sys.stderr)
        stderrhandler.setFormatter(logging.Formatter("    %(message)s"))
        stderrhandler.setLevel(logging.DEBUG)
        if log_file:
            stderrhandler.addFilter(debug_filter())
        logger.addHandler(stderrhandler)
        has_handler = True

    # no logging
    if not has_handler:
        logger.addHandler(NullHandler())


        
#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#   Testing


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888


#
#   debug code not run if called as a module
#     
if __name__ == '__main__':
    import unittest, os
    #       self.assertEqual(self.seq, range(10))
    #       self.assert_(element in self.seq)
    #       self.assertRaises(ValueError, random.sample, self.seq, 20)
    class Test_adjacent_iter(unittest.TestCase):
        def setUp (self):
            if not os.path.exists("temp_logging"):
                os.makedirs("temp_logging")

        def test_function(self):
            """
                test 
            """
            #
            #   file output is indiscriminate
            # 
            #   debug output depends on verbose
            # 
            # 
            
            
            loggers = []
            for i in range(4):
                loggers.append(logging.getLogger("test%d" % i))
            setup_default_log (loggers[0], "temp_logging/test1.log", 1)
            setup_default_log (loggers[1], None, 0)
            setup_default_log (loggers[2], None, 1)
            setup_default_log (loggers[3], "temp_logging/test2.log", 0)
            for i in range(4):
                loggers[i].error("logger %d error test" % i)
                loggers[i].debug("logger %d debug test" % i)
                loggers[i].log(MESSAGE, "logger %d message test" % i)


    if sys.argv.count("--debug"):
        sys.argv.remove("--debug")
    unittest.main()





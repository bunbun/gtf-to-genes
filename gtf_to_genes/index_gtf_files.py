#!/usr/bin/env python
"""
    index_gtf_files.py
    
===========================    
    To index all gtf files    
===========================
            
        ::

            index_file           = "/directory/where/index/lives"
            search_path_root     = "/databases/ftp.ensembl.org/"
            regex_input          = "(.+\/)(([^.]+)\..+\.(.+)\.gtf(?:\.gz)?)$"
            regex_output_pattern = r"\1\2.cache"
                                   #or  r"{INDEX_FILE_PATH}/\2.cache"
            identifier1_pattern  = r"\3"
            identifier1_pattern  = r"\4"

            index_gtf_files(index_file, 
                            search_path_root, 
                            regex_input, 
                            regex_output_pattern,
                            identifier1_pattern, 
                            identifier2_pattern)
=================================================    
    To read indexed genes from a particular species
=================================================
            
        ::
    
            >> species, version, genes = get_indexed_genes(options.index_file,  logger,  "Macaca_mulatta")
            >> print id1, id2
            >> if genes:
            >>     print genes.keys()
            >>     print "# of protein coding genes = ", len(genes['protein_coding'])

            Macaca_mulatta 56
            ['pseudogene', 'snRNA', 'Mt_tRNA', 'protein_coding', 'Mt_rRNA', 'rRNA', 'miRNA', 'misc_RNA', 'snoRNA']
            # of protein coding genes = 21905

==============================================================
    To read indexed genes matching a regular expression
==============================================================

    ::
        
        >> data = get_indexed_genes_from_regex_matching_gtf(options.index_file, logger, "Ciona")
        >> for id1, id2, file_name, genes in data:
        >>     print os.path.basename(file_name)
        >>     print id1, id2, 
        >>     print genes.keys()
        >>     print "# of protein coding genes = ", len(genes['protein_coding'])
        
        Ciona_savignyi.CSAV2.0.56.gtf.gz
        Ciona_savignyi 56 ['pseudogene', 'snRNA', 'protein_coding', 'rRNA', 'miRNA', 'misc_RNA', 'snoRNA']
        # of protein coding genes = 11604
        Ciona_intestinalis.JGI2.56.gtf.gz
        Ciona_intestinalis 56 ['rRNA', 'snRNA', 'protein_coding', 'miRNA', 'misc_RNA', 'snoRNA']
        # of protein coding genes = 14180
        
        
"""

################################################################################
#
#   index_gtf_files
#
#
#   Copyright (c) 3/12/2010 Leo Goodstadt
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
import re
import gene

# add self to search path for testing
if __name__ == '__main__':
    exe_path = os.path.split(os.path.abspath(sys.argv[0]))[0]
    module_name = os.path.split(sys.argv[0])[1]
    module_name = os.path.splitext(module_name)[0];
else:
    module_name = __name__



#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#   options        


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888


if __name__ == '__main__':
    def check_mandatory_options (options, mandatory_options, helpstr):
        """
        Check if specified mandatory options have b een defined
        """
        missing_options = []
        for o in mandatory_options:
            if not getattr(options, o):
                missing_options.append("--" + o)

        if not len(missing_options):
            return

        raise Exception("Missing mandatory parameter%s: %s.\n\n%s\n\n" % 
                        ("s" if len(missing_options) > 1 else "",
                         ", ".join(missing_options),
                         helpstr))

        
    from optparse import OptionParser
    import StringIO
    
    parser = OptionParser(version="%prog 1.0", usage = "\n\n    %progs [options]")
    parser.add_option("-r", "--search_path_root", dest="search_path_root",
                      metavar="FILE", 
                      type="string",
                      help="Root of directory tree which contain GTF files.")
    parser.add_option("-i", "--index_file", dest="index_file",
                      metavar="FILE", 
                      type="string",
                      help="Name and path of index file. "
                          "Index from release/species to gtf file name.")
    parser.add_option("--ignore_cache", dest="ignore_cache",
                        action="store_true", default=False,
                        help="Re-cache genes from GTF.")
    parser.add_option("--regex_input", dest="regex_input",
                        metavar="REGEX", 
                        type="string",
                        default = "(.+\/)(([^.]+)\..+\.(.+)\.gtf(?:\.gz)?)$",
                        help="Regular expression to recognise GTF files and "
                             "to construct the file name of GTF "
                             "cache and two identifiers from the path of the GTF file. "
                             "Defaults to matching the file name part of the path.")
    parser.add_option("--regex_output_pattern", dest="regex_output_pattern",
                        metavar="STRING", 
                        type="string",
                        default = r"\1\2.cache",
                        #default = r"{INDEX_FILE_PATH}/\2.cache",
                        help="Regular expression used to construct file name of GTF cache. "
                             "Defaults to adding a '.cache' suffix")
    parser.add_option("--identifier1_pattern", dest="identifier1_pattern",
                        metavar="STRING", 
                        type="string",
                        default = r"\3",
                        help="Regular expression used to construct the first identifier "
                             "in the GTF file index.")
    parser.add_option("--identifier2_pattern", dest="identifier2_pattern",
                        metavar="STRING", 
                        type="string",
                        default = r"\4",
                        help="Regular expression used to construct the second identifier "
                             "in the GTF file index.")
    
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
    (options, remaining_args) = parser.parse_args()
    
    if options.debug:
        if not options.verbose:
            options.verbose = 1
        if not options.search_path_root:
            options.search_path_root = "/net/cpp-mirror/databases/ftp.ensembl.org"
        if not options.index_file:
            options.index_file = os.path.join(options.search_path_root, "gtf.index")
        if not options.log_file:
            options.log_file = os.path.splitext(options.index_file)[0] + ".log"
            if options.log_file == options.index_file:
                index_file_path = os.path.split(os.path.abspath(options.index_file))[0]
                options.log_file = os.path.join(index_file_path, "cache_genes_from_gtf.log")
                if options.log_file == options.index_file:
                    options.log_file = None
    
    # 
    #   mandatory options
    # 
    mandatory_options = ["log_file", "index_file", "search_path_root"]
    check_mandatory_options (options, mandatory_options, helpstr)
    
    
#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#   imports        


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888


#from json import dumps
#from collections import defaultdict



#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#   Functions        


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#_________________________________________________________________________________________
# 
#    get_indexed_genes
#_________________________________________________________________________________________
def get_indexed_genes(index_file_name, logger, identifier1, identifier2 = None):
    """
    Get gene structures contained in a GTF file whose file name matches identifier1 and identifier2
    
    :param index_file_name: path to index file created by ``index_gtf_files(...)``
    :type index_file_name: string
    :param identifier1: identifier1 parsed from the GTF file name by ``index_gtf_files(...)``
    :type identifier1: string
    :param identifier2: identifier2 parsed from the GTF file name by ``index_gtf_files(...)``
                        If this is None, then only identifier1 will be used to match
    :type identifier2: string
    :rtype tuple of (<matching identifier1>, <identifier2>, <dictionary of lists of genes>)
    """
    index_data = sorted(read_index_file(index_file_name),  reverse = True)
    
    # go through in reverse order of id2 so that higher version numbers are retrieved first
    for id1, id2, original_gtf_path, gtf_cache_path in index_data:
        if id1 == identifier1:
            if identifier2 == None or id2 == identifier2:
                gene_structures = gene.t_parse_gtf(id1 + "." + id2)
                logger.info("Get indexed genes for %s.%s from %s" % 
                                    (id1, id2, original_gtf_path))
                return (id1, id2, gene_structures.get_genes(original_gtf_path, gtf_cache_path, 
                                                            logger))
    logger.info("Identifier %s.%s was not found in the index file %s" % 
                        (identifier1, identifier2, index_file_name))
    return (None,None,None)
    
#_________________________________________________________________________________________
# 
#    get_indexed_genes_from_regex_matching_gtf
#_________________________________________________________________________________________
def get_indexed_genes_from_regex_matching_gtf(index_file_name, logger, regex_str):
    """
    Get gene structures contained in a GTF file whose file name matches regex
    because more than one species may match, returns a list of file names / genes

    :param index_file_name: path to index file created by ``index_gtf_files(...)``
    :type index_file_name: string
    :param regex_str: regular expression used to match GTF file name by ``index_gtf_files(...)``
    :type regex_str: string 
    :rtype list of tuple of (<matching identifier1>, <identifier2>, <dictionary of lists of genes>)
    """
    index_data = sorted(read_index_file(index_file_name),  reverse = True)
    file_name_regex = re.compile(regex_str)

    # go through in reverse order of id2 so that higher version numbers are retrieved first
    results = []
    for id1, id2, original_gtf_path, gtf_cache_path in index_data:
        m = file_name_regex.search(original_gtf_path)
        if not m:
            continue
        gene_structures = gene.t_parse_gtf(id1 + "." + id2)
        logger.info("Get indexed genes for %s.%s from %s" % 
                            (id1, id2, original_gtf_path))
        results.append((id1, id2, original_gtf_path, 
                        gene_structures.get_genes(original_gtf_path, gtf_cache_path, 
                                                  logger)))
    if not len(results):
        logger.info("Regular expression %s did not match any entries in the index file %s" % 
                            (regex_str, index_file_name))
    return results
#_________________________________________________________________________________________
# 
#   read_index_file
#_________________________________________________________________________________________
def read_index_file (index_file_name):
    """
    Reads the index file created by ``index_gtf_files(...)``
    
    :param index_file_name: Path to the index file
    :type index_file_name: string
    :rtype List of [<id1><id2><original_gtf_path><gtf_cache_path>]
    
    """
    index_file = open(index_file_name)
    
    fields_list = []
    for line in index_file:
        if not len(line) or line[0] == '#':
            continue
        fields = line.rstrip().split("\t")
        if len(fields) != 4:
            continue
        fields_list.append(fields)
    return fields_list
    
#_________________________________________________________________________________________
# 
#    index_gtf_files
#_________________________________________________________________________________________

def index_gtf_files(index_file_name, search_path_root, regex_input, regex_output_pattern,
                    identifier1_pattern, identifier2_pattern):
    """
    Iterate through a directory, looking for all GTF files matching a regular expression.
    Cache the GTF data to a file and write the location of the file to an index file
    """
    
    index_file = open(index_file_name,  "w")
    index_file_path = os.path.dirname(index_file_name)
    regex = re.compile(regex_input)
    for root, dirs, files in os.walk(options.search_path_root):
        for file_name in files:
            file_path = os.path.abspath(os.path.join(root, file_name))
            m = regex.search(file_path)
            if not m:
                continue
            cache_file_path  = regex.sub(regex_output_pattern.format(INDEX_FILE_PATH = index_file_path), 
                                        file_path)
            identifier1      = regex.sub(identifier1_pattern , file_path) 
            identifier2      = regex.sub(identifier2_pattern , file_path) 
            index_file.write("%s\t%s\t%s\t%s\n" % 
                                (identifier1, identifier2, file_path, cache_file_path))

            gene_structures = gene.t_parse_gtf(identifier1 + "." + identifier2)
            logger.debug("%s\t%s\t%s\t%s" % 
                                (identifier1, identifier2, file_path, cache_file_path))
            genes_by_type = gene_structures.get_genes(file_path, cache_file_path, logger, ignore_cache = options.ignore_cache)


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
        programme_name = os.path.split(sys.argv[0])[1]
        logger.info("%s" % (" ".join(sys.argv)))

#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#   Functions


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#   Main logic


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

    
    
if __name__ == '__main__':
    
    index_gtf_files(    options.index_file, 
                        options.search_path_root, 
                        options.regex_input, 
                        options.regex_output_pattern,
                        options.identifier1_pattern, 
                        options.identifier2_pattern)

    

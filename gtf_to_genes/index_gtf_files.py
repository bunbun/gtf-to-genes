#!/usr/bin/env python
"""
    index_gtf_files.py

===========================
    To index all gtf files
===========================

        ::
            # index_file = "/net/cpp-mirror/databases/ftp.ensembl.org/gtf.index"

            from gtf_to_genes import *

            index_file           = "/directory/where/index/lives"
            search_path_root     = "/databases/ftp.ensembl.org/"
            regex_input          = "(.+\/)(([^.]+)\..+\.(.+)\.gtf(?:\.gz)?)$"
            cache_file_pattern   = r"\1\2.cache"
                                   #or  r"{INDEX_FILE_PATH}/\2.cache"
            identifier_pattern   = r"\3:\4"

            index_gtf_files(index_file,
                            search_path_root,
                            regex_input,
                            cache_file_pattern,
                            identifier_pattern)
=================================================
    To read indexed genes from a particular species
=================================================

        ::

            >> index_file = "/net/cpp-mirror/databases/ftp.ensembl.org/gtf.index"
            >> from gtf_to_genes import *
            >> import logging
            >> logger = logging.getLogger("test")
            >>
            >> genes_with_identifiers = get_indexed_genes_matching_identifier(index_file,  logger,  "Ciona")
            >> for id, file_path, genes in genes_with_identifiers:
            >>     print id
            >>     print file_path
            >>     print genes.keys()
            >>     print "# of protein coding genes = ", len(genes['protein_coding']), "\n"

            Ciona_intestinalis:56
            /net/cpp-mirror/databases/ftp.ensembl.org/pub/release-56/gtf/ciona_intestinalis/Ciona_intestinalis.JGI2.56.gtf.gz
            ['rRNA', 'snRNA', 'protein_coding', 'miRNA', 'misc_RNA', 'snoRNA']
            # of protein coding genes =  14180

            Ciona_savignyi:56
            /net/cpp-mirror/databases/ftp.ensembl.org/pub/release-56/gtf/ciona_savignyi/Ciona_savignyi.CSAV2.0.56.gtf.gz
            ['pseudogene', 'snRNA', 'protein_coding', 'rRNA', 'miRNA', 'misc_RNA', 'snoRNA']
            # of protein coding genes =  11604

            >> index_file = "/net/cpp-mirror/databases/ftp.ensembl.org/gtf.index"
            >> from gtf_to_genes import *
            >> import logging
            >> logger = logging.getLogger("test")
            >>
            >> species_id, gtf_path, genes = get_indexed_genes_for_identifier(index_file,  logger,  "Ciona_savignyi:56")
            >> print species
            >> if genes:
            >>     print genes.keys()
            >>     print "# of protein coding genes = ", len(genes['protein_coding'])

            Ciona_savignyi:56
            ['pseudogene', 'snRNA', 'protein_coding', 'rRNA', 'miRNA', 'misc_RNA', 'snoRNA']
            # of protein coding genes =  11604

==============================================================
    To read indexed genes matching a regular expression
==============================================================

        ::

            >> index_file = "/net/cpp-mirror/databases/ftp.ensembl.org/gtf.index"
            >> from gtf_to_genes import *
            >> import logging
            >> logger = logging.getLogger("test")
            >>
            >> data = get_indexed_genes_matching_gtf_file_name(index_file, logger, "Ciona")
            >> for id, file_name, genes in data:
            >>     print os.path.basename(file_name)
            >>     print id, genes.keys()
            >>     print "# of protein coding genes = ", len(genes['protein_coding'])

            Ciona_intestinalis.JGI2.56.gtf.gz
            Ciona_intestinalis:56 ['rRNA', 'snRNA', 'protein_coding', 'miRNA', 'misc_RNA', 'snoRNA']
            # of protein coding genes =  14180

            Ciona_savignyi.CSAV2.0.56.gtf.gz
            Ciona_savignyi:56 ['pseudogene', 'snRNA', 'protein_coding', 'rRNA', 'miRNA', 'misc_RNA', 'snoRNA']
            # of protein coding genes =  11604

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
    import argparse
    parser = argparse.ArgumentParser(description='Parses and indices gtf files for fast future '
                                     'retrieval?')

    #
    #   Multiple input files
    #
    # at least one argument, otherwise nargs='*',
    parser.add_argument("-r", "--search_path_root",
                        required = True,
                        metavar="FILE",
                        help="Root of directory tree which contain GTF files.")

    parser.add_argument("-i", "--index_file",
                        required = True,
                        metavar="FILE",
                        help="Name and path of index file. "
                             "Index from release/species to gtf file name.")
    parser.add_argument("--ignore_cache",
                        action="store_true",
                        help="Re-cache genes from GTF.")
    parser.add_argument("--regex_input",
                        metavar="REGEX",
                                                # Ensembl naming scheme
                        default = "((?:.+\/)?)" # path
                                  "(([^.]+)\."  # species
                                  ".+\."        # assembly name
                                  "(.+)"        # version name
                                  "\.gtf(?:\.gz)?)$",
                        help="Regular expression to recognise GTF files and "
                             "to construct the file name of GTF "
                             "cache and two identifiers from the path of the GTF file. "
                             "Defaults to matching the file name part of the path.")
    parser.add_argument("--cache_file_pattern",
                        default = r"\1\2.cache",
                        #default = r"{INDEX_FILE_PATH}/\2.cache",
                        help="Regular expression used to construct file name of GTF cache. "
                             "Defaults to adding a '.cache' suffix")
    parser.add_argument("--identifier_pattern", dest="identifier_pattern",
                        default = r"\3:\4",
                        help="Regular expression used to construct the name identifying the GTF cached data "
                             "in the GTF file index.")

    common_group = parser.add_argument_group('Common arguments')
    common_group.add_argument('--verbose', "-v", const=1,
                              metavar="VERBOSITY", default=0, nargs='?', type= int,
                              help="Print more verbose messages for each additional verbose level.")
    common_group.add_argument('--version', action='version', version='%(prog)s 1.0')
    common_group.add_argument("-L", "--log_file", metavar="FILE", type=str,
                                  help="Name and path of log file")

    options = parser.parse_args()

    if not options.log_file:
        options.log_file            = os.path.join("index_gtf_files.log")



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
def get_indexed_genes_for_identifier(index_file_name, logger, identifier):
    """
    Get gene structures contained in a GTF file whose file name matches identifier

    :param index_file_name: path to index file created by ``index_gtf_files(...)``
    :type index_file_name: string
    :param identifier: identifier parsed from the GTF file name by ``index_gtf_files(...)``
    :type identifier: string
    :rtype tuple of (<matching identifier>, <original GTF path>, <dictionary of lists of genes>)
    """
    index_data = sorted(_read_index_file(index_file_name),  reverse = True)
    all_ids = []
    # go through in reverse order of id2 so that higher version numbers are retrieved first
    for id, original_gtf_path, gtf_cache_path in index_data:
        all_ids.append(id)
        if id == identifier or id == identifier.replace(' ', '_'):
            gene_structures = gene.t_parse_gtf(id)
            logger.info("Get indexed genes for %s from %s" % (id, original_gtf_path))
            return (id, original_gtf_path, gene_structures.get_genes(original_gtf_path, gtf_cache_path,
                                                            logger))
    logger.info("Identifier %s was not found in the index file %s" %
                        (identifier, index_file_name))
    return (None, all_ids, None)

#_________________________________________________________________________________________
#
#    get_indexed_genes
#_________________________________________________________________________________________
def get_indexed_genes_matching_identifier(index_file_name, logger, regex_str):
    """
    Get gene structures contained in a GTF file whose identifier matches regex_str
    Because more than one species may match, returns a list of file names / genes

    :param index_file_name: path to index file created by ``index_gtf_files(...)``
    :type index_file_name: string
    :param regex_str: regular expression used to match identifiers parsed from the GTF file name by ``index_gtf_files(...)``
    :type regex_str: string
    :rtype list of tuples of (<matching identifier>, <original GTF path>, <dictionary of lists of genes>)
    """
    index_data = _read_index_file(index_file_name)
    regex = re.compile(regex_str)

    results = []
    #
    #   fields = <identifier><original_path><cache_path>
    #
    for id, original_gtf_path, gtf_cache_path in index_data:
        m = regex.search(id)
        if not m:
            continue
        gene_structures = gene.t_parse_gtf(id)
        logger.info("Get indexed genes for %s from %s" %
                            (id, original_gtf_path))
        results.append((id, original_gtf_path,
                        gene_structures.get_genes(original_gtf_path, gtf_cache_path,
                                                  logger)))
    if not len(results):
        logger.info("Regular expression %s did not match any entries in the index file %s" %
                            (regex_str, index_file_name))
    return results


#_________________________________________________________________________________________
#
#    get_indexed_genes_matching_gtf_file_name
#_________________________________________________________________________________________
def list_indexed_species(index_file_name):
    """
    Return a list of indexed species
    """
    return [ii[0 ] for ii in _read_index_file (index_file_name)]

#_________________________________________________________________________________________
#
#    get_indexed_genes_matching_gtf_file_name
#_________________________________________________________________________________________
def get_indexed_genes_matching_gtf_file_name(index_file_name, logger, regex_str):
    """
    Get gene structures contained in a GTF file whose file name matches regex_str
    Because more than one species may match, returns a list of file names / genes

    :param index_file_name: path to index file created by ``index_gtf_files(...)``
    :type index_file_name: string
    :param regex_str: regular expression used to match GTF file name by ``index_gtf_files(...)``
    :type regex_str: string
    :rtype list of tuples of (<matching identifier>, <original GTF path>, <dictionary of lists of genes>)
    """
    index_data = _read_index_file(index_file_name)
    regex = re.compile(regex_str)

    results = []
    #
    #   fields = <identifier><original_path><cache_path>
    #
    for id, original_gtf_path, gtf_cache_path in index_data:
        m = regex.search(original_gtf_path)
        if not m:
            continue
        gene_structures = gene.t_parse_gtf(id)
        logger.info("Get indexed genes for %s from %s" %
                            (id, original_gtf_path))
        results.append((id, original_gtf_path,
                        gene_structures.get_genes(original_gtf_path, gtf_cache_path,
                                                  logger)))
    if not len(results):
        logger.info("Regular expression %s did not match any entries in the index file %s" %
                            (regex_str, index_file_name))
    return results
#_________________________________________________________________________________________
#
#   _read_index_file
#_________________________________________________________________________________________
def _read_index_file (index_file_name):
    """
    Reads the index file created by ``index_gtf_files(...)``

    :param index_file_name: Path to the index file
    :type index_file_name: string
    :rtype List of [<id><original_gtf_path><gtf_cache_path>]

    """
    index_file = open(index_file_name)

    fields_list = []
    for line in index_file:
        if not len(line) or line[0] == '#':
            continue
        #
        #   fields = <identifier><original_path><cache_path>
        #
        fields = line.rstrip().split("\t")
        if len(fields) != 3:
            continue
        fields_list.append(fields)
    return fields_list

#_________________________________________________________________________________________
#
#    index_gtf_files
#_________________________________________________________________________________________

def index_gtf_files(index_file_name, search_path_root, regex_input, cache_file_pattern,
                    identifier_pattern, ignore_cache = False, logger = None):
    """
    Iterate through a directory, looking for all GTF files matching a regular expression.
    Cache the GTF data to a file and write the location of the file to an index file
    """
    if options.verbose > 1:
        logger.debug("Using regular expression %s" % regex_input)

    index_file = open(index_file_name,  "w")
    index_file_path = os.path.dirname(index_file_name)
    regex = re.compile(regex_input)
    for root, dirs, files in os.walk(search_path_root):
        for file_name in files:
            file_path = os.path.abspath(os.path.join(root, file_name))
            m = regex.search(file_path)
            if not m:
                continue
            if options.verbose > 1:
                logger.debug("Indexing %s" % file_name)
            cache_file_path  = regex.sub(cache_file_pattern.format(INDEX_FILE_PATH = index_file_path),
                                        file_path)
            identifier      = regex.sub(identifier_pattern , file_path)
            index_file.write("%s\t%s\t%s\n" % (identifier, file_path, cache_file_path))

            gene_structures = gene.t_parse_gtf(identifier)
            if logger:
                logger.debug("%s\t%s\t%s" %       (identifier, file_path, cache_file_path))
            gene_structures.index_genes (file_path, cache_file_path, logger, ignore_cache = ignore_cache)


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
    logger.info(" ".join(sys.argv))


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
                        options.cache_file_pattern,
                        options.identifier_pattern,
                        options.ignore_cache,
                        logger)

    #if options.debug:
    #    print "\n\nGetting genes from all species matching 'Ciona'"
    #    data = get_indexed_genes_matching_identifier(options.index_file,  logger,  "Ciona")
    #    for id, file_name, genes in data:
    #        print os.path.basename(file_name)
    #        print id, genes.keys()
    #        print "# of protein coding genes = ", len(genes['protein_coding'])
    #
    #    print "\n\nGetting genes from 'Ciona_savignyi:56'"
    #    id, file_name, genes = get_indexed_genes_for_identifier(options.index_file,  logger,  "Ciona_savignyi:56")
    #
    #    if genes:
    #        print os.path.basename(file_name)
    #        print id, genes.keys()
    #        print "# of protein coding genes = ", len(genes['protein_coding'])
    #
    #    print "\n\nGetting genes from all files matching 'Ciona'"
    #    data = get_indexed_genes_matching_gtf_file_name(options.index_file, logger, "Ciona")
    #    for id, file_name, genes in data:
    #        print os.path.basename(file_name)
    #        print id, genes.keys()
    #        print "# of protein coding genes = ", len(genes['protein_coding'])

    print "\n\nDone\n"


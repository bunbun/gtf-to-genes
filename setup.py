#!/usr/bin/env python
import ez_setup
ez_setup.use_setuptools()

import sys, os
if not sys.version_info[0:2] >= (2,4):
    sys.stderr.write("Requires Python later than 2.4\n")
    sys.exit(1)
    
# quickly import the latest version of ruffus
sys.path.insert(0, os.path.abspath("."))
import gtf_to_genes.gtf_to_genes_version
sys.path.pop(0)
    
    
if sys.version_info[0:2] >= (2,6):
    module_dependencies = []
else:
    module_dependencies = []


from setuptools import setup, find_packages
setup(
        name='gtf_to_genes',
        version=gtf_to_genes.gtf_to_genes_version.__version, #major.minor[.patch[.sub]]
        description='Fast GTF parser',
        long_description=\
"""     
***************************************
Overview
***************************************
    We want an extremely fast, lightweight way to access gene data stored in GTF format.
    
    The parsed data is held in an intuitive
        Gene
             -> transcript
             -> transcript
        with exons being stored as intervals
    
    Our aim is to 
       * cache data in binary format, which can be 
       * re-read in < 10s for even the largest genomes
    
    Currently initial parsing Ensembl Homo sapiens release 56 takes around 4.5 minutes.
    The binary data can be reloaded in < 10s.
    This contains *all* of the data structure in the original GTF file

    Note that we sacrifice memory usage for speed. This is seldom a problem for modern computers
    and genome sizes (There are around ~400,000 exons but there are stored as intervals / int pairs)

***************************************
A Simple example
***************************************
    ::
        gene_structures = t_parse_gtf("Mus musculus")

        #
        #   used cached data for speed
        #    
        ignore_cache = False
    
        # 
        #   get all protein coding genes only
        # 
        genes_by_type = gene_structures.get_genes(gtf_file, logger, ["protein_coding"], ignore_cache = ignore_cache)
    
        #
        #   print out gene counts
        #
        t_parse_gtf.log_gene_types (logger, genes_by_type)
    
        return genes_by_type        


""",
        author='Leo Goodstadt',
        author_email='gtf_to_genes@llew.org.uk',
        url='http://code.google.com/p/gtf-to-genes/',
    
        install_requires = module_dependencies, 
        setup_requires   = module_dependencies, 

        
        classifiers=[
                    'Intended Audience :: End Users/Desktop',
                    'Development Status :: 5 - Production/Stable',
                    'Intended Audience :: Developers',
                    'Intended Audience :: Science/Research',
                    'Intended Audience :: Information Technology',
                    'License :: OSI Approved :: MIT License',
                    'Programming Language :: Python',
                    'Topic :: Scientific/Engineering',
                    'Topic :: Scientific/Engineering :: Bio-Informatics',
                    'Environment :: Console',
                    ],
        license = "MIT",
        keywords = "GTF Ensembl gene transcript parser GFF bioinformatics science",


        packages=['gtf_to_genes'],
        package_dir={'gtf_to_genes': 'gtf_to_genes'},
        include_package_data = True,    # include everything in source control
        #package_data = {
        #    # If any package contains *.txt files, include them:
        #    '': ['*.TXT'],                                \
        #}


     )

#
#  http://pypi.python.org/pypi
#  http://docs.python.org/distutils/packageindex.html
#   
# 
# 
# python setup.py register
# python setup.py sdist --format=gztar upload

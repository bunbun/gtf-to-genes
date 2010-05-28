##################################################################################
:mod:`gtf_to-genes` -- Cache read of entire GTF files getting all gene structures
##################################################################################


    First index GTF file, then read all genes from a species a super high speed


##################################################################################
    Example:
##################################################################################

=========================================================
    To index entire directory tree
=========================================================            

    ::

        from gtf_to_genes import *
        import logging
        logger = logging.getLogger("test")

        index_file       = "/net/cpp-mirror/databases/ftp.ensembl.org/gtf.index"
        search_path_root = "/net/cpp-mirror/databases/ftp.ensembl.org"
        regex_input          = r"(.+\/)(([^.]+)\..+\.(.+)\.gtf(?:\.gz)?)$"

        # put cache file in same directory as GTF file
        cache_file_pattern   = r"\1\2.cache"

        # 
        # uncomment this line to put cache file in same directory index file
        #
        #cache_file_pattern   = r"{INDEX_FILE_PATH}/\2.cache"

        #
        # Unique identifier per GTF file
        # e.g. "Anolis_carolinensis:56"
        #
        identifier_pattern   = r"\3:\4"


        index_gtf_files(index_file,
                        search_path_root,
                        regex_input,
                        cache_file_pattern,
                        identifier_pattern,
                        logger)

=========================================================
    To read indexed genes from a particular species
=========================================================            

        ::

            index_file = "/net/cpp-mirror/databases/ftp.ensembl.org/gtf.index"
            from gtf_to_genes import *
            import logging
            logger = logging.getLogger("test")
            
            species, gtf_file_name, genes = get_indexed_genes_for_identifier(index_file,  logger,  "Ciona_savignyi:56")
            print species
            if genes:
                print genes.keys()
                print "# of protein coding genes = ", len(genes['protein_coding'])

        Produces::

            Ciona_savignyi:56
            ['pseudogene', 'snRNA', 'protein_coding', 'rRNA', 'miRNA', 'misc_RNA', 'snoRNA']
            # of protein coding genes =  11604
            
=========================================================
    To read indexed genes matching a regular expression
=========================================================
        ::
    
            index_file = "/net/cpp-mirror/databases/ftp.ensembl.org/gtf.index"
            from gtf_to_genes import *
            import logging
            logger = logging.getLogger("test")
            
            genes_with_identifiers  = get_indexed_genes_matching_identifier(index_file,  logger,  "Ciona")
            #genes_with_identifiers = get_indexed_genes_matching_gtf_file_name(index_file, logger, "Ciona")
            for id, file_path, genes in genes_with_identifiers:
                print "%30s (%d protein coding genes)" % (id, len(genes['protein_coding']))
                print file_path, "\n", genes.keys()

        Produces::

            Ciona_intestinalis:56         (14180 protein coding genes)
            /net/cpp-mirror/databases/ftp.ensembl.org/pub/release-56/gtf/ciona_intestinalis/Ciona_intestinalis.JGI2.56.gtf.gz
            ['rRNA', 'snRNA', 'protein_coding', 'miRNA', 'misc_RNA', 'snoRNA']
            
            Ciona_savignyi:56             (11604 protein coding genes)
            /net/cpp-mirror/databases/ftp.ensembl.org/pub/release-56/gtf/ciona_savignyi/Ciona_savignyi.CSAV2.0.56.gtf.gz
            ['pseudogene', 'snRNA', 'protein_coding', 'rRNA', 'miRNA', 'misc_RNA', 'snoRNA']


=========================================================
    Gene data
=========================================================
    Each gene (type = t_gene) looks like this::

        >>> genes = genes_with_identifiers[0][2]
        >>> genes['protein_coding'][20]
                beg          = 9203853
                coding_exons = [(9203910, 9204859), (9203919, 9204859)]
                contig       = '19'
                end          = 9204889
                exons        = [(9203853, 9204889)]
                gene_id      = 'ENSG00000170929'
                gene_type    = 'protein_coding'
                names        = ['OR1M1']
                strand       = True
                transcripts  = {   'ENST00000305465': 
                                       beg                 = 9203853
                                       cdna_id             = 'ENST00000305465'
                                       coding_beg          = 9203910
                                       coding_end          = 9204859
                                       coding_exon_indices = '\x00\x00'
                                       coding_frames       = '\x00\x00'
                                       end                 = 9203853
                                       exon_indices        = '\x00\x00'
                                       gene                = 'ENSG00000170929'
                                       names               = ['OR1M1-201']
                                       prot_ids            = ['ENSP00000303195']
                                       start_codons        = ()
                                       stop_codons         = (((9204858, 9204862), 0),),
                                   'ENST00000429566': 
                                       beg                 = 9203853
                                       cdna_id             = 'ENST00000429566'
                                       coding_beg          = 9203919
                                       coding_end          = 9204859
                                       coding_exon_indices = '\x01\x00'
                                       coding_frames       = '\x00\x00'
                                       end                 = 9203853
                                       exon_indices        = '\x00\x00'
                                       gene                = 'ENSG00000170929'
                                       names               = ['OR1M1-202']
                                       prot_ids            = ['ENSP00000401966']
                                       start_codons        = (((9203919, 9203923), 0),)
                                       stop_codons         = (((9204858, 9204862), 0),)}


#########################################################
    Data structures for genes/transcripts
#########################################################

==================================================================================================
Genes
==================================================================================================
.. class:: t_gene

	All co-ordinates use 0-based [) convention

	.. attribute::  gene_id

					Ensembl Gene identifier
	.. attribute::  contig

					Chromosome or contig
	.. attribute::  strand

					Strand
	.. attribute::  gene_type

					``protein_coding`` etc

                    Human genome includes: ``rRNA``, ``protein_coding``, ``snoRNA_pseudogene``, ``snRNA_pseudogene``, ``IG_V_gene``, ``misc_RNA``, ``misc_RNA_pseudogene``, ``IG_J_gene``, ``IG_C_gene``, ``Mt_tRNA``, ``Mt_rRNA``, ``scRNA_pseudogene``, ``pseudogene``, ``snRNA``, ``tRNA_pseudogene``, ``rRNA_pseudogene``, ``miRNA``, ``IG_D_gene``, ``processed_transcript``, ``Mt_tRNA_pseudogene``, ``snoRNA``, ``miRNA_pseudogene``

	.. attribute::  names

					Comma separated name of gene names
	.. attribute::  exons

					List of begin/ends 
	.. attribute::  coding_exons

					List of begin/ends
	.. attribute::  beg

					Begin of entire gene
	.. attribute::  end

					End of entire gene
	.. attribute::  transcripts

					Dictionary of ``t_transcript`` keyed by transcript identifier

    .. staticmethod:: load(data_file)

		Load data for gene from open file ``data_file``


    .. method:: dump(self, dump_file)

		Save data for gene to open file ``dump_file``


==================================================================================================
Transcript
==================================================================================================
.. class:: t_transcript

	.. attribute::  gene

					Link to parent gene
	.. attribute::  cdna_id

					Ensembl Transcript identifier (e.g. "ENST0000000004547")
	.. attribute::  prot_ids

					Ensembl Translation identifiers (e.g. "ENSP0000000004547")
	.. attribute::  names

					Comma separated name of transcript names
	.. attribute::  exon_indices

					Indices into gene.exons
	.. attribute::  coding_exon_indices

					Indices into gene.coding_exons
	.. attribute::  coding_frames

					List of coding frames
	.. attribute::  start_codons

					Indices of start codon(s)
	.. attribute::  stop_codons

					Indices of stop codon(s)
	.. attribute::  beg

					Genomic start coordinates
	.. attribute::  end

					Genomic end coordinates
	.. attribute::  coding_beg

					Genomic coding start coordinates
	.. attribute::  coding_end

					Genomic coding end coordinates

    .. staticmethod:: load(data_file)

		Load data for gene from open file ``data_file``


    .. method:: dump(self, dump_file)

		Save data for gene to open file ``dump_file``


#########################################################
Function Reference
#########################################################
=========================================================    
Indexing entire nested directory containing GTF files
=========================================================
.. function:: index_gtf_files(index_file, search_path_root, regex_input, cache_file_pattern, identifier_pattern)

    Iterate through a directory, looking for all GTF files matching a regular expression.
    Cache the GTF data to a file and write the location of the file to an index file
    
    :param index_file: Index file to hold list of parsed GTF files
    :param search_path_root: Root directory to start scanning for GTF files
    :param regex_input: Regular expression to match GTF or gzipped GTF files.
                       Brackets can be used to construct the corresponding cache file name.
                       E.g. ``r"(.+\/)(([^.]+)\..+\.(.+)\.gtf(?:\.gz)?)$"``
    :param cache_file_pattern: Pattern used to construct cache file name after regular expression substitution.
                       Brackets can be used to construct the corresponding cache file name
                       E.g. ``r"\1\2.cache" or r"{INDEX_FILE_PATH}/\2.cache"``
                       might give "/path/to/gtf/homo_sapiens.cache"
    :param identifier_pattern: Pattern used to constuct the (e.g. species) identifier for this GTF file
                       E.g. ``r"\1:\2" might give "homo_sapiens:47"``

=========================================================    
Genes for matching species
=========================================================

++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Get genes for specified species
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
.. function:: get_indexed_genes_for_identifier(index_file_name, logger, identifier)

    Get gene structures contained in a GTF file whose file name matches ``identifier``

    Returns None,None,None if no identifier matches
    
    :param index_file_name: path to index file created by ``index_gtf_files(...)``
    :type index_file_name: string
    :param identifier: identifier parsed from the GTF file name by ``index_gtf_files(...)``
    :type identifier: string
    :rtype: tuple of (``<matching identifier>``, ``<original GTF path>``, ``<dictionary of lists of genes>``)


++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Get genes for all species which match
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    Use a regular expression to specify either the species or the GTF file name

    More than one set of genes (species) may match

    Returns a list of ``identifier``, ``GTF file name``, ``genes``

.. function:: get_indexed_genes_matching_identifier(index_file_name, logger, regex_str)

    Get gene structures contained in a GTF file whose identifier matches regex_str

    :param index_file_name: path to index file created by ``index_gtf_files(...)``
    :type index_file_name: string
    :param regex_str: regular expression used to match identifiers parsed from the GTF file name by ``index_gtf_files(...)``
    :type regex_str: string 
    :rtype: list of tuples of (``<matching identifier>``, ``<original GTF path>``, ``<dictionary of lists of genes>``)

.. function:: get_indexed_genes_matching_gtf_file_name(index_file_name, logger, regex_str)

    Get gene structures contained in a GTF file whose file name matches regex_str

    :param index_file_name: path to index file created by ``index_gtf_files(...)``
    :type index_file_name: string
    :param regex_str: regular expression used to match GTF file name by ``index_gtf_files(...)``
    :type regex_str: string 
    :rtype: list of tuples of (``<matching identifier>``, ``<original GTF path>``, ``<dictionary of lists of genes>``)



License is MIT.



Contents:
==================

.. toctree::
   :maxdepth: 2

   missing_start_stop_codons.rst



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


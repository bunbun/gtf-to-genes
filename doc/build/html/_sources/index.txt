##################################################################################
:mod:`gtf_to_genes` Fast access to genes in GTF files
##################################################################################

==================================================================================================
Overview
==================================================================================================
    gtf_to_genes is a Python parser which caches all the genes / transcripts from a GTF file and
    caches the data into python classes for high speed access.

    The initial parsing and indexing of a GTF file is a slow process, taking over a minute. 
    But once that is done, the genes and transcripts can be retrieved in seconds.

==================================================================================================
Quick start
==================================================================================================
    The following code retrieves all the cached genes for Ensembl v.56 gene build of the Pacific sea squirt:
        ::

            species, gtf_file_name, genes = get_indexed_genes_for_identifier(index_file,  logger,  "Ciona_savignyi:56")

    | ``species`` and ``gtf_file_name`` are useful in checking that we are looking at the right data.
    | However ``genes`` is what we are interested in. These include not only protein coding but also rRNA and miRNA 
      genes etc.:

        ::

            >>> print genes.keys()
            ['pseudogene', 'snRNA', 'protein_coding', 'rRNA', 'miRNA', 'misc_RNA', 'snoRNA']
    
            >>> print "# of protein coding genes = ", len(genes['protein_coding'])
            # of protein coding genes =  11604

    A list of :ref:`t_gene <t_gene_class>` objects is stored for each gene type (e.g. ``protein_coding``). 
    This includes all exons / coding exons of each gene, as well as its constituent transcripts.

=========================================================
    Gene data
=========================================================
    Data for each gene is kept in a :ref:`t_gene <t_gene_class>`.


    Each :ref:`t_gene <t_gene_class>` looks like this::

        # 
        #   t_gene
        #
        {
            'gene_id'       : 'ENSG00000053438',
            'gene_type'     : 'protein_coding',
            'contig'        : '20',
            'beg'           : 36149615,
            'end'           : 36152092,
            'strand'        : True,
            'names'         : ['NNAT'],
            'exons'         : [(36149615, 36149805), (36149619, 36149805), (36150758, 36150840), (36151067, 36152092)],
            'coding_exons'  : [(36149732, 36149805), (36150758, 36150840), (36151067, 36151158)],
            'virtual_exons' : [(36149615, 36149805), (36150758, 36150840), (36151067, 36152092)],
            'transcripts'   : []
        }

    With the constituent :ref:`transcripts <t_transcript_class>` looking like this:

        Transcript 1: 
        ::

                  # 
                  #   t_transcript  1
                  #
                  {
                      'gene'                        : 'Back reference to <t_gene> object',
                      'cdna_id'                     : 'ENST00000062104',
                      'prot_id'                     : 'ENSP00000062104',
                      'beg'                         : 36149615,
                      'end'                         : 36152092,
                      'coding_beg'                  : 36149732,
                      'coding_end'                  : 36151158,
                      'names'                       : ['NNAT-001'],
                      'exon_indices'                : array('H', [0, 2, 3]),
                      'coding_exon_indices'         : array('H', [0, 1, 2]),
                      'coding_frames'               : array('H', [0, 0, 0]),
                      'virtual_exon_indices'        : array('H', [0, 1, 2]),
                      'virtual_coding_exon_indices' : array('H', [0, 1, 2]),
                      'start_codons'                : ((36149732, 36149736),),
                      'stop_codons'                 : ((36151157, 36151161),)
                  },

        Transcript 2: 
        ::

                  # 
                  #   t_transcript
                  #
                  {
                      'gene'                        : 'Back reference to <t_gene> object',
                      'cdna_id'                     : 'ENST00000346199',
                      'prot_id'                     : 'ENSP00000335497',
                      'beg'                         : 36149619,
                      'end'                         : 36152092,
                      'coding_beg'                  : 36149732,
                      'coding_end'                  : 36151158,
                      'names'                       : ['NNAT-002'],
                      'exon_indices'                : array('H', [1, 3]),
                      'coding_exon_indices'         : array('H', [0, 2]),
                      'coding_frames'               : array('H', [0, 0]),
                      'virtual_exon_indices'        : array('H', [0, 2]),
                      'virtual_coding_exon_indices' : array('H', [0, 2]),
                      'start_codons'                : ((36149732, 36149736),),
                      'stop_codons'                 : ((36151157, 36151161),)
                  }


=========================================================
    Exons
=========================================================
    Exon use for each transcript are stored by reference to the exon lists in its parent gene.
    This is both for efficiency and for ease of comparison. Exons that are actually shared 
    across transcripts may nonetheless have different coding / transcript start / stop sites.
    To facilitate comparisons, each gene also contains a "virtual" list with overlapping adjacent exons merged.

    Let us look at human gene ``ENSG00000053438`` in Ensembl v.58:

        ::

            species, gtf_file_name, genes = get_indexed_genes_for_identifier(index_file,  logger,  "Homo_sapiens:56")
            gi = index_genes_by_gene_id(genes)
            gene = gi['ENSG00000053438']

    It uses 4 different exons:

        ::

            >>> print gene.coding_exons
            [(36149615, 36149805), (36149619, 36149805), (36150758, 36150840), (36151067, 36152092)]

    Of these, only the third exon appears to be shared across the two transcripts:

        ::
    
            >>> print gene.transcripts[0].exon_indices
            >>> print gene.transcripts[1].exon_indices
            array('H', [0, 2, 3]),
            array('H', [1, 3]),

    However, ``exons[0]`` and ``exons[1]`` represent overlapping loci with different predicted transcription start sites.
    This is clear when we look at the ``virtual`` exon use (after merging overlaps):

        ::
    
            >>> print gene.transcripts[0].virtual_exon_indices
            >>> print gene.transcripts[1].virtual_exon_indices
            array('H', [0, 1, 2])
            array('H', [0,    2])

    Only the middle exon is actually alternatively spliced.


=========================================================
    Coordinates
=========================================================
    Following the standard python convention, and for ease of coordinate
    arithmetic, all coordinates are "zero-based, half-open".

    In other words, the first base in a contig or chromosome is base ``0``.
    If a sequence is 10 bases long, it may have the following loci:

        ::

            {
                'beg'           : 10,
                'end'           : 20,
            }

    The length is ``end - beg`` = ``20 - 10`` = ``10``

    N.B. the position of the last base in the sequence is base ``19`` (not ``20``).
    The ``end`` position is always one past the end of the sequence.

    This is the standard python convention, for example, in list slicing (``a_list[beg:end]``)

#########################################################
Classes for genes and transcripts
#########################################################

Gene and transcript data are stored in the following classes.
(All co-ordinates use 0-based [) convention.)

==================================================================================================
t_gene
==================================================================================================
.. _t_gene_class:

.. class:: t_gene

    The :ref:`t_gene <t_gene_class>` class supports the following methods and attributes:

.. _t_gene.gene_id:

    .. attribute::  gene_id

                    Ensembl Gene identifier
    .. attribute::  gene_type

                    ``protein_coding`` etc

                    Human genome includes: ``rRNA``, ``protein_coding``, ``snoRNA_pseudogene``, ``snRNA_pseudogene``, ``IG_V_gene``, ``misc_RNA``, ``misc_RNA_pseudogene``, ``IG_J_gene``, ``IG_C_gene``, ``Mt_tRNA``, ``Mt_rRNA``, ``scRNA_pseudogene``, ``pseudogene``, ``snRNA``, ``tRNA_pseudogene``, ``rRNA_pseudogene``, ``miRNA``, ``IG_D_gene``, ``processed_transcript``, ``Mt_tRNA_pseudogene``, ``snoRNA``, ``miRNA_pseudogene``

    .. attribute::  contig

                    Chromosome or contig
    .. attribute::  beg

                    Begin of entire gene
    .. attribute::  end

                    End of entire gene
    .. attribute::  strand

                    Strand
    .. attribute::  names

                    Comma separated name of gene names

.. _t_gene.exons:

    .. attribute::  exons

                    List of begins/ends of all exons used by transcripts. Note that some of these exons may overlap.
                    
.. _t_gene.coding_exons:

    .. attribute::  coding_exons

                    List of begins/ends of all exon coding sequence used by transcripts. . Note that some of these exons may overlap with each other. Begins and Ends are counted from the first coding base (including any start codon) to the last coding base (including any stop codon).

.. _t_gene.virtual_exons:

    .. attribute::  virtual_exons

                    List of begins/ends of loci from :ref:`t_gene.exons <t_gene.exons>` with overlaps merged together. 
                    These are useful for looking at transcript exon usage. 
                    For example, two transcripts may have overlapping first exons which differ by the length of the
                    5' UTR. These will have two separate entries in :ref:`t_gene.exons <t_gene.exons>` but will map to the same
                    "virtual exon" because they overlap.
    .. attribute::  transcripts

                    List of :ref:`t_transcript <t_transcript_class>` 

    .. staticmethod:: load(data_file)

        Load data for gene from open file ``data_file``


    .. method:: dump(self, dump_file)

        Save data for gene to open file ``dump_file``


==================================================================================================
t_transcript
==================================================================================================

.. _t_transcript_class:

.. class:: t_transcript

    The :ref:`t_transcript <t_transcript_class>` class supports the following methods and attributes:

    .. attribute::  gene

                    Link to parent :ref:`gene <t_gene_class>` 

.. _t_transcript.cdna_id:

    .. attribute::  cdna_id

                    Ensembl Transcript identifier (e.g. ``ENST0000000004547``)

.. _t_transcript.prot_id:

    .. attribute::  prot_id

                    Ensembl Translation identifiers (e.g. ``ENSP0000000004547``)
    .. attribute::  beg

                    Genomic start coordinates
    .. attribute::  end

                    Genomic end coordinates
    .. attribute::  coding_beg

                    Genomic coding start coordinates
    .. attribute::  coding_end

                    Genomic coding end coordinates

    .. attribute::  names

                    Comma separated name of transcript names
    .. attribute::  exon_indices

                    Indices into gene.exons
    .. attribute::  coding_exon_indices

                    Indices into :ref:`t_gene.coding_exons <t_gene.coding_exons>`

    .. attribute::  virtual_exon_indices

                    Indices into :ref:`t_gene.virtual_exons <t_gene.virtual_exons>`
    .. attribute::  virtual_coding_exon_indices

                    Indices into :ref:`t_gene.virtual_exons <t_gene.virtual_exons>`
    .. attribute::  coding_frames

                    List of coding frames
    .. attribute::  start_codons

                    Indices of start codon(s)
    .. attribute::  stop_codons

                    Indices of stop codon(s)
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




=====================================================================
Indexing genes / transcripts by identifier
=====================================================================
.. function:: index_transcripts(item, by_prot_id = True)

    Return a dictionary of transcripts indexed by either prot_id or cdna_id

    Recurses into a dictionary or list of genes or transcripts to find each 
    :ref:`t_transcript <t_transcript_class>`. Thus the ``genes`` returned by the GTF cache functions 
    (dictionaries of lists of genes) will be handled correctly.
    Objects other than :ref:`t_transcript <t_transcript_class>` or :ref:`t_gene <t_gene_class>` are otherwise ignored.

    :param item: Transcript/Gene (:ref:`t_transcript <t_transcript_class>` or :ref:`t_gene <t_gene_class>`) or ``list`` / ``dict`` or ``set`` of genes or transcripts.
    :param by_prot_id: Index by peptide :ref:`prot_id <t_transcript.prot_id>` instead of transcript identifier :ref:`cdna_id <t_transcript.cdna_id>`.
    :type by_prot_id: boolean
    :rtype: ``dict`` of :ref:`t_transcript <t_transcript_class>` indexed by :ref:`t_transcript.prot_id <t_transcript.prot_id>` or :ref:`t_transcript.cdna_id <t_transcript.cdna_id>`

.. function:: index_genes_by_gene_id(item)

    Return a dictionary of genes indexed by either :ref:`t_gene.gene_id <t_gene.gene_id>` 

    Recurses into a dictionary or list of genes to find each 
    :ref:`t_gene <t_gene_class>`. Thus the ``genes`` returned by the GTF cache functions 
    (dictionaries of lists of genes) will be
    handled correctly.
    Objects other than :ref:`t_transcript <t_transcript_class>` or :ref:`t_gene <t_gene_class>` are otherwise ignored.

    :param item: Gene (:ref:`t_gene <t_gene_class>`) or ``list`` / ``dict`` or ``set`` of genes.
    :rtype: ``dict`` of :ref:`t_gene <t_gene_class>` indexed by :ref:`t_gene.gene_id <t_gene.gene_id>`




##################################################################################
    Examples:
##################################################################################

=========================================================
    To index entire directory tree
=========================================================            
    To download the gtf files from an entire release of Ensembl:

    ::

        cd /net/cpp-mirror/databases/
        wget ftp://ftp.ensembl.org/pub/release-58/gtf -r -S 

    
    Index all downloaded GTF files

    ::

        from gtf_to_genes import *
        import logging
        logger = logging.getLogger("test")

        # Index  file
        index_file       = "/net/cpp-mirror/databases/ftp.ensembl.org/gtf.index"

        # look for GTF or gzipped GTF files
        regex_input          = r"(.+\/)(([^.]+)\..+\.(.+)\.gtf(?:\.gz)?)$"
        search_path_root = "/net/cpp-mirror/databases/ftp.ensembl.org"

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
            
=====================================================================
    To read indexed genes whose species match a regular expression
=====================================================================
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



##################################################################################
More examples
##################################################################################

.. toctree::
   :maxdepth: 2

   missing_start_stop_codons.rst
   find_overlaps.rst



##################################################################################
Indices and tables
##################################################################################

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


gtf_to_genes uses an MIT license.



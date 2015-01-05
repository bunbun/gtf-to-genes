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
    The following code retrieves all the cached genes for Ensembl v.77 human gene build:
        ::

            import sys
            import os
            from gtf_to_genes import *
            import logging
            logger = logging.getLogger("test")
            index_file = "/path/to/your/gtf.index"
            species_id, gtf_path, genes = get_indexed_genes_for_identifier(index_file, logger,
                                                                           "Homo_sapiens:77")

    | ``species`` and ``gtf_file_name`` are useful in checking that we are looking at the right data.
    | However ``genes`` is what we are interested in. These include not only protein coding but also rRNA and miRNA
      genes etc.:

        ::

            >>> print genes.keys()
            ['unitary_pseudogene', 'rRNA', 'lincRNA', 'IG_C_pseudogene',
            'translated_processed_pseudogene', 'Mt_tRNA', 'antisense',
            'IG_V_gene', 'misc_RNA', 'polymorphic_pseudogene', 'known_ncrna',
            'IG_J_gene', 'TR_J_pseudogene', 'IG_J_pseudogene', 'TEC',
            'protein_coding', 'Mt_rRNA', 'TR_V_pseudogene',
            '3prime_overlapping_ncrna', 'TR_J_gene', 'TR_D_gene',
            'IG_V_pseudogene', 'pseudogene', 'snRNA', 'unprocessed_pseudogene',
            'TR_V_gene', 'transcribed_unprocessed_pseudogene', 'sense_intronic',
            'miRNA', 'translated_unprocessed_pseudogene', 'non_coding',
            'IG_C_gene', 'sense_overlapping', 'IG_D_gene', 'TR_C_gene',
            'processed_transcript', 'transcribed_processed_pseudogene',
            'transcribed_unitary_pseudogene', 'snoRNA', 'processed_pseudogene']

            >>> print "# of protein coding genes = ", len(genes['protein_coding'])
            # of protein coding genes =  21983

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
            'gene_id'              : 'ENSMUSG00000073565',
            'gene_type'            : 'protein_coding',
            'genome_locus'         : 'chr18:51117897-51304641 :+',
            'gene_source'          : 'ensembl',
            'names'                : ['Prr16'],
            'exons'                : [(51117897, 51118089), (51302609, 51304641)],
            'coding_exons'         : [(51117930, 51118089), (51302609, 51303365)],
            'exon_ids'             : ['ENSMUSE00000707417', 'ENSMUSE00000707416'],
            'coding_exon_ids'      : ['ENSMUSE00000707417', 'ENSMUSE00000707416'],
            'virtual_exons'        : [(51117897, 51118089), (51302609, 51304641)],
            'virtual_coding_exons' : [(51117930, 51118089), (51302609, 51303365)],
            'transcripts'          : [],
            'contig'               : '18',
            'beg'                  : 51117897,
            'end'                  : 51304641,
            'strand'               : True


    With the constituent :ref:`transcripts <t_transcript_class>` looking like this:

        Transcript 1:
        ::

            #
            #   t_transcript  1
            #
            {
                'cdna_id'                     : 'ENSMUST00000116639',
                'prot_id'                     : 'ENSMUSP00000112338',
                'gene'                        : '<gene_id=ENSMUSG00000073565>',
                'names'                       : ['Prr16-201'],
                'cdna_source'                 : 'ensembl',
                'transcript_type'             : 'protein_coding',
                'beg'                         : 51117897,
                'end'                         : 51304641,
                'exons'                       : [(51117897, 51118089), (51302609, 51304641)],
                'exon_ids'                    : ['ENSMUSE00000707417', 'ENSMUSE00000707416'],
                'coding_beg'                  : 51117930,
                'coding_end'                  : 51303365,
                'coding_exons'                : [(51117930, 51118089), (51302609, 51303365)],
                'coding_exon_ids'             : ['ENSMUSE00000707417', 'ENSMUSE00000707416'],
                'start_codons'                : ((51117930, 51117933),),
                'stop_codons'                 : ((51303362, 51303365),),
                'coding_frames'               : array('H', [0, 0]),
                'exon_indices'                : array('H', [0, 1]),
                'coding_exon_indices'         : array('H', [0, 1]),
                'virtual_exon_indices'        : array('H', [0, 1]),
                'virtual_coding_exon_indices' : array('H', [0, 1])
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

            species, gtf_file_name, genes = get_indexed_genes_for_identifier(index_file,  logger,  "Homo_sapiens:77")
            gi = index_genes_by_gene_id(genes)
            gene = gi['ENSG00000053438']

    It uses 5 different exons:

        ::

            >>> print gene.coding_exons
            [(37521331, 37521403), (37521331, 37521403), (37522357, 37522438), (37522666, 37522759), (37522666, 37522759)]

    None of these are shared across the two transcripts

        ::

            >>> for transcript in gene.transcripts:
            ...     print transcript.exon_indices
            ...
            array('H', [0, 2, 4])
            array('H', [1, 3])

    However, ``exons[0]``, ``exons[1]`` and ``exons[3]``, ``exons[4]`` represent overlapping exons
    with different predicted transcription start sites.
    This is clear when we look at the ``virtual`` (after merging overlaps) exons:

        ::

            >>> for transcript in gene.transcripts:
            ...     print transcript.virtual_exon_indices
            ...
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

                    Human genome currently includes: ``unitary_pseudogene``, ``rRNA``, ``lincRNA``, ``IG_C_pseudogene``, ``translated_processed_pseudogene``, ``Mt_tRNA``, ``antisense``, ``IG_V_gene``, ``misc_RNA``, ``polymorphic_pseudogene``, ``known_ncrna``, ``IG_J_gene``, ``TR_J_pseudogene``, ``IG_J_pseudogene``, ``TEC``, ``protein_coding``, ``Mt_rRNA``, ``TR_V_pseudogene``, ``3prime_overlapping_ncrna``, ``TR_J_gene``, ``TR_D_gene``, ``IG_V_pseudogene``, ``pseudogene``, ``snRNA``, ``unprocessed_pseudogene``, ``TR_V_gene``, ``transcribed_unprocessed_pseudogene``, ``sense_intronic``, ``miRNA``, ``translated_unprocessed_pseudogene``, ``non_coding``, ``IG_C_gene``, ``sense_overlapping``, ``IG_D_gene``, ``TR_C_gene``, ``processed_transcript``, ``transcribed_processed_pseudogene``, ``transcribed_unitary_pseudogene``, ``snoRNA``, ``processed_pseudogene``

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
    .. attribute::  gene_source

                    Which method was used to predict the gene. E.g. ``ensembl_havana``

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

    .. method:: get_genome_locus(self)

        Return locus in the form of e.g. ``chr20:1234567-3456789 +``
        (Note the ``chr`` UCSC format)

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
    .. attribute::  transcript_type

                    index into list of transcript_types
                    For example::

                        >>> transcript.transcript_type_names[transcript.transcript_type]
                        protein_coding

    .. attribute::  beg

                    Genomic start coordinates
    .. attribute::  end

                    Genomic end coordinates
    .. attribute::  coding_beg

                    Genomic coding start coordinates
    .. attribute::  coding_end

                    Genomic coding end coordinates

    .. attribute::  cdna_source

                    Which method was used to predict the gene. E.g. ``ensembl_havana``
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

    .. method:: get_coding_exons(self)

        Returns list of coding exon bounds

    .. method:: get_coding_exon_ids(self)

        Returns list of coding exon identifiers

    .. method:: get_exons(self)

        Returns list of exon bounds

    .. method:: get_exon_ids(self)

        Returns list of coding exon identifiers

    .. method:: get_transcript_type_name(self)

        Returns the name of the transcript type, e.g. ``protein_coding``, ``retained_intron``, ``nonsense_mediated_decay``


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

    For example:
    ::

        cdna_by_cdna_id = index_transcripts(genes, False)
        prot_by_prot_id = index_transcripts(genes, True)


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
    To index entire directory tree in python
=========================================================
    To download the gtf files from an entire release of Ensembl:

    ::

        cd /your/gtf/path
        wget --recursive --timestamping --server-response ftp://ftp.ensembl.org/pub/current_gtf/homo_sapiens/Homo_sapiens.GRCh38.78.gtf.gz

        # file is now here:
        # /your/gtf/path/ftp.ensembl.org/pub/current_gtf/homo_sapiens/Homo_sapiens.GRCh38.78.gtf.gz


    Index all downloaded GTF files

    ::

        from gtf_to_genes import *
        import logging
        logger = logging.getLogger("test")

        # Index  file
        index_file       = "/your/gtf/path/gtf.index"

        # look for GTF or gzipped GTF files
        regex_input          = r"(.+\/)(([^.]+)\..+\.(.+)\.gtf(?:\.gz)?)$"
        search_path_root = "/your/gtf/path"

        # put cache file in same directory as GTF index file
        cache_file_pattern   = r"/your/gtf/path/\2.cache"

        #
        # uncomment this line to put cache file in same directory index file
        #
        #cache_file_pattern   = r"{INDEX_FILE_PATH}/\2.cache"

        #
        # Unique identifier per GTF file
        # e.g. "Anolis_carolinensis:77"
        #
        identifier_pattern   = r"\3:\4"


        index_gtf_files(index_file,
                        search_path_root,
                        regex_input,
                        cache_file_pattern,
                        identifier_pattern,
                        False,
                        logger)

=========================================================
    To index entire directory tree on the command line
=========================================================

    Download GTF file as above using ``wget``
    ::

        ../index_gtf_files.py -r /your/gtf/path --index /your/gtf/path/gtf.index -v -L index_gtf_files.log --ignore -v 5

=========================================================
    To read indexed genes from a particular species
=========================================================

    Let us read genes from the Pacific Sea Squirt:
        ::

            from gtf_to_genes import *
            import logging
            logger = logging.getLogger("test")

            index_file = "/your/gtf/path/gtf.index"
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



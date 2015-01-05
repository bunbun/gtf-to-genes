============================================================================
Example 1. Count number of genes with no start or stop codon
============================================================================

    Index *Homo sapiens* GTF file at ``"/your/gtf/path/ftp.ensembl.org/pub/release-77/gtf/homo_sapiens/"``
    (Change to your own location)
    ::

        from gtf_to_genes import *
        import logging
        logger = logging.getLogger("test")

        index_file       = "/your/gtf/path/ftp.ensembl.org/gtf.index"
        search_path_root = "/your/gtf/path/ftp.ensembl.org"
        search_path_root = "/your/gtf/path/ftp.ensembl.org/pub/release-77/gtf/homo_sapiens/"
        regex_input          = r"(.+\/)(([^.]+)\..+\.(.+)\.gtf(?:\.gz)?)$"

        # put cache file in same directory as GTF file
        cache_file_pattern   = r"\1\2.cache"

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
                        True,
                        logger)


    Get genes for Homo sapiens and count all transcripts with no stop or start codon
    ::

        from gtf_to_genes import *
        import logging
        logger = logging.getLogger("test")

        index_file       = "/your/gtf/path/ftp.ensembl.org/gtf.index"

        species, gtf_file_name, genes = get_indexed_genes_for_identifier(index_file,  logger,  "Homo_sapiens:77")
        print species
        if genes:
            print genes.keys()
            print "# of protein coding genes = ", len(genes['protein_coding'])



        #
        #   counts
        #
        cnt_gene = 0
        cnt_transcript = 0
        cnt_no_start_codon_transcript = 0
        cnt_no_stop__codon_transcript = 0
        cnt_neither__codon_transcript = 0
        cnt_no_start_codon_gene = 0
        cnt_no_stop__codon_gene = 0
        cnt_neither__codon_gene = 0
        cnt_no_start_codon_gene_any = 0
        cnt_no_stop__codon_gene_any = 0
        cnt_neither__codon_gene_any = 0

        #
        #   iterate through each gene
        #
        for g in genes['protein_coding']:
            cnt_gene += 1
            cnt_no_start_codon = 0
            cnt_no_stop__codon = 0
            cnt_neither__codon = 0
            #
            #   iterate through each transcript
            #
            for t in g.transcripts:
                cnt_transcript += 1
                if len(t.start_codons) == 0:
                    if len(t.stop_codons) == 0:
                        cnt_neither__codon += 1
                    else:
                        cnt_no_start_codon += 1
                elif len(t.stop_codons) == 0:
                    cnt_no_stop__codon += 1
            #
            #   Save counts for this gene
            #
            cnt_no_start_codon_transcript += cnt_no_start_codon
            cnt_no_stop__codon_transcript += cnt_no_stop__codon
            cnt_neither__codon_transcript += cnt_neither__codon
            cnt_no_start_codon_gene_any += 1 if cnt_no_start_codon else 0
            cnt_no_stop__codon_gene_any += 1 if cnt_no_stop__codon else 0
            cnt_neither__codon_gene_any += 1 if cnt_neither__codon else 0
            cnt_no_start_codon_gene_any += 1 if cnt_no_start_codon == len(g.transcripts) else 0
            cnt_no_stop__codon_gene_any += 1 if cnt_no_stop__codon == len(g.transcripts) else 0
            cnt_neither__codon_gene_any += 1 if cnt_neither__codon == len(g.transcripts) else 0


        #
        #   print summary
        #
        for i in range(1):
            print "%6d genes"                                                 % cnt_gene
            print "%6d transcripts"                                           % cnt_transcript
            print "%6d transcripts with no start codon"                       % cnt_no_start_codon_transcript
            print "%6d transcripts with no stop codon"                        % cnt_no_stop__codon_transcript
            print "%6d transcripts with no start or stop codon"               % cnt_neither__codon_transcript
            print "%6d genes with no start codon"                             % cnt_no_start_codon_gene
            print "%6d genes with no stop codon"                              % cnt_no_stop__codon_gene
            print "%6d genes with no start or stop codon"                     % cnt_neither__codon_gene
            print "%6d genes with any transcript with no start codon"         % cnt_no_start_codon_gene_any
            print "%6d genes with any transcript with no stop codon"          % cnt_no_stop__codon_gene_any
            print "%6d genes with any transcript with no start or stop codon" % cnt_neither__codon_gene_any


    Results:

        55.1% of transcripts have a missing stop or stop codon``(232 + 229 + 83681) / 152637.0``
        74.1% of genes have at least 1 transcript with a missing stop or stop codon``16286 / 21983.0``

    ::

         21983 genes
        152637 transcripts
           232 transcripts with no start codon
           229 transcripts with no stop codon
         83681 transcripts with no start or stop codon
             0 genes with no start codon
             0 genes with no stop codon
             0 genes with no start or stop codon
           207 genes with any transcript with no start codon
           253 genes with any transcript with no stop codon
         16286 genes with any transcript with no start or stop codon


    These are the equivalent numbers for Ensembl v.53:

        45.8% of transcripts have a missing stop or stop codon``(5993 + 9776 + 34401) / 109199.0``
        47.0% of genes have at least 1 transcript with a missing stop or stop codon``11097 / 23621.0``

    ::

         23621 genes
        109199 transcripts
          5993 transcripts with no start codon
          9776 transcripts with no stop codon
         34401 transcripts with no start or stop codon
             0 genes with no start codon
             0 genes with no stop codon
             0 genes with no start or stop codon
          4475 genes with any transcript with no start codon
          5070 genes with any transcript with no stop codon
         11097 genes with any transcript with no start or stop codon


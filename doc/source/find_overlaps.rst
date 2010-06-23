##################################################################################
Example 2. Find overlaps with any particular gene:
##################################################################################
=========================================================
1) Get list of genes for your desired species (as above)
=========================================================

        ::

            species_version = "Mus_musculus:58"

            index_file = "/net/cpp-mirror/databases/ftp.ensembl.org/gtf.index"
            from gtf_to_genes import *
            import logging
            logger = logging.getLogger("test")

            species, gtf_file_name, genes = get_indexed_genes_for_identifier(index_file,  logger,  species_version)
            print species
            if genes:
                print genes.keys()
                print "# of protein coding genes = ", len(genes['protein_coding'])

=================================================================
2) Use quicksect to build an overlap/intersection dictionary
=================================================================
    Enable fast lookups from loci to a particular gene.

    Coding / UTR go separately

        ::

            from quicksect import Feature, IntervalNode, IntervalTree
            def create_gene_quicksect (genes):
                """
                Create quicksect dictionary for looking up genes
                """
                gene_quicksect = {"coding": IntervalTree(), "utr": IntervalTree()}
                coding_or_utrs = sorted(gene_quicksect.keys())
                for genes in genes.itervalues():
                    for gene in genes:
                        # repeat for both coding and all (utr+coding) exons
                        for coding_or_utr, exons in zip(coding_or_utrs, [gene.coding_exons, gene.exons]):
                            for exon in exons:
                                gene_quicksect[coding_or_utr].insert(Feature(exon[0], exon[1],
                                                                         chr = gene.contig,
                                                                         #strand = gene.strand,
                                                                         info = {"gene_id"   : gene.gene_id,
                                                                                 "gene_type" : gene.gene_type,
                                                                                 "gene_name" : ", ".join(gene.names)
                                                                                 }))
                return gene_quicksect

            def create_transcript_quicksect (genes):
                """
                Create quicksect dictionary for looking up genes
                """
                transcript_quicksect = {"coding": IntervalTree(), "utr": IntervalTree()}
                coding_or_utrs = sorted(transcript_quicksect.keys())
                for genes in genes.itervalues():
                    for gene in genes:
                        for transcript in gene.transcripts.values():
                            for i in transcript.coding_exon_indices:
                                exon = gene.coding_exons[i]
                                transcript_quicksect["coding"].insert(Feature(exon[0], exon[1],
                                                                         chr = gene.contig,
                                                                         #strand = gene.strand,
                                                                         info = {"gene_id"       : gene.gene_id,
                                                                                 "transcript_id" : transcript.cdna_id,
                                                                                 "gene_type"     : gene.gene_type,
                                                                                 "gene_name"     : ", ".join(gene.names)
                                                                                 }))
                            for i in transcript.exon_indices:
                                exon = gene.exons[i]
                                transcript_quicksect["utr"].insert(Feature(exon[0], exon[1],
                                                                         chr = gene.contig,
                                                                         #strand = gene.strand,
                                                                         info = {"gene_id"       : gene.gene_id,
                                                                                 "transcript_id" : transcript.cdna_id,
                                                                                 "gene_type"     : gene.gene_type,
                                                                                 "gene_name"     : ", ".join(gene.names)
                                                                                 }))
                return transcript_quicksect


=================================================================
3) Function to find overlap
=================================================================

    ::

        from collections import defaultdict
        from quicksect import Feature

        def get_overlapping_gene (contig, beg, end, gene_quicksect_by_coding_or_utr):
            if contig[0:3] == "chr":
                contig = contig[3:]
            overlapping_genes = defaultdict(list)
            all_overlaps = set()
            for coding_or_utr, gene_quicksect in gene_quicksect_by_coding_or_utr.iteritems():
                overlaps = gene_quicksect.find(Feature(beg, end, chr = contig))
                if not overlaps:
                    continue
                for e in overlaps:
                    #print (e.info["gene_id"],
                    #       e.info["gene_type"],
                    #       e.info["gene_name"],
                    #       e.start,
                    #       e.stop)
                    all_overlaps.add((
                                  e.info["gene_type"],
                                  e.info["gene_id"],
                                  e.info["gene_name"]))
            for gene_type, gene_id, gene_name in sorted(all_overlaps):
                print "%20s %25s %s" % (gene_type, gene_id, gene_name)


        def get_overlapping_transcript (contig, beg, end, transcript_quicksect_by_coding_or_utr):
            if contig[0:3] == "chr":
                contig = contig[3:]
            overlapping_transcripts = defaultdict(list)
            all_overlaps = set()
            for coding_or_utr, transcript_quicksect in transcript_quicksect_by_coding_or_utr.iteritems():
                overlaps = transcript_quicksect.find(Feature(beg, end, chr = contig))
                if not overlaps:
                    continue
                for e in overlaps:
                    all_overlaps.add((
                                  e.info["gene_type"],
                                   e.info["gene_id"],
                                  e.info["transcript_id"],
                                  e.info["gene_name"]))
            for gene_type, gene_id, transcript_id, gene_name in sorted(all_overlaps):
                print "%20s %25s %25s %s" % (gene_type, gene_id, transcript_id, gene_name)




=================================================================
4) Does this work
=================================================================

    Make the overlap lookup quicksect

        ::

            gene_quicksect = create_gene_quicksect(genes)
            transcript_quicksect = create_transcript_quicksect(genes)

    Try some same regions

    http://genome.ucsc.edu/cgi-bin/hgTracks?db=mm9&position=chr2:50100000-52000000

        ::

            >>> get_overlapping_gene ("chr2", 50100000, 52000000, gene_quicksect)
            processed_transcript        ENSMUSG00000085014 AL805965.2
            processed_transcript        ENSMUSG00000085429 AL844592.3
            processed_transcript        ENSMUSG00000085862 BX649293.2
            processed_transcript        ENSMUSG00000086007 BX682229.1
            processed_transcript        ENSMUSG00000086349 AL935175.1
            processed_transcript        ENSMUSG00000086523 AL935175.2
                  protein_coding        ENSMUSG00000017144 Rnd3
                  protein_coding        ENSMUSG00000026766 Mmadhc
                  protein_coding        ENSMUSG00000026946 Nmi
                  protein_coding        ENSMUSG00000026950 Neb
                  protein_coding        ENSMUSG00000036202 Rif1
                  protein_coding        ENSMUSG00000036249 Rbm43
                  protein_coding        ENSMUSG00000053475 Tnfaip6
                  protein_coding        ENSMUSG00000056115 Tas2r134
                      pseudogene        ENSMUSG00000080782 AL929026.1
                      pseudogene        ENSMUSG00000081173 AL844893.1
                      pseudogene        ENSMUSG00000081457 AL844592.2
                      pseudogene        ENSMUSG00000081484 AL805965.1
                      pseudogene        ENSMUSG00000082483 BX649293.1
                      pseudogene        ENSMUSG00000082846 AL844592.1
                      pseudogene        ENSMUSG00000083270 AL844550.1
                      pseudogene        ENSMUSG00000083449 AL844592.4
                      pseudogene        ENSMUSG00000083472 BX005304.3
                      pseudogene        ENSMUSG00000084334 BX005304.1
                      pseudogene        ENSMUSG00000087701 RP23-332K7.1
                          snoRNA        ENSMUSG00000089614 SCARNA14

    http://genome.ucsc.edu/cgi-bin/hgTracks?db=mm9&position=chr10:69500000-69900000

        ::

            >>> get_overlapping_gene ("chr10", 69500000, 69900000, gene_quicksect)
                  protein_coding        ENSMUSG00000019933 2310015B20Rik
                  protein_coding        ENSMUSG00000037762 Slc16a9
                  protein_coding        ENSMUSG00000048701 Ccdc6
                      pseudogene        ENSMUSG00000052426 AC122923.1

    Transcripts in
    http://genome.ucsc.edu/cgi-bin/hgTracks?db=mm9&position=chr10:69500000-69900000

        ::

            >>> get_overlapping_transcript ("chr10", 69500000, 69900000, transcript_quicksect)
                  protein_coding        ENSMUSG00000019933        ENSMUST00000020090 2310015B20Rik
                  protein_coding        ENSMUSG00000037762        ENSMUST00000046807 Slc16a9
                  protein_coding        ENSMUSG00000037762        ENSMUST00000155933 Slc16a9
                  protein_coding        ENSMUSG00000048701        ENSMUST00000063086 Ccdc6
                  protein_coding        ENSMUSG00000048701        ENSMUST00000135607 Ccdc6
                  protein_coding        ENSMUSG00000048701        ENSMUST00000145990 Ccdc6
                  protein_coding        ENSMUSG00000048701        ENSMUST00000147545 Ccdc6
                  protein_coding        ENSMUSG00000048701        ENSMUST00000156001 Ccdc6
                      pseudogene        ENSMUSG00000052426        ENSMUST00000064271 AC122923.1


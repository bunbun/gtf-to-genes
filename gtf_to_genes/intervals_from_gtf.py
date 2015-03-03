#!/usr/bin/env python
"""

    intervals_from_gtf.py
    [--log_file PATH]
    [--verbose]

"""

################################################################################
#
#   intervals_from_gtf
#
#
#   Copyright (c) 03 December 2014 Leo Goodstadt
#
#################################################################################

import sys
import os

# add self to search path for testing
if __name__ == '__main__':
    exe_path = os.path.split(os.path.abspath(sys.argv[0]))[0]
    sys.path.insert(0, os.path.join(exe_path, ".."))


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#   imports


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

from gtf_to_genes import get_indexed_genes_for_identifier, __version__, list_indexed_species
from intervals import t_intervals, t_interval
from collections import defaultdict

#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#   options


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888


import argparse


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        # we do our own wrapping
        # argument descriptions start at column 31
        formatter_class=lambda prog: argparse.RawTextHelpFormatter(prog, max_help_position=31),
        description="""
Parses gene structure data from GTF files and prints lists of intervals.
""",
epilog="""
________________________________________________________________________

=============
Output files:
=============

  For convenience, intervals are output to files with default names in
  the --output OUTPUT directory unless the path is specified explicitly.

    1. Default name and directory
       --genes

          => './gene.loci'

    2. Default name, Specified directory
       --genes  --output '/my/path'

          => '/my/path/gene.loci'

    3. Specified path
       --genes 'genes.intervals' --output  '/my/path'

          => 'genes.intervals'

________________________________________________________________________

""")

    standard_options = argparse._ArgumentGroup(parser, 'Standard arguments')
    parser._action_groups.insert(0, standard_options)


    standard_options.add_argument("-i", "--index", metavar="PATH",
                        default = os.path.join(exe_path, "gtf.index"),
                        help="PATH to the index of all gtf files.")
    standard_options.add_argument("--print_species", action = "store_true",
                        help="Print list of indexed species and exit.")


    standard_options.add_argument("-s", "--species", metavar="NAME",
                      default = "Mus_musculus:60",
                      help="NAME and Ensembl version of species.\n"
                           "(default= 'Mus_musculus:77')")
    standard_options.add_argument("--genomic_bounds", metavar="PATH",
                      help="Tab-delimited file <CONTIG><START><END> in zero\n"
                           "based coordinates. Required for --intergenic,\n"
                           "--upstream,--downstream or --gene_territories")
    standard_options.add_argument("--gene_types", metavar="TYPES", action = 'append',
                      help="Comma separated list of valid gene types. All\n"
                            "other gene types will be ignored, e.g. for\n"
                            "the purposes of calculating gene territories.\n"
                            "E.g. 'protein_coding'. See --print_gene_types.")
    standard_options.add_argument('--print_gene_types', "-p", action="store_true",
                                help="Print a list of valid gene types and exit")


    output_files = parser.add_argument_group('Output files')
    output_files.add_argument("-G", "--genes", metavar="PATH", nargs='?', const = "DEFAULT",
                      help="Gene spans.\n(default= 'OUTPUT/genes.loci')")
    output_files.add_argument("-C", "--coding_exons", metavar="PATH", nargs='?', const = "DEFAULT",
                      help="(default= 'OUTPUT/coding_exons.loci')")
    output_files.add_argument("-E", "--exons", metavar="PATH", nargs='?', const = "DEFAULT",
                      help="Exons incl. UTR. (default= 'OUTPUT/exons.loci')")
    output_files.add_argument("-I", "--introns", metavar="PATH", nargs='?', const = "DEFAULT",
                      help="(default= 'OUTPUT/introns.loci')")
    output_files.add_argument("--coding_segments", metavar="PATH", nargs='?', const = "DEFAULT",
                      help="Coding regions with overlaps combined per gene.\n"
                           "(default= 'OUTPUT/coding_segments.loci')")
    output_files.add_argument("--exonic_segments", metavar="PATH", nargs='?', const = "DEFAULT",
                      help="Exonic regions with overlaps combined per gene.\n"
                           "(default= 'OUTPUT/exonic_segments.loci')")
    output_files.add_argument("--upstream", metavar="PATH", nargs='?', const = "DEFAULT",
                      help="5' FLANK. (default= 'OUTPUT/upstream.loci')")
    output_files.add_argument("--downstream", metavar="PATH", nargs='?', const = "DEFAULT",
                      help="3' FLANK. (default= 'OUTPUT/downstream.loci')")
    output_files.add_argument("--flanking_size", metavar="FLANK", type=int,
                                  default = 20000,
                                help="FLANK intervals (default= 20000).")
    output_files.add_argument("--intergenic", metavar="PATH", nargs='?', const = "DEFAULT",
                      help="(default= 'OUTPUT/intergenic.loci')")
    output_files.add_argument("--gene_territories", metavar="PATH", nargs='?', const = "DEFAULT",
                      help="Divide intergenic intervals between neighbours.\n"
                           "(default= 'OUTPUT/gene_territories.loci')")
    output_files.add_argument("-o", "--output", metavar="OUTPUT", default = ".",
                      help="Default output directory. (default= '.')")

    advanced_options = parser.add_argument_group('Advanced arguments')

    advanced_options.add_argument("--zero_based_coordinates", action="store_true",
                                   help="Output intervals using sensible UCSC count-\n"
                                        "from-zero rather than the Ensembl count-\n"
                                        "from-one coordinates.")
    advanced_options.add_argument("--exon_identifier",
                                  metavar= "FORMAT",
                                    default = "{gene_id}:{exon_id}:{chromosome}:{beg}-{end}:{strand}:{length}:{type}:{gene_type}",
                                  help = "Last column for --exons, --coding_exons\n"
                                         "Default= '{gene_id}:{exon_id}:{chromosome}:\n"
                                         "          {beg}-{end}:{strand}:{length}:{type}:{gene_type}'")
    advanced_options.add_argument("--segment_identifier",
                                  metavar= "FORMAT",
                                  default = "{gene_id}:{chromosome}:{beg}-{end}:{strand}:{length}:{type}:{gene_type}",
                                  help = "Last column for --coding_segments, --exonic_segments, --introns\n"
                                         "Default= '{gene_id}:{chromosome}:{beg}-{end}:\n"
                                         "          {strand}:{length}:{type}:{gene_type}'.")
    advanced_options.add_argument("--gene_identifier",
                                  metavar= "FORMAT",
                                  default = "{gene_id}:{chromosome}:{beg}-{end}:{strand}:{length}:{type}:{gene_type}",
                                  help = "Last column for --genes, --gene_territories,\n"
                                         "                --upstream and --downstream.\n"
                                         "Default= '{gene_id}:{chromosome}:{beg}-{end}:\n"
                                         "          {strand}:{length}:{type}:{gene_type}'.")


    common_group = parser.add_argument_group('Common arguments')
    common_group.add_argument('--verbose', "-v", const=1, metavar="VERBOSITY", default=0, nargs='?', type= int,
                                help="Print more verbose messages for each\nadditional verbose level.")
    common_group.add_argument('--version', action='version', version='%%(prog)s %s' % __version__)
    common_group.add_argument("-L", "--log_file", metavar="FILE", type=str,
                                  help="Name and path of log file.\n(Default= 'intervals_from_gtf.log')")


    options = parser.parse_args()

    #
    #   Defaults
    #
    if not options.log_file:
        options.log_file            = os.path.join("intervals_from_gtf.log")
    if options.coding_exons == "DEFAULT":
        options.coding_exons = os.path.join(options.output, "coding_exons.loci")
    if options.coding_segments == "DEFAULT":
        options.coding_segments = os.path.join(options.output, "coding_segments.loci")
    if options.exonic_segments == "DEFAULT":
        options.exonic_segments = os.path.join(options.output, "exonic_segments.loci")
    if options.genes == "DEFAULT":
        options.genes = os.path.join(options.output, "genes.loci")
    if options.exons == "DEFAULT":
        options.exons = os.path.join(options.output, "exons.loci")
    if options.introns == "DEFAULT":
        options.introns = os.path.join(options.output, "introns.loci")
    if options.upstream == "DEFAULT":
        options.upstream = os.path.join(options.output, "upstream.loci")
    if options.downstream == "DEFAULT":
        options.downstream = os.path.join(options.output, "downstream.loci")
    if options.intergenic == "DEFAULT":
        options.intergenic = os.path.join(options.output, "intergenic.loci")
    if options.gene_territories == "DEFAULT":
        options.gene_territories = os.path.join(options.output, "gene_territories.loci")
    if (options.gene_territories or options.intergenic or
        options.upstream or options.downstream) and not options.genomic_bounds:
        raise Exception("You need to specify the --genomic_bounds, i.e. the Start and End points "
                        "of each contig / chromosome for --intergenic or --gene_territories")


    if options.zero_based_coordinates:
        zero_based_coord = 0
    else:
        zero_based_coord = 1


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#   Logger


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

import logging
import logging.handlers

MESSAGE = 15
logging.addLevelName(MESSAGE, "MESSAGE")

def setup_std_logging (module_name, log_file, verbose):
    """
    set up logging using programme options
    """

    logger = logging.getLogger(module_name)


    # We are interesting in all messages
    logger.setLevel(logging.DEBUG)
    has_handler = False

    #
    #   log to file if that is specified
    #
    if log_file:
        handler = logging.FileHandler(log_file, delay=False)
        class stripped_down_formatter(logging.Formatter):
            def format(self, record):
                prefix = ""
                if not hasattr(self, "first_used"):
                    self.first_used = True
                    prefix = "\n" + self.formatTime(record, "%Y-%m-%d")
                    prefix += " %(name)s\n" % record.__dict__
                if record.levelname in ("INFO", "MESSAGE", "DEBUG"):
                    self._fmt = " %(asctime)s - %(message)s"
                else:
                    self._fmt = " %(asctime)s - %(levelname)-7s - %(message)s"
                return prefix + logging.Formatter.format(self, record)
        handler.setFormatter(stripped_down_formatter("%(asctime)s - %(name)s - %(levelname)6s - %(message)s", "%H:%M:%S"))
        handler.setLevel(MESSAGE)
        logger.addHandler(handler)
        has_handler = True

    #
    #   log to stderr if verbose
    #
    if verbose:
        stderrhandler = logging.StreamHandler(sys.stderr)
        stderrhandler.setFormatter(logging.Formatter("    %(message)s"))
        stderrhandler.setLevel(logging.DEBUG)
        if log_file:
            class debug_filter(logging.Filter):
                """
                Ignore INFO messages
                """
                def filter(self, record):
                    return logging.INFO != record.levelno
            stderrhandler.addFilter(debug_filter())
        logger.addHandler(stderrhandler)
        has_handler = True

    #
    #   no logging
    #
    if not has_handler:
        class NullHandler(logging.Handler):
            """
            for when there is no logging
            """
            def emit(self, record):
                pass
        logger.addHandler(NullHandler())


    return logger

if __name__ == '__main__':

    #
    #   set up log: name = script name sans extension
    #
    module_name = os.path.splitext(os.path.basename(sys.argv[0]))[0];
    logger = setup_std_logging(module_name, options.log_file, options.verbose)


    #
    #   log programme parameters
    #
    logger.info(" ".join(sys.argv))



#__________________________________________________________________________________________________
#
#   read_contig_extents
#__________________________________________________________________________________________________
def read_contig_extents():
    """
    Read extents of each chromosome from file
    Expecting 0-based coordinates
    """
    if not options.genomic_bounds:
        return {}

    contig_extents = dict()
    with open(options.genomic_bounds) as genomic_file:
        for ii, line in enumerate(genomic_file):
            if not len(line) or line[0] == '#':
                continue
            line = line.rstrip()
            fields = line.split("\t")
            try:
                if len(fields) != 3:
                    raise Exception("bad")
                chrom, beg, end = fields
                beg = int(beg)
                end = int(end)
                contig_extents[chrom] = t_interval(beg, end)
                if chrom[:3] == "chr":
                    contig_extents[chrom[3:]] = t_interval(beg, end)
            except:
                raise Exception("Non-conforming line # %d in --genomic_bounds %s. Expecting "
                                "tab delimited <CHROM><BEG><END>. Found %d fields.\n(%s)"
                                % (ii + 1,
                                   options.genomic_bounds,
                                   len(fields), line))


    return contig_extents

#__________________________________________________________________________________________________
#
#   write_exons
#__________________________________________________________________________________________________
def write_exons(file_handle, gene, iterator, exon_identifier_format, interval_type):
    if file_handle:
        for exon, exon_id in iterator:
            file_handle.write(exon_identifier_format.format(
                              gene    = gene,
                              gene_id = gene.gene_id,
                              chromosome = gene.contig,
                              exon_id = exon_id,
                              beg     = exon[0] + zero_based_coord,
                              end     = exon[1],
                              strand  = "-+"[gene.strand],
                              type = interval_type,
                              gene_type= gene.gene_type,
                              length  = exon[1] - exon[0]))

#__________________________________________________________________________________________________
#
#   write_virtual_exons
#__________________________________________________________________________________________________
def write_virtual_exons(file_handle, gene, iterator, exon_identifier_format, interval_type):
    # no exon_id for combined exonic regions
    if file_handle:
        for exon in iterator:
            file_handle.write(exon_identifier_format.format(
                              gene    = gene,
                              gene_id = gene.gene_id,
                              chromosome = gene.contig,
                              beg     = exon[0] + zero_based_coord,
                              end     = exon[1],
                              strand  = "-+"[gene.strand],
                              type    = interval_type,
                              gene_type= gene.gene_type,
                              length  = exon[1] - exon[0]))

#__________________________________________________________________________________________________
#
#   write_genes
#__________________________________________________________________________________________________
def write_gene(file_handle, gene, beg, end, gene_identifier_format, interval_type):
    file_handle.write(gene_identifier_format.format(chromosome=gene.contig,
                                                   beg = beg + zero_based_coord,
                                                   end = end,
                                                   length = end - beg,
                                                   gene_id = gene.gene_id,
                                                   gene_type = gene.gene_type,
                                                   strand  = "-+"[gene.strand],
                                                   type=interval_type))


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#   Main logic


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
if __name__ == '__main__':

    if options.print_species:
        for species in list_indexed_species(options.index):
            print species
        sys.exit()


    #
    #   get genes previously parsed by index_gtf_files.py
    #
    if options.verbose >= 2:
        logger.log(MESSAGE, "Reading GTF data for %s from %s " % (options.species, options.index))
    species, gtf_file_name, genes = get_indexed_genes_for_identifier(options.index,
                                                                     logger,
                                                                     options.species)
    if not genes:
        raise Exception("No genes for %s from %s" % (options.species, options.index))
    if not "protein_coding" in genes:
        raise Exception("No protein coding genes for %s" % (options.species,))

    if options.print_gene_types:
        logger.info("Valid gene types:")
        sys.stdout.write("Valid gene types:\n")
        for gene_type in sorted(genes.keys()):
            sys.stdout.write("    %s\n" % gene_type)
            logger.info("    %s\n" % gene_type)
        sys.exit()

    open_file_names = set()

    #
    #   output files
    #
    def open_file(file_name, open_files):
        # do not output
        if not file_name:
            return None
        # Shared output
        if file_name not in open_files:
            open_files[file_name] = open(file_name, "w")
        return open_files[file_name]


    open_files = dict()
    coding_exons_file           = open_file(options.coding_exons         , open_files)
    coding_segments_file        = open_file(options.coding_segments      , open_files)
    exonic_segments_file        = open_file(options.exonic_segments      , open_files)
    genes_file                  = open_file(options.genes                , open_files)
    exons_file                  = open_file(options.exons                , open_files)
    introns_file                = open_file(options.introns              , open_files)
    upstream_file               = open_file(options.upstream             , open_files)
    downstream_file             = open_file(options.downstream           , open_files)
    gene_territories_file       = open_file(options.gene_territories     , open_files)
    intergenic_file             = open_file(options.intergenic           , open_files)

    if options.verbose >= 2:
        logger.log(MESSAGE, "Using gene structures from %s..." % (gtf_file_name,))
        for option, name in ( (options.coding_exons    , "coding exons"),
                              (options.coding_segments , "combined coding segments"),
                              (options.exonic_segments , "combined exonic segments"),
                              (options.genes           , "genes"),
                              (options.exons           , "exons incl. UTR"),
                              (options.introns         , "merged intronic segments"),
                              (options.upstream        , "5' %d bp regions" % options.flanking_size),
                              (options.downstream      , "3' %d bp regions" % options.flanking_size),
                              (options.gene_territories, "gene territories"),
                              (options.intergenic      , "intergenic regions"),
                              ):
            if option:
                logger.log(MESSAGE, "Writing intervals for %-24s to %s" % (name, option  ))






    contig_extents = read_contig_extents()



    contigs_without_extents = set()

    logger.debug("Writing GTF data")
    genic_intervals_per_contig = defaultdict(t_intervals)

    # how exons / segments are described
    exon_identifier_format    = ("{chromosome}\t{beg}\t{end}\t%s\n" % options.exon_identifier   ).decode("string_escape")
    segment_identifier_format = ("{chromosome}\t{beg}\t{end}\t%s\n" % options.segment_identifier).decode("string_escape")
    gene_identifier_format    = ("{chromosome}\t{beg}\t{end}\t%s\n" % options.gene_identifier   ).decode("string_escape")

    if options.gene_types is None or not len(options.gene_types):
        valid_gene_types = genes.keys()
    else:
        import re
        valid_gene_types = set(re.split("[,\s]+", ",".join(options.gene_types)))
        if not len(options.gene_types):
            raise Exception("No valid gene types specified in --gene_types")

    # only protein_coding genes
    for gene_type in valid_gene_types:
        if not gene_type in genes:
            logger.log(MESSAGE, "Warning: Gene type '%s' was not found in the GTF file "
                                "and will be ignored..." % (gene_type,))
        for gene in genes[gene_type]:
            # gene spans
            if genes_file:
                # CHROM BEG END GENE_ID
                write_gene(genes_file, gene, gene.beg, gene.end,
                           gene_identifier_format, "gene")
                genic_intervals_per_contig[gene.contig].append((gene.beg, gene.end, gene))


            # (overlapping) coding exons / exon + UTR
            # CHROM BEG END GENE_ID:EXON_ID:CHROM:BEG:END:LENGTH
            write_exons(exons_file, gene, gene.get_exons(), exon_identifier_format, "exon")
            write_exons(coding_exons_file, gene, gene.get_coding_exons(), exon_identifier_format, "coding_exon")

            # (combined) coding exonic / intronic regions per gene
            # CHROM BEG END GENE_ID:CHROM:BEG:END:LENGTH
            write_virtual_exons(coding_segments_file, gene, gene.get_virtual_coding_exons(),
                                segment_identifier_format, "CDS")
            write_virtual_exons(exonic_segments_file, gene, gene.get_virtual_exons(),
                                segment_identifier_format, "exonic")
            write_virtual_exons(introns_file, gene, gene.get_virtual_introns(),
                                segment_identifier_format, "intronic")


            if gene.contig in contig_extents:
                contig_beg, contig_end = contig_extents[gene.contig]
            else:
                contigs_without_extents.add(gene.contig)
                contig_beg, contig_end = 0, gene.end + options.flanking_size


            # make sure that bounds are within contig bounds, but not so far that the upstream downstream extents become negative
            if gene.strand:
                if upstream_file    :
                    write_gene(upstream_file, gene,
                               min((gene.beg, max((contig_beg, gene.beg - options.flanking_size)))),
                               gene.beg,
                               gene_identifier_format, "5flank")
                if downstream_file  :
                    write_gene(downstream_file, gene,
                               gene.end,
                               max((gene.end, min((contig_end, gene.end + options.flanking_size)))),
                               gene_identifier_format, "3flank")
            else:
                if downstream_file  :
                    write_gene(downstream_file, gene,
                               min([gene.beg, max([contig_beg, gene.beg - options.flanking_size])]),
                               gene.beg,
                               gene_identifier_format, "3flank")
                if upstream_file    :
                    write_gene(upstream_file, gene,
                               gene.end,
                               max((gene.end, min((contig_end, gene.end + options.flanking_size)))),
                               gene_identifier_format, "5flank")




    if gene_territories_file or intergenic_file:
        for contig in genic_intervals_per_contig:

            # known  contig bounds
            if contig in contig_extents:
                contig_beg, contig_end = contig_extents[contig]

            # end at last gene
            else:
                contig_beg, contig_end = 0, genic_intervals_per_contig[contig].data[-1][1]
                contigs_without_extents.add(contig)

            # complement of intergenic
            if intergenic_file:
                intergenic_regions = genic_intervals_per_contig[contig].complemented(contig_beg, contig_end).combine_overlapping().remove_empty()
                for beg, end in intergenic_regions.data:
                    file_handle.write(segment_identifier.format(chromosome=contig,
                                                                   beg = beg + zero_based_coord,
                                                                   end = end,
                                                                   length = end - beg,
                                                                   gene_id = "",
                                                                   gene_type = "",
                                                                   strand  = "+",
                                                                   type="intergenic"))


                    intergenic_file.write("{contig}\t{beg}\t{end}\t\n".format(beg = beg, end = end, contig=contig))

            # gene territories
            if gene_territories_file:
                genic_intervals_per_contig[contig].extend_into_gaps(contig_beg, contig_end)
                for beg, end, gene in genic_intervals_per_contig[contig].data:
                    write_gene(gene_territories_file, gene,
                               beg, end,
                               gene_identifier_format, "GENE_TERRITORY")

    # invert to get intronic
    if len(contigs_without_extents):
        logger.log(MESSAGE, "")
        logger.log(MESSAGE, "=" * 80)
        logger.log(MESSAGE, "WARNING: The following contigs do not have known bounds: ")
        contigs_without_extents  = sorted(contigs_without_extents, key = lambda x: (len(x), x))
        logger.log(MESSAGE, ", ".join(contig for contig in contigs_without_extents))
        logger.log(MESSAGE, "=" * 80)
        logger.log(MESSAGE, "")



    #
    #   Expand to gene territorities
    #

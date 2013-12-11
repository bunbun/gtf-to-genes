#!/usr/bin/env python
################################################################################
#
#   gene
#
#       ALL coordinates are stored 0-based
#
#
#   Copyright (c) 3/1/2010 Leo Goodstadt
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
FILE_VERSION_MAJ = 5
FILE_VERSION_MIN = 1234


import sys, os
from collections import defaultdict
import minimal_gtf_iterator
import time
from dump_object import dump_object
from itertools import izip
from array import array
import struct
import gzip
import marshal
do_dump = marshal.dump
do_load = marshal.load

from random_access_file_by_sections import fill_directory_of_sections, write_directory_of_sections, read_directory_of_sections


# how exonic data is stored
EXON_NUMBER = 0
GENOMIC_POS = 1
EXON_DATA   = 2
EXON_ID     = 3

INTERVAL_BEG = 0
INTERVAL_END = 1


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#   Functions


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
#_____________________________________________________________________________________
#
#   median
#_____________________________________________________________________________________
def median (l):
    if not len(l):
        return None
    median_pos = int(len(l) / 2)
    return sorted(l)[median_pos]


#_____________________________________________________________________________________
#
#   non_str_sequence
#_____________________________________________________________________________________
def non_str_sequence (arg):
    """
    Whether arg is a sequence.
    We treat strings / dicts however as a singleton not as a sequence

    """
    if (isinstance(arg, (basestring, dict))):
        return False
    try:
        test = iter(arg)
        return True
    except TypeError:
        return False

#_____________________________________________________________________________________
#
#   overlapping_combined
#_____________________________________________________________________________________
def overlapping_combined( orig_data, reverse = False):
    """
    Return list of intervals with overlapping neighbours merged together
    Assumes sorted intervals unless reverse is set

    """
    if not orig_data or not len(orig_data): return []
    if len(orig_data) == 1:
        return orig_data

    new_data = []

    if reverse:
        data = orig_data[:]
        data.reverse()
    else:
        data = orig_data

    if not data[0][0] <= data[1][0]:
        print data, reverse
    assert(data[0][0] <= data[1][0])

    # start with the first interval
    prev_beg, prev_end = data[0]

    # check if any subsequent intervals overlap
    for beg, end in data[1:]:
        if beg - prev_end > 0:
            new_data.append((prev_beg, prev_end))
            prev_beg = beg
        prev_end = max(end, prev_end)

    new_data.append((prev_beg, prev_end))

    if reverse:
        new_data.reverse()
    return new_data

#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#   t_gene

#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
class t_gene(object):
    """
    """

    #_____________________________________________________________________________________
    #
    #   __init__
    #_____________________________________________________________________________________
    def __init__(self, gene_id, contig, strand, gene_type, names, exons, coding_exons, exon_ids, coding_exon_ids):
        self.gene_id      = gene_id
        self.contig       = contig
        self.strand       = strand
        self.gene_type    = gene_type
        self.names        = names
        self.exons        = exons
        self.coding_exons = coding_exons
        self.exon_ids        = exon_ids
        self.coding_exon_ids = coding_exon_ids
        self.beg          = min(e[INTERVAL_BEG] for e in self.exons)
        self.end          = max(e[INTERVAL_END] for e in self.exons)
        self.transcripts  = []

        #self.exons        = sorted(self.exons       , reverse = not self.strand)
        #self.coding_exons = sorted(self.coding_exons, reverse = not self.strand)
        self.virtual_exons= overlapping_combined(exons, reverse = not self.strand)

    __slots__ = [
                         "gene_id",
                         "gene_type",
                         "contig",
                         "beg",
                         "end",
                         "strand",
                         "names",
                         "exons",
                         "coding_exons",
                         "exon_ids",
                         "coding_exon_ids",
                         "virtual_exons",
                         "transcripts",
                         # unused and undefined
                         #"transcripts_types",
                        ]
    def __eq__(self, other):
        if type(other) == t_gene:
            return cmp(self, other) == 0
        else:
            return False

    def __cmp__(self, other):
        if type(other) != t_gene:
            return -1
        return (
                cmp(self.contig           , other.contig          ) or
                cmp(self.beg              , other.beg             ) or
                cmp(self.end              , other.end             ) or
                cmp(self.strand           , other.strand          ) or
                cmp(self.gene_type        , other.gene_type       ) or
                cmp(self.gene_id          , other.gene_id         ) or
                cmp(self.names            , other.names           ) or
                cmp(self.exons            , other.exons           ) or
                cmp(self.coding_exons     , other.coding_exons    ) or
                cmp(self.exon_ids         , other.exon_ids        ) or
                cmp(self.coding_exon_ids  , other.coding_exon_ids ) or
                cmp(self.transcripts      , other.transcripts     ))



    #_____________________________________________________________________________________
    #
    #   add_transcript
    #_____________________________________________________________________________________
    def add_transcript (self, transcript):
        self.transcripts.append(transcript)


    #_____________________________________________________________________________________
    #
    #   __repr__
    #_____________________________________________________________________________________
    def __repr__ (self):
        return dump_object(self, sorted_attribute_names = self.__slots__)

    #_____________________________________________________________________________________
    #
    #   dump
    #_____________________________________________________________________________________
    def dump(self, data_file):
        """
        dump
        """
        do_dump(self.gene_id        , data_file)
        do_dump(self.contig         , data_file)
        do_dump(self.strand         , data_file)
        do_dump(self.gene_type      , data_file)
        do_dump(self.names          , data_file)
        do_dump(self.exons          , data_file)
        do_dump(self.coding_exons   , data_file)
        do_dump(self.exon_ids       , data_file)
        do_dump(self.coding_exon_ids, data_file)
        do_dump(len(self.transcripts), data_file)
        for transcript in self.transcripts:
            transcript.dump(data_file)

    #_____________________________________________________________________________________
    #
    #   load
    #_____________________________________________________________________________________
    @staticmethod
    def load(data_file):
        """
        load
        """
        gene_id            = do_load(data_file)
        contig             = do_load(data_file)
        strand             = do_load(data_file)
        gene_type          = do_load(data_file)
        names              = do_load(data_file)
        exons              = do_load(data_file)
        coding_exons       = do_load(data_file)
        exon_ids           = do_load(data_file)
        coding_exon_ids    = do_load(data_file)

        gene = t_gene(gene_id, contig, strand, gene_type, names, exons, coding_exons, exon_ids, coding_exon_ids)
        gene.virtual_exons= overlapping_combined(exons, reverse = not gene.strand)
        cnt_transcripts = do_load(data_file)
        for i in range(cnt_transcripts):
            transcript = t_transcript.load(gene, data_file)
            gene.add_transcript(transcript)
        return gene

    #_____________________________________________________________________________________
    #
    #   get_transcript_type_names
    #_____________________________________________________________________________________
    def get_transcript_type_names (self):
        names = set(t.get_transcript_type_name() for t in self.transcripts)
        return tuple(sorted(names))




#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#   t_transcript


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
class t_transcript(object):
    transcript_type_names = []
    unique_transcript_type_names = dict()

    def get_indices_into_virtual_exons (self, exons, exon_indices):
        """
        Translate from coding /non-coding exons into indices of exons after overlapping
        loci have been combined. This allows exon usage across transcripts to be compared
        easily
        """
        virtual_exon_indices = []
        for ei in exon_indices:
            exon = exons[ei]
            for i, ve in enumerate(self.gene.virtual_exons):
                if exon[INTERVAL_BEG] >= ve[INTERVAL_BEG] and exon[INTERVAL_END] <= ve[INTERVAL_END]:
                    virtual_exon_indices.append(i)
                    break
            else:
                raise Exception("Exon [%d-%d] not found in parent exons (%s) of transcript %s" %
                                    (exon[INTERVAL_BEG], exon[INTERVAL_END], self.cdna_id, self.gene.gene_id))
        return array('H', virtual_exon_indices)



    #_____________________________________________________________________________________
    #
    #   __init__
    #_____________________________________________________________________________________
    def __init__(self, gene, cdna_id, names, prot_id, exon_indices, coding_exon_indices,
                    coding_frames, start_codons, stop_codons, transcript_type,
                    beg = None, end = None, coding_beg = None, coding_end = None):
        """
        Initialise
        """
        self.gene                   = gene
        self.cdna_id                = cdna_id
        self.names                  = names
        self.prot_id                = prot_id
        self.exon_indices           = exon_indices
        self.coding_exon_indices    = coding_exon_indices
        self.coding_frames          = coding_frames
        self.start_codons           = start_codons
        self.stop_codons            = stop_codons
        self.transcript_type        = transcript_type
        if not len(self.gene.exons):
            raise Exception("No exons for transcript %s (%s) in gene %s" % (cdna_id, prot_id, gene.gene_id))
        self.beg                    = min(self.gene.exons[e][INTERVAL_BEG] for e in self.exon_indices) if beg == None else beg
        self.end                    = max(self.gene.exons[e][INTERVAL_END] for e in self.exon_indices) if end == None else end
        if coding_beg != None:
            self.coding_beg = coding_beg
        else:
            self.coding_beg      = min(self.gene.coding_exons[e][INTERVAL_BEG] for e in self.coding_exon_indices) if len(self.coding_exon_indices) else None
        if coding_end != None:
            self.coding_end = coding_end
        else:
            self.coding_end      = max(self.gene.coding_exons[e][INTERVAL_END] for e in self.coding_exon_indices) if len(self.coding_exon_indices) else None

        #print >>sys.stderr, "gene.exons=", self.gene.exons
        #print >>sys.stderr, "exon_indices=", self.exon_indices


        self.virtual_exon_indices = self.get_indices_into_virtual_exons (self.gene.exons, self.exon_indices)
        self.virtual_coding_exon_indices = self.get_indices_into_virtual_exons (self.gene.coding_exons, self.coding_exon_indices)
    __slots__ = [
                         "gene",
                         "cdna_id",
                         "prot_id",
                         "beg",
                         "end",
                         "coding_beg",
                         "coding_end",
                         "names",
                         "exon_indices",
                         "coding_exon_indices",
                         "coding_frames",
                         "virtual_exon_indices"       ,
                         "virtual_coding_exon_indices",
                         "start_codons",
                         "stop_codons",
                         "transcript_type"
                        ]

    def __eq__(self, other):
       if type(other) == t_transcript:
           return cmp(self, other) == 0
       else:
           return False

    def __cmp__(self, other):
       if type(other) != t_transcript:
           return -1
       return (
               cmp(self.cdna_id             , other.cdna_id            ) or
               cmp(self.names               , other.names              ) or
               cmp(self.prot_id             , other.prot_id            ) or
               cmp(self.exon_indices        , other.exon_indices       ) or
               cmp(self.coding_exon_indices , other.coding_exon_indices) or
               cmp(self.coding_frames       , other.coding_frames      ) or
               cmp(self.start_codons        , other.start_codons       ) or
               cmp(self.stop_codons         , other.stop_codons        ) or
               cmp(self.beg                 , other.beg                ) or
               cmp(self.end                 , other.end                ) or
               cmp(self.coding_beg          , other.coding_beg         ) or
               cmp(self.coding_end          , other.coding_end         ) or
               cmp(self.transcript_type     , other.transcript_type         )
               )
    #_____________________________________________________________________________________
    #
    #   dump
    #_____________________________________________________________________________________
    def dump(self, data_file):
        """
        dump
        """
        do_dump(self.cdna_id                 , data_file)
        do_dump(self.names                   , data_file)
        do_dump(self.prot_id                 , data_file)
        do_dump(self.exon_indices            , data_file)
        do_dump(self.coding_exon_indices     , data_file)
        do_dump(self.coding_frames           , data_file)
        do_dump(self.start_codons            , data_file)
        do_dump(self.stop_codons             , data_file)
        do_dump(self.beg                     , data_file)
        do_dump(self.end                     , data_file)
        do_dump(self.coding_beg              , data_file)
        do_dump(self.coding_end              , data_file)
        do_dump(self.transcript_type         , data_file)

    #_____________________________________________________________________________________
    #
    #   load
    #_____________________________________________________________________________________
    @staticmethod
    def load(gene, data_file):
        """
        load
        """
        cdna_id                 = do_load(data_file)
        names                   = do_load(data_file)
        prot_id                 = do_load(data_file)
        exon_indices            = array('H', do_load(data_file))
        coding_exon_indices     = array('H', do_load(data_file))
        coding_frames           = array('H', do_load(data_file))
        start_codons            = do_load(data_file)
        stop_codons             = do_load(data_file)
        beg                     = do_load(data_file)
        end                     = do_load(data_file)
        coding_beg              = do_load(data_file)
        coding_end              = do_load(data_file)
        transcript_type         = do_load(data_file)

        transcript = t_transcript(gene               ,
                                  cdna_id            ,
                                  names              ,
                                  prot_id            ,
                                  exon_indices       ,
                                  coding_exon_indices,
                                  coding_frames      ,
                                  start_codons       ,
                                  stop_codons        ,
                                  transcript_type    ,
                                  beg                ,
                                  end                ,
                                  coding_beg         ,
                                  coding_end         )
        transcript.virtual_exon_indices = transcript.get_indices_into_virtual_exons (transcript.gene.exons, transcript.exon_indices)
        transcript.virtual_coding_exon_indices = transcript.get_indices_into_virtual_exons (transcript.gene.coding_exons, transcript.coding_exon_indices)
        return transcript



    #_____________________________________________________________________________________
    #
    #   __repr__
    #_____________________________________________________________________________________
    def __repr__ (self):
        saved_gene = self.gene
        self.gene = "<Reference to gene with gene_id=%s>" % self.gene.gene_id
        self_str = dump_object(self, sorted_attribute_names = self.__slots__)
        self.gene = saved_gene
        return self_str


    #_____________________________________________________________________________________
    #
    #   get_transcript_type_name
    #_____________________________________________________________________________________
    def get_transcript_type_name (self):
        return t_transcript.transcript_type_names[self.transcript_type]



#_____________________________________________________________________________________
#
#   sort_intervals_by_exon_number_then_pos
#_____________________________________________________________________________________
def sort_by_exon_number_then_pos (intervals):
    if not intervals:
        return []
    # None
    if iter(intervals).next()[0] != None:
        return [data[EXON_DATA:] for data in sorted(intervals)]
    else:
        return [data[EXON_DATA:] for data in sorted(intervals, reverse = not strand, key = lambda frame: frame[GENOMIC_POS:])]



class t_parse_gtf(object):
    #_____________________________________________________________________________________
    #
    #   construct_gene_list_from_parsed_data
    #_____________________________________________________________________________________
    def construct_gene_list_from_parsed_data (self, logger):
        """
        Use parse data to construct list of genes
        """
        genes_by_type = defaultdict(list)
        for unique_gene_id in self.gene_id_to_contig_strands:
            gene_id, contig, strand = unique_gene_id

            #
            #   construct gene
            #
            gene_names          = list(self.gene_id_to_names[unique_gene_id])
            gene_type           = self.gene_id_to_gene_type[unique_gene_id]

            gene_exons          = set()
            gene_coding_exons   = set()

            # get all exons and coding exons from transcripts
            for cdna_id in self.gene_id_to_cdna_ids[unique_gene_id]:
                gene_exons        |= self.cdna_id_to_exons[cdna_id]
                gene_coding_exons |= self.cdna_id_to_coding_exons[cdna_id]



            # sort by position
            gene_exons                  = sorted(gene_exons,        reverse = not strand, key = lambda exon: exon[1:])
            gene_coding_exons           = sorted(gene_coding_exons, reverse = not strand, key = lambda exon: exon[1:])
            gene_exon_intervals         = [interval for exon_number, pos, interval, exon_id in gene_exons]
            gene_exon_ids               = [exon_id  for exon_number, pos, interval, exon_id in gene_exons]
            gene_coding_exon_intervals  = [interval for exon_number, pos, interval, exon_id, frame in gene_coding_exons]
            gene_coding_exon_ids        = [exon_id  for exon_number, pos, interval, exon_id, frame in gene_coding_exons]

            new_gene = t_gene(gene_id, contig, strand, gene_type,
                              gene_names, gene_exon_intervals, gene_coding_exon_intervals, gene_exon_ids, gene_coding_exon_ids)

            #print >>sys.stderr, "new_gene=", new_gene

            #
            # turn gene_exons into dictionary of indices
            #
            gene_exon_intervals        = dict(izip(gene_exon_intervals,        xrange(len(gene_exons))))
            gene_coding_exon_intervals = dict(izip(gene_coding_exon_intervals, xrange(len(gene_coding_exons))))

            #print >>sys.stderr, gene_exon_intervals, gene_coding_exon_intervals

            genes_by_type[gene_type].append(new_gene)

            #
            #   construct transcripts
            #
            for cdna_id in self.gene_id_to_cdna_ids[unique_gene_id]:
                names                =  self.cdna_id_to_names[cdna_id]
                if len(names) > 1:
                    logger.warning("Transcript %s has %d (>1) names (%s)" %
                                     cdna_id, len(names), ", ".join(names))
                prot_id             = self.cdna_id_to_prot_id[cdna_id] if cdna_id in self.cdna_id_to_prot_id else None
                exons               = sort_by_exon_number_then_pos (self.cdna_id_to_exons       [cdna_id])
                coding_exons        = sort_by_exon_number_then_pos (self.cdna_id_to_coding_exons[cdna_id])
                exon_indices        = array('H', [gene_exon_intervals[e[0]] for e in exons])
                coding_exon_indices = array('H', [gene_coding_exon_intervals[e[0]] for e in coding_exons])

                #
                #   sort by start position
                #   Each coding frame / start / stop codon is tagged with the beg position of that exon
                #   Assume each transcript does not have overlapping exons!!!
                #
                # strip all but frame
                coding_frames =  array('H', [e[2] for e in coding_exons])


                start_codons  =  tuple(sort_by_exon_number_then_pos (self.cdna_id_to_start_codons[cdna_id]))
                stop_codons   =  tuple(sort_by_exon_number_then_pos (self.cdna_id_to_stop_codons [cdna_id]))

                # check exon numbers correct
                if len(coding_frames) != len(coding_exons):
                    raise Exception("Not all CDS Exons have frames for %s/%s" %
                                     (gene_id, cdna_id))


                #
                # hash transcript_type to same string
                #
                transcript_type_name = self.cdna_id_to_transcript_type[cdna_id]
                if transcript_type_name not in t_transcript.unique_transcript_type_names:
                    t_transcript.unique_transcript_type_names[transcript_type_name] = len(t_transcript.transcript_type_names)
                    t_transcript.transcript_type_names.append(transcript_type_name)
                transcript_type = t_transcript.unique_transcript_type_names[transcript_type_name]
                transcript = t_transcript(  new_gene,
                                            cdna_id,
                                            list(names),
                                            prot_id,
                                            exon_indices,
                                            coding_exon_indices,
                                            coding_frames,
                                            start_codons,
                                            stop_codons,
                                            transcript_type)

                new_gene.add_transcript(transcript)

        return genes_by_type

    #_____________________________________________________________________________________
    #
    #   log_gene_types
    #_____________________________________________________________________________________
    @staticmethod
    def log_gene_types (logger, genes_by_type):
        """
        Log the different sorts of genes
        """
        #
        # Count the different types of genes
        #

        gene_counts_by_types                   = dict()
        gene_counts_by_transcript_types        = defaultdict(lambda: defaultdict(int))
        transcript_counts_by_types             = defaultdict(int)
        exon_counts_by_types                   = defaultdict(int)
        transcript_exon_counts_by_types        = defaultdict(int)
        transcript_coding_exon_counts_by_types = defaultdict(int)
        all_merged_exon_counts_by_types        = defaultdict(list)
        all_merged_coding_exon_counts_by_types = defaultdict(list)

        try:
            for gene_type, genes in genes_by_type.iteritems():
                gene_counts_by_types[gene_type] =len(genes)
                for g in genes:
                    gene_counts_by_transcript_types[gene_type][g.get_transcript_type_names()]+=1

                    transcript_counts_by_types[gene_type] += len(g.transcripts)
                    exon_counts_by_types[gene_type] += len(g.exons)
                    all_merged_exon_counts_by_types[gene_type].append(len(overlapping_combined(g.exons, not g.strand)))
                    all_merged_coding_exon_counts_by_types[gene_type].append(len(overlapping_combined(g.coding_exons, not g.strand)))
                    for t in g.transcripts:
                        transcript_exon_counts_by_types[gene_type] += len(t.exon_indices)
                        transcript_coding_exon_counts_by_types[gene_type] += len(t.coding_exon_indices)



        except:
            print g
            raise




        #
        # print out protein-coding and pseudogene counts first
        #
        priority_types = ["protein_coding", "pseudogene", "retrotransposed"]
        for gt in priority_types:
            if gt in gene_counts_by_types:
                if "gene" in gt:
                    logger.info("   %ss." % (gt.capitalize()))
                else:
                    logger.info("   %s genes." % (gt.capitalize()))

                if gt in gene_counts_by_transcript_types:
                    for transcript_types in sorted(gene_counts_by_transcript_types[gt].keys()):
                        logger.info("       %8d genes with %s transcripts" % (gene_counts_by_transcript_types[gt][transcript_types], transcript_types))

                logger.info("     %8d genes." % (gene_counts_by_types[gt]))
                logger.info("     %8d transcripts." % (transcript_counts_by_types[gt]))
                logger.info("     %8d exons." % (exon_counts_by_types[gt]))

                logger.info("     %8d exon regions."           % (sum(all_merged_exon_counts_by_types[gt])))
                logger.info("     %8d exon regions at median." % (median(all_merged_exon_counts_by_types[gt])))
                logger.info("     %8d exon regions at max."    % (max(all_merged_exon_counts_by_types[gt])))
                logger.info("     %8d exon regions at min."    % (min(all_merged_exon_counts_by_types[gt])))

                if len(all_merged_coding_exon_counts_by_types[gt]):
                    logger.info("     %8d exon coding regions."           % (sum(all_merged_coding_exon_counts_by_types[gt])))
                    logger.info("     %8d exon coding regions at median." % (median(all_merged_coding_exon_counts_by_types[gt])))
                    logger.info("     %8d exon coding regions at max."    % (max(all_merged_coding_exon_counts_by_types[gt])))
                    logger.info("     %8d exon coding regions at min."    % (min(all_merged_coding_exon_counts_by_types[gt])))

                logger.info("     %8d exons in all transcripts." % (transcript_exon_counts_by_types[gt]))
                if transcript_coding_exon_counts_by_types[gt]:
                    logger.info("     %8d coding exons in all transcripts." % (transcript_coding_exon_counts_by_types[gt]))

        #
        # print out other counts
        #
        for gene_type, cnts in sorted(gene_counts_by_types.items()):
            if gene_type in priority_types:
                continue
            logger.info("     %5d %-20s genes." % (cnts, gene_type))
            if gene_type in gene_counts_by_transcript_types:
                for transcript_types in sorted(gene_counts_by_transcript_types[gene_type].keys()):
                    logger.info("       %8d genes with %s transcripts" % (gene_counts_by_transcript_types[gene_type][transcript_types], transcript_types))

    #_____________________________________________________________________________________
    #
    #   save_genes_to_cache
    #_____________________________________________________________________________________
    def save_genes_to_cache (self, genes_by_type, CACHED_RESULTS, logger):
        start_time = time.time()

        # make sure directory exists
        CACHED_RESULTS_DIR = os.path.split(CACHED_RESULTS)[0]
        if not os.path.exists(CACHED_RESULTS_DIR):
            os.makedirs(CACHED_RESULTS_DIR)

        data_file = open(CACHED_RESULTS, 'wb', 5)
        data_file.write(struct.pack("q", FILE_VERSION_MAJ))
        data_file.write(struct.pack("q", FILE_VERSION_MIN))

        do_dump(t_transcript.transcript_type_names, data_file)
        do_dump(t_transcript.unique_transcript_type_names, data_file)

        # use placeholder of the correct size to hold file positions
        file_pos_placeholder = long(data_file.tell())


        #
        #   save "directory" recording where features of each type are in the files so we
        #       can jump directly to that type
        section_name_directory_entry_pos = write_directory_of_sections(data_file, genes_by_type.keys())

        #
        #   save actual data for each gene_type in a separate section
        #
        file_pos_by_section = {}
        for gene_type in genes_by_type.keys():
            file_pos_by_section[gene_type] = data_file.tell()

            marshal.dump(len(genes_by_type[gene_type]), data_file)
            for gene in genes_by_type[gene_type]:
                gene.dump(data_file)

        #
        #   go back and write section starts in the directory
        #
        fill_directory_of_sections(data_file, section_name_directory_entry_pos, file_pos_by_section)


        end_time = time.time()
        logger.info("  Cached in %ds" % (end_time - start_time))

    #_____________________________________________________________________________________
    #
    #   cache_is_valid
    #_____________________________________________________________________________________
    def cache_is_valid (self, CACHED_RESULTS, logger):
        """
        Load genes from cache file
        """
        start_time = time.time()
        try:
            if not os.path.exists(CACHED_RESULTS):
                return False

            data_file = open(CACHED_RESULTS, 'rb')

            #
            #   check version
            #
            file_version1 = struct.unpack("q", data_file.read(8))[0]
            file_version2 = struct.unpack("q", data_file.read(8))[0]
            if file_version1 != FILE_VERSION_MAJ or file_version2 != FILE_VERSION_MIN:
                logger.info("  Cache file wrong version = %d.%d (should be %d.%d)" %
                            (file_version1, file_version2,
                             FILE_VERSION_MAJ, FILE_VERSION_MIN))
                return False

            #
            #   read directory of gene type sections
            #
            file_pos_by_section = read_directory_of_sections (data_file)

        except:
            return False

        return True

    #_____________________________________________________________________________________
    #
    #   load_genes_from_cache
    #_____________________________________________________________________________________
    def load_genes_from_cache (self, CACHED_RESULTS, logger, valid_gene_types = None):
        """
        Load genes from cache file
        """
        start_time = time.time()
        try:
            if not os.path.exists(CACHED_RESULTS):
                return None

            logger.debug("  Loading genes from cache")
            logger.info("  Cache file = %s" % CACHED_RESULTS)
            data_file = open(CACHED_RESULTS, 'rb')
            file_version1 = struct.unpack("q", data_file.read(8))[0]
            file_version2 = struct.unpack("q", data_file.read(8))[0]
            if file_version1 != FILE_VERSION_MAJ or file_version2 != FILE_VERSION_MIN:
                logger.info("  Cache file wrong version = %d.%d (should be %d.%d)" %
                            (file_version1, file_version2,
                             FILE_VERSION_MAJ, FILE_VERSION_MIN))
                return None

            t_transcript.transcript_type_names          = do_load(data_file)
            t_transcript.unique_transcript_type_names   = do_load(data_file)

            #
            #   read directory of gene type sections
            #
            file_pos_by_section = read_directory_of_sections (data_file)

            if valid_gene_types == None:
                valid_gene_types = file_pos_by_section.keys()


            #
            #   Read data for all valid gene type
            #
            genes_by_type = defaultdict(list)
            for gene_type in valid_gene_types:
                if gene_type not in file_pos_by_section:
                    continue

                data_file.seek(file_pos_by_section[gene_type], os.SEEK_SET)

                num_genes = marshal.load(data_file)
                for i in range(num_genes):
                    genes_by_type[gene_type].append(t_gene.load(data_file))

            #
            #   Save a log of the different gene types
            #
            self.log_gene_types (logger, genes_by_type)

            end_time = time.time()
            cnt_genes = sum(len(genes) for genes in genes_by_type.values())
            logger.info("  Loaded %d genes in %ds" % (cnt_genes,  end_time - start_time))
            return genes_by_type
        except:
            raise
            return None

    #_____________________________________________________________________________________
    #
    #   get_genes
    #_____________________________________________________________________________________
    def index_genes (self, gtf_file_path, cached_gtf_file_path = None, logger = None, ignore_cache = False):
        """
        get_genes
        """
        logger.info(self.species_name)

        if cached_gtf_file_path == None:
            ignore_cache = True

        #
        #   load from cache
        #
        if not ignore_cache:
            if not os.path.exists(cached_gtf_file_path):
                logger.debug("%s has not been cached will be reparsed" %
                             (  os.path.basename(gtf_file_path, )))
            elif self.cache_is_valid(cached_gtf_file_path, logger):
                return
            else:
                logger.debug("%s is not a valid cache file. %s will be reparsed" %
                             (os.path.basename(cached_gtf_file_path),
                                os.path.basename(gtf_file_path)))


        #
        #   construct from file
        #
        genes_by_type = self.construct_gene_list_from_file (gtf_file_path, logger)

        #
        #   save to cache
        #
        if cached_gtf_file_path != None:
            self.save_genes_to_cache (genes_by_type, cached_gtf_file_path, logger)

        return

    #_____________________________________________________________________________________
    #
    #   get_genes
    #_____________________________________________________________________________________
    def get_genes (self, gtf_file_path, cached_gtf_file_path = None, logger = None, valid_gene_types = None, ignore_cache = False):
        """
        get_genes
        """
        logger.info(self.species_name)

        if cached_gtf_file_path == None:
            ignore_cache = True

        #
        #   load from cache
        #
        if not ignore_cache:
            genes_by_type = self.load_genes_from_cache (cached_gtf_file_path, logger, valid_gene_types)
            if genes_by_type:
                return genes_by_type

        #
        #   construct from file
        #
        genes_by_type = self.construct_gene_list_from_file (gtf_file_path, logger)

        #
        #   save to cache
        #
        if cached_gtf_file_path != None:
            self.save_genes_to_cache (genes_by_type, cached_gtf_file_path, logger)


        if valid_gene_types != None:
            for gene_type in genes_by_type.keys():
                if gene_type not in valid_gene_types:
                    del genes_by_type[gene_type]

        return genes_by_type



    #_____________________________________________________________________________________
    #
    #   construct_gene_list_from_file
    #_____________________________________________________________________________________
    def construct_gene_list_from_file (self, gtf_file_path, logger):
        """
        construct genes from GTF file
        """

        #
        #   parse into useful structures and check all is well
        #
        logger.debug("  Parsing genes data from GTF file.")
        logger.info("  GTF file = %s" % gtf_file_path)
        start_time = time.time()
        self.parse(gtf_file_path, logger)
        self.check_number_of_names (logger)
        self.check_split_codons (logger)
        end_time = time.time()
        logger.info("  Parsed in %ds" % (end_time - start_time))
        start_time = time.time()

        #
        #   build genes
        #
        logger.debug("  Creating data structures for genes.")
        genes_by_type = self.construct_gene_list_from_parsed_data(logger)
        end_time = time.time()

        cnt_all_genes = sum(len(g) for g in genes_by_type.itervalues())
        logger.info("  Constructed %d genes in %ds" % (cnt_all_genes,  end_time - start_time))
        self.log_gene_types (logger, genes_by_type)

        return genes_by_type







    #_____________________________________________________________________________________
    #
    #   __init__
    #_____________________________________________________________________________________
    def __init__(self, species_name):
        # list of all features
        self.feature_set = set()
        self.gene_types  = set()
        self.species_name= species_name


        #
        #   gene_id -> beg/end/cdna_id/exon_index
        #   gene_id -> gene_name
        #           -> contig / strand
        #
        self.gene_id_to_contig_strands = set()
        self.gene_id_to_names          = defaultdict(set)
        self.gene_id_to_cdna_ids       = defaultdict(set)
        self.gene_id_to_gene_type      = dict()

        #
        #   cdna_id
        #
        self.cdna_id_to_names        = defaultdict(set)
        self.cdna_id_to_prot_id      = dict()
        self.cdna_id_to_start_codons = defaultdict(set)
        self.cdna_id_to_stop_codons  = defaultdict(set)
        self.cdna_id_to_exons        = defaultdict(set)
        self.cdna_id_to_coding_exons = defaultdict(set)
        self.cdna_id_to_transcript_type= dict()


    #    Ensemble Mus_musculus.GRCm38.74.gtf
    #
    #    gtf_entry:
    #        mContig
    #        mSource
    #        mFeature,
    #        mStart
    #        mEnd
    #        mScore
    #        mStrand
    #        mFrame
    #    mSource
    #            protein_coding, TR_V_gene, miRNA, Mt_rRNA, Mt_tRNA, rRNA, snoRNA, snRNA
    #            IG_C_gene, IG_D_gene, IG_J_gene, IG_LV_gene, IG_V_gene
    #            lincRNA, 3prime_overlapping_ncrna
    #            misc_RNA,
    #            nonsense_mediated_decay, non_stop_decay, retained_intron, sense_intronic, sense_overlapping, antisense
    #            pseudogene, processed_transcript, processed_pseudogene, unprocessed_pseudogene, polymorphic_pseudogene, unitary_pseudogene
    #            transcribed_processed_pseudogene, transcribed_unprocessed_pseudogene
    #            translated_processed_pseudogene, translated_unprocessed_pseudogene
    #            IG_V_pseudogene, TR_V_pseudogene
    #    mContig     =   1-19, MT, X, Y, GL456210.1, GL456211.1 ...
    #    mFeature    =   CDS, exon, start_codon, stop_codon
    #    mScore      =   .
    #    mStrand     =   + -
    #    mFrame      =   0 1 2
    #    gtf_entry.mAttributes by mFeature
    #               common Fields    =   exon_number
    #                                    gene_id
    #                                    gene_name
    #                                    transcript_id
    #                                    transcript_name
    #                CDS     mSource =   protein_coding
    #                                    IG_C_gene, IG_LV_gene, IG_V_gene, TR_V_gene
    #                                    non_stop_decay, nonsense_mediated_decay, polymorphic_pseudogene
    #                        Fields  =   gene_biotype
    #                                    protein_id
    #                exon    mSource =   everything
    #                        Fields  =   exon_id
    #                                    gene_biotype
    #                start_codon
    #                stop_codon
    #                        mSource =   protein_coding
    #                                    TR_V_gene,  IG_C_gene, IG_LV_gene
    #                                    non_stop_decay, nonsense_mediated_decay, polymorphic_pseudogene
    #                        Fields =    gene_biotype

    #_____________________________________________________________________________________
    #
    #   parse
    #_____________________________________________________________________________________
    def parse (self, gtf_file_path, logger):
        """
        parse file into regular data structures
        """
        if gtf_file_path[-3:] == '.gz':
            gtf_file = gzip.open(gtf_file_path)
        else:
            gtf_file = open(gtf_file_path)

        logger.info("  Parsing %s" % gtf_file_path)

        try:
            for line_num, gtf_entry in enumerate(minimal_gtf_iterator.iterator(gtf_file)):

                mAttributes = gtf_entry.mAttributes

                # common features

                exon_number = int(mAttributes["exon_number"]) - 1 if "exon_number" in mAttributes else None
                strand      = gtf_entry.mStrand in ['1', '+']
                gene_id     = gtf_entry.mGeneId, gtf_entry.mContig, strand
                cdna_id     = gtf_entry.mTranscriptId
                transcript_type= gtf_entry.mSource
                gene_type   = mAttributes.get("gene_biotype", "MISSING")
                beg         = int(gtf_entry.mStart) - 1
                end         = int(gtf_entry.mEnd  )
                interval    = (beg, end)


                #
                #   add data to self
                #
                self.feature_set.add(gtf_entry.mFeature)

                # transcript_name can be a list of strings
                if "transcript_name" in mAttributes:
                    if non_str_sequence(mAttributes["transcript_name"]):
                        self.cdna_id_to_names[cdna_id].mAttributes["transcript_name"]
                    else:
                        self.cdna_id_to_names[cdna_id].add(mAttributes["transcript_name"])


                # gene_name can be a list of strings
                if "gene_name" in mAttributes:
                    if non_str_sequence(mAttributes["gene_name"]):
                        self.gene_id_to_names[gene_id] |= mAttributes["gene_name"]
                    else:
                        self.gene_id_to_names[gene_id].add(mAttributes["gene_name"])


                self.gene_id_to_contig_strands.add(gene_id)
                self.gene_id_to_gene_type[gene_id] = gene_type

                self.cdna_id_to_transcript_type[cdna_id] =transcript_type

                self.gene_types.add(gene_type)
                self.gene_id_to_cdna_ids[gene_id].add(cdna_id)


                #
                #   exon
                #
                if gtf_entry.mFeature == "exon":
                    self.cdna_id_to_exons   [cdna_id].add((exon_number, beg, interval, mAttributes.get("exon_id", "MISSING")))

                #
                #   CDS
                #
                elif gtf_entry.mFeature == "CDS":
                    prot_id = mAttributes["protein_id"]
                    if (cdna_id in self.cdna_id_to_prot_id and
                        self.cdna_id_to_prot_id[cdna_id] != prot_id):
                        logger.warning("Transcript %s has two peptide IDs (%s and %s)"
                                       % (cdna_id, self.cdna_id_to_prot_id[cdna_id], prot_id))
                    self.cdna_id_to_prot_id[cdna_id] = prot_id
                    #
                    #   sort by start position
                    #   Assume each transcript does not have overlapping exons!!!
                    #
                    self.cdna_id_to_coding_exons [cdna_id].add((exon_number, beg, interval, mAttributes.get("exon_id", "MISSING"), int(gtf_entry.mFrame)))

                #
                # start /stop codon
                #
                #
                #   sort by exon number then start position
                #   Assume each transcript does not have overlapping exons!!!
                #
                elif gtf_entry.mFeature == "start_codon":
                    self.cdna_id_to_start_codons[cdna_id].add((exon_number, beg, interval))
                elif gtf_entry.mFeature == "stop_codon":
                    self.cdna_id_to_stop_codons[cdna_id].add((exon_number, beg, interval))
        except:
            #print line_num, "    " + dump_object(gtf_entry)
            print line_num, "    ", gtf_entry
            raise

        self.gene_types = list(self.gene_types)
        self.feature_set = list(self.feature_set)
        logger.info("  Features = %s" % (", ".join(self.feature_set)))
        logger.info("  Gene types = %s" % (", ".join(self.gene_types)))
        #assert(self.feature_set == ['start_codon', 'exon', 'stop_codon', 'CDS'])


    #_____________________________________________________________________________________
    #
    #   check_number_of_names
    #_____________________________________________________________________________________
    def check_number_of_names (self, logger):
        """
        Count the number genes / transcripts with > 1 name
        """
        for cdna_id, transcript_names in self.cdna_id_to_names.iteritems():
            if len(transcript_names) != 1:
                logger.warning("Transcripts with > 1 name: %23s %s" % (cdna_id, str(transcript_names)))

        for gene_id, gene_names in self.gene_id_to_names.iteritems():
            if len(gene_names) != 1:
                logger.warning("Genes with > 1 name: %23s %s" % (gene_id, str(gene_names)))

    #_____________________________________________________________________________________
    #
    #   check_split_codons
    #_____________________________________________________________________________________
    def check_split_codons (self, logger):
        """
        Count the number of split codons
        """
        def do_count (cdna_id_to_codons, description, logger):
            """
            abstract out counting for start/stop codons
            """
            cnt_split_codons  = 0
            for cdna_id, codons in cdna_id_to_codons.iteritems():
                if len(codons) > 1:
                    cnt_split_codons += 1
                if len(codons) > 2:
                    logger.error("Transcripts with %s codons " % description +
                                    "split over > 2 exons: "+
                                    "%23s %s" % (cdna_id, str(codons)))

            #logger.info("%3d transcripts with %s codons split over 2 exons: " %
            #                (cnt_split_codons, description))


        do_count (self.cdna_id_to_start_codons, "start", logger)
        do_count (self.cdna_id_to_stop_codons, "stop", logger)



#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#   Testing


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
if __name__ == '__main__':
    exe_path = os.path.split(os.path.abspath(sys.argv[0]))[0]
    #sys.path.append(os.path.abspath(os.path.join(exe_path,"gpipe")))

import unittest

class t_stderr_logger:
    """
    Everything to stderr
    """
    def __init__ (self):
        self.unique_prefix = ""
    def add_unique_prefix (self):
        import random
        random.seed()
        self.unique_prefix= str(random.randint(0,1000)) + " "
    def info (self, message):
        sys.stderr.write(self.unique_prefix + message + "\n")
    def warning (self, message):
        sys.stderr.write("\n\n" + self.unique_prefix + "WARNING:\n    " + message + "\n\n")
    def debug (self, message):
        sys.stderr.write(self.unique_prefix + message + "\n")

class t_stream_logger:
    """
    Everything to stderr
    """
    def __init__ (self, stream):
        self.stream = stream
    def info (self, message):
        self.stream.write(message + "\n")
    def warning (self, message):
        sys.stream.write("\n\nWARNING:\n    " + message + "\n\n")
    def debug (self, message):
        self.stream.write(message + "\n")


def run_function ():
    import cStringIO

    def iter (gtf_file, species, ignore_cache):

        ignore_strm = cStringIO.StringIO()
        output_strm = cStringIO.StringIO()
        ignore_logger = t_stream_logger(ignore_strm)
        logger = t_stream_logger(output_strm)
        gene_structures = t_parse_gtf(species)
        genes_by_type = gene_structures.get_genes(gtf_file, gtf_file + ".cache", ignore_logger, #["protein_coding"],
                                                                                                None, ignore_cache = ignore_cache)
        t_parse_gtf.log_gene_types (logger, genes_by_type)
        return genes_by_type, output_strm.getvalue()


    #
    # Use test data file
    #
    species = "test"
    exe_path = os.path.split(os.path.abspath(sys.argv[0]))[0]
    gtf_file = os.path.join(exe_path, "test_data", "test.shortish")

    #
    # Uncomment Use data files from ensembl
    #


    species  = "Mus musculus"
    gtf_file = "test_data/Mus_musculus.GRCm38.74.gtf"

    #gtf_file = "/data/mus/mirror/ftp.ensembl.org/pub/release-64/gtf/mus_musculus/Mus_musculus.NCBIM37.64.gtf.gz"

    #
    # determine heap size / memory usage
    #
    h = None
    try:
        #from guppy import hpy;
        #h=hpy()
        #h.setrelheap()
        pass
    except:
        pass

    genes1, output_str =iter(gtf_file, species, False)
    if h != None:
        heap_size = h.heap().size
        if heap_size > 1024:
            if heap_size > (1024 * 1024):
                heap_size_str = "%.1f MB" % (heap_size / 1024 / 1024)
            else:
                heap_size_str = "%.1f kB" % (heap_size / 1024)
        else:
            heap_size_str = "%d Bytes" % (heap_size)
        print ("\n\n" + "*" * 80 + "\n\n"   +
               "    Total memory needed for %s genes = %s\n\n" % (species, heap_size_str)  +
               "*" * 80 + "\n\n")

        #
        #   print heap details
        #
        #   print h.heap()

    correct_output_str = \
    """   Protein_coding genes.
              4 genes with ('nonsense_mediated_decay', 'processed_transcript', 'protein_coding', 'retained_intron') transcripts
              1 genes with ('nonsense_mediated_decay', 'protein_coding', 'retained_intron') transcripts
              2 genes with ('processed_transcript', 'protein_coding') transcripts
              9 genes with ('protein_coding',) transcripts
              1 genes with ('protein_coding', 'retained_intron') transcripts
           17 genes.
           80 transcripts.
          414 exons.
          145 exon regions.
            5 exon regions at median.
           29 exon regions at max.
            1 exon regions at min.
          113 exon coding regions.
            4 exon coding regions at median.
           22 exon coding regions at max.
            0 exon coding regions at min.
          544 exons in all transcripts.
          367 coding exons in all transcripts.
   Pseudogenes.
              5 genes with ('processed_pseudogene',) transcripts
              1 genes with ('unprocessed_pseudogene',) transcripts
            6 genes.
            6 transcripts.
            6 exons.
            6 exon regions.
            1 exon regions at median.
            1 exon regions at max.
            1 exon regions at min.
            0 exon coding regions.
            0 exon coding regions at median.
            0 exon coding regions at max.
            0 exon coding regions at min.
            6 exons in all transcripts.
         3 antisense            genes.
              3 genes with ('antisense',) transcripts
         1 lincRNA              genes.
              1 genes with ('lincRNA',) transcripts
         1 miRNA                genes.
              1 genes with ('miRNA',) transcripts
         4 snRNA                genes.
              4 genes with ('snRNA',) transcripts
         1 snoRNA               genes.
              1 genes with ('snoRNA',) transcripts
"""
    if species == "test":
        if output_str != correct_output_str:
            print "EXPECTING:\n", correct_output_str, "\n\n"
            print "FOUND:\n", output_str, "\n\n"
        assert (output_str == correct_output_str)
    else:
        print output_str


class Test_gene(unittest.TestCase):

    #       self.assertEqual(self.seq, range(10))
    #       self.assert_(element in self.seq)
    #       self.assertRaises(ValueError, random.sample, self.seq, 20)


    def test_function(self):
        #import pstats, cProfile
        #cProfile.runctx("run_function()", globals(), locals(), "Profile.prof")
        #s = pstats.Stats("Profile.prof")
        #s.strip_dirs().sort_stats("time").print_stats()
        run_function()




    def print_to_file (self, file_name, genes_per_gene_type):
        output_file = open(file_name, "w")
        gene_types = sorted(genes_per_gene_type.keys())
        for gene_type in gene_types:
            for g in genes_per_gene_type[gene_type]:
                output_file.write(str(g)+"\n")

    def iter(self, gtf_file, ignore_cache):
        """
            test
        """
        logger = t_stderr_logger()
        gene_structures = t_parse_gtf("Mus musculus")
        genes_by_type = gene_structures.get_genes(gtf_file, logger, ["protein_coding"], ignore_cache = ignore_cache)
        t_parse_gtf.log_gene_types (logger, genes_by_type)
        return genes_by_type




#
#   debug code not run if called as a module
#
if __name__ == '__main__':
    if sys.argv.count("--debug"):
        sys.argv.remove("--debug")
    unittest.main()





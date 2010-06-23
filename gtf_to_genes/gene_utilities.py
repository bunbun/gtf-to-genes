#!/usr/bin/env python
################################################################################
#
#   gene_utilities
#
#
#   Copyright (c) 20/6/2010 Leo Goodstadt
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
from collections import defaultdict

from gene import t_gene, t_transcript


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#   Functions


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#_____________________________________________________________________________________
#
#   index_transcripts
#_____________________________________________________________________________________
def index_transcripts(item, by_prot_id = True):
    transcript_index = dict()
    return do_index_transcripts(item, by_prot_id, transcript_index)

def do_index_transcripts(item, by_prot_id, transcript_index):
    """
    Go through item and find all t_gene objects and return a dictionary indexed by gene_id
    """
    if by_prot_id:
        attribute_name = "prot_id"
    else:
        attribute_name = "cdna_id"



    if isinstance(item, t_transcript):
        transcript_id = getattr(item,attribute_name)
        if transcript_id:
            transcript_index[transcript_id] = item
        return transcript_index

    if isinstance(item, t_gene):
        for t in item.transcripts:
            transcript_id = getattr(t,attribute_name)
            if transcript_id:
                transcript_index[transcript_id] = t
        return transcript_index

    # iterate through dictionary / list values
    if isinstance(item, dict):
        iteritem = item.itervalues()
    elif isinstance(item, (list, set, tuple)):
        iteritem = iter(item)
    else:
        return transcript_index

    for v in iteritem:
        if isinstance(v, t_gene):
            for t in v.transcripts:
                transcript_id = getattr(t,attribute_name)
                if transcript_id:
                    transcript_index[transcript_id] = t
        elif isinstance(v, t_transcript):
            transcript_id = getattr(v,attribute_name)
            if transcript_id:
                transcript_index[transcript_id] = v

        #   Don't know what this is : recurse
        #   Could special case string int here but if you are going to pass a big
        #       array of strings to this function and not expect some slow down...
        #
        else:
            do_index_transcripts(v, by_prot_id, transcript_index)

    return transcript_index

#_____________________________________________________________________________________
#
#   index_genes_by_gene_id
#_____________________________________________________________________________________
def index_genes_by_gene_id(item):
    gene_index = dict()
    return do_index_genes_by_gene_id(item, gene_index)

def do_index_genes_by_gene_id(item, gene_index):
    """
    Go through item and find all t_gene objects and return a dictionary indexed by gene_id
    """
    if isinstance(item, t_gene):
        gene_index[item.gene_id] = item
        return gene_index

    # iterate through dictionary / list values
    if isinstance(item, dict):
        itergenes = item.itervalues()
    elif isinstance(item, (list, set, tuple)):
        itergenes = iter(item)
    else:
        return gene_index

    for v in itergenes:
        if isinstance(v, t_gene):
            gene_index[v.gene_id] = v
        else:
            do_index_genes_by_gene_id(v, gene_index)

    return gene_index





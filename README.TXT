gtf_to_genes is a Python parser which caches all the genes / transcripts 
from a GTF file and caches the data into python classes for high speed 
access.

The initial parsing and indexing of a GTF file is a slow process, taking 
over a minute. But once that is done, the genes and transcripts can be 
retrieved in seconds.

To install run:

	python ./setup.py install

	or

	easy_install -U gtf_to-genes


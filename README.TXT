gtf_to_genes is a Python parser which caches all the genes / transcripts 
from a GTF file and caches the data into python classes for high speed 
access.

The initial parsing and indexing of a GTF file is a slow process, taking 
over a minute. But once that is done, the genes and transcripts can be 
retrieved in seconds.

Exons are stored efficiently per gene with each transcript only having indices into
these exon lists

array.array are used as much as possible for speed and efficiency.

To install run:

	python ./setup.py install

	or

	easy_install -U gtf_to-genes


TODO:

	* __str()__ should show exons rather than indices into parent's exons

# gtf\_to\_genes #
is a Python parser which caches all the genes / transcripts from a GTF file and caches the data into python classes for high speed access.



We want to access GTF format gene structure data with the maximum speed and with the minimum of fuss.

The initial parsing and indexing of a GTF file is a slow process, taking over a minute. But once the data can cached in binary format, the genes and transcripts can be retrieved efficiently in seconds..

Currently times to load the complete data set for Ensembl release 56 gene models are:

> Homo sapiens      139.0 Mb 9s

> Canis familiaris   52.0 Mb 3s

Homo sapiens is the largest GTF file available from Ensembl.

The numbers are for a duo core (the code is single threaded) computer with a SATA 7200rpm hard drives.

# Project Documentation #
Project documentation is available here:
http://gtf-to-genes.googlecode.com/hg/doc/build/html/index.html


# Installation #
Run:
```
	python ./setup.py install
```

or
```

	easy_install -U gtf_to-genes
```
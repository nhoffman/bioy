==========================================
bioy: a collection of bioinformatics tools
==========================================

Bio-y
    (pronounced "bio-ee") The adjective form of the noun "Bio"

.. contents:: Table of Contents

authors
=======

* Noah Hoffman
* Chris Rosenthal
* Tyler Land

dependencies
============

* A unix-like system; tested primarily on Ubuntu 12.04
* Python 2.7.x
* setuptools

Some functions require the following python packages:

* numpy
* pandas
* biopython

And other require external programs, including:

* ssearch36 (http://faculty.virginia.edu/wrpearson/fasta/fasta36/)
* usearch6 (http://www.drive5.com/usearch/download.html)
* muscle (http://www.drive5.com/muscle/)

installation
============

To install bioy and python dependencies, run setup.py or pip from the
project directory::

  % cd bioy
  % python setup.py install
  # or
  % pip install -U .

If you don't want to install the dependencies (numpy and pandas take a
while to compile), use::

  % pip install --no-deps -U .

Numpy and pandas require many dependencies to compile (and you'll
likely need to compile them because versions in package managers are
typically out of date). Fortunately, these can pretty easily be
installed on Ubuntu 12.04 by running::

  % sudo apt-get build-dep python-numpy python-pandas

A virtualenv containing a complete python execution environment can be
created using `dev/bootstrap.sh`::

  % dev/bootstrap.sh -h
  Create a virtualenv and install all pipeline dependencies
  Options:
  --venv            - path of virtualenv [bioy-env]
  --python          - path to the python interpreter [/usr/local/bin/python]
  --wheelstreet     - path to directory containing python wheels; wheel files will be
  in a subdirectory named according to the python interpreter version;
  uses WHEELSTREET if defined.
  (a suggested location is ~/wheelstreet) []
  --requirements    - a file listing python packages to install [requirements.txt]

execution
=========

The ``bioy`` script provides the user interface, and uses standard
UNIX command line syntax. Note that for development, it is convenient
to run ``bioy`` from within the project directory by specifying the
relative path to the script::

  % ./bioy

  Commands are constructed as follows. Every command starts with the
  name of the script, followed by an "action" followed by a series of
  required or optional "arguments". The name of the script, the action,
  and options and their arguments are entered on the command line
  separated by spaces. Help text is available for both the ``bioy``
  script and individual actions using the ``-h`` or ``--help`` options::

  usage: bioy [-h] [-V] [-v] [-q]

	      {help,align_clusters,all_pairwise,blast,children,classifier,classify,cmscores,consensus,csv2fasta,csv2hdf5,csvmod,dedup,denoise,errors,fasta,fasta2csv,fastq_stats,gb2fa,index,map_clusters,primer_trim,pull_reads,repl,reshape,reverse_complement,rldecode,rlencode,split_barcodes,split_reads,ssearch,ssearch2csv,ssearch_count,tree_edit,tsv2csv,usearch}
	      ...

  Tools for microbial sequence analysis and classification.

  positional arguments:
    {help,align_clusters,all_pairwise,blast,children,classifier,classify,cmscores,consensus,csv2fasta,csv2hdf5,csvmod,dedup,denoise,errors,fasta,fasta2csv,fastq_stats,gb2fa,index,map_clusters,primer_trim,pull_reads,repl,reshape,reverse_complement,rldecode,rlencode,split_barcodes,split_reads,ssearch,ssearch2csv,ssearch_count,tree_edit,tsv2csv,usearch}
      help                Detailed help for actions using `help <action>`
      align_clusters      Align reads contributing to a denoised cluster.
      all_pairwise        Calculate all Smith-Waterman pairwise distances among
			  sequences.
      blast               Run blastn and produce classify friendly output
      children            Return the children of a taxtable given a list of
			  taxids
      classifier          Classify sequences by grouping blast output by
			  matching taxonomic names
      classify            Classify sequences by grouping blast output by
			  matching taxonomic names
      cmscores            Convert raw cmalign alignment scores to csv format.
      consensus           Calculate the consensus for a multiple aignment
      csv2fasta           Turn a csv file into a fasta file specifying two
			  columns
      csv2hdf5            Convert a csv file to HDF5
      csvmod              Add or rename columns in a csv file.
      dedup               Fast deduplicate sequences by coalescing identical
			  substrings
      denoise             Denoise a fasta file of clustered sequences
      errors              Tally and classify errors given ./ion rlaligns
			  reference and query sequences
      fasta               Run the fasta pairwise aligment tool and output in csv
			  format.
      fasta2csv           Turn a fasta file into a csv
      fastq_stats         Describe distributions of sequencing quality scores
      gb2fa               Outputs a standard Genbank Record File into fasta file
			  format and optional seqinfo file in format ['seqname',
			  'tax_id','accession','description','length','ambig_cou
			  nt','is_type','rdp_lineage']
      index               Add simple indices to an sqlite database
      map_clusters        Create a readmap and specimenmap and/or weights file
			  from a
      primer_trim         Parse region between primers from fasta file
      pull_reads          Parse barcode, primer, and read from a fastq file
      repl                Replace strings in one or more files.
      reshape             convert a tsv file to a csv with an optional split/add
			  columns feature
      reverse_complement  reverse complement rle and non-rle sequences
      rldecode            Run-length decode a fasta file
      rlencode            Run-length encode a fasta file
      split_barcodes      Partition reads in a fastq file by barcode and write
			  an annotated fasta file
      split_reads         Parse reads from a fasta file by read to specimen csv
			  map file
      ssearch             Run the ssearch (Smith-Waterman) pairwise aligment
			  tool and output in csv format.
      ssearch2csv         Parse ssearch36 -m10 output and print specified
			  contents
      ssearch_count       Tally ssearch base count by position
      tree_edit           Tree leaf name editor that wraps BioPython.
      tsv2csv             convert a tsv file to a csv with an optional split/add
			  columns feature
      usearch             Run usearch global and produce classify friendly
			  output

  optional arguments:
    -h, --help            show this help message and exit
    -V, --version         Print the version number and exit
    -v, --verbose         Increase verbosity of screen output (eg, -v is
			  verbose, -vv more so)
    -q, --quiet           Suppress output

versions
========

We use abbrevited git sha hashes to identify the software version::

  % ./bioy --version
  0128.9790c13

The version information is saved in ``bioy_pkg/data`` when ``setup.py``
is run (on installation, or even by executing ``python setup.py
-h``).

unit tests
==========

Unit tests are implemented using the ``unittest`` module in the Python
standard library. The ``tests`` subdirectory is itself a Python
package that imports the local version (ie, the version in the project
directory, not the version installed to the system) of the
package. All unit tests can be run like this::

    % ./testall
    ...........
    ----------------------------------------------------------------------
    Ran 11 tests in 0.059s

    OK

A single unit test can be run by referring to a specific module,
class, or method within the ``tests`` package using dot notation::

    % ./testone -v tests.test_utils

license
=======

Copyright (c) 2012 Noah Hoffman

Released under the GPLv3 License:

TODO: include license text

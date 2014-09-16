==================
 Changes for bioy
==================


1.8.4-next release
============

* classify is now deprecated for the new faster classifier with test cases.  See bioy classifier --help for details.

* Known bugs: Tax_ids of valid Blast hits (hits that meet their rank thresholds) may be assigned 
              tax_ids of a higher threshold that *could* represent invalid tax_ids (tax_ids that may
              *not* have passed the rank threshold).

1.8.3
=====

 * bioy ssearch checks for empty files before calling ssearch, circumventing a race condition bug

1.8.2
=====

 * hdf5 support

1.8.0
=====

 * remove scipy dependency
 * add pandas dependency
 * dev/bootstrap.sh creates a development environment
 * further simplification of classify --details - only the top pident hit per assignment tax_id is reported
 * tests of subcommands use bioy_pkg.scrips.main as entry point


1.7.7
=====

 * bioy classify details output includes rows unique by sseqid,pident and coverage
 * new coverage column to blast output which is calculated in bioy blast
 * coverage filtering is done in bioy blast
 * bioy primer_trim has new option --keep-all-seqs which keeps sequences who's primers did not pass thresholds
 * bioy map_clusters takes a fasta file as an optional argument.  Optional arguments have been renamed and reorganized.
 * bugfix: formatting error in bioy map_clusters

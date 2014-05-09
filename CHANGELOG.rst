==================
 Changes for bioy
==================


upcoming
========

(features to be included in the next release)

 * remove scipy dependency
 * further simplification of details - only the top pident hit per assignment tax_id is reported

1.7.7
=====

 * bioy classify details output includes rows unique by sseqid,pident and coverage
 * new coverage column to blast output which is calculated in bioy blast
 * coverage filtering is done in bioy blast
 * bioy primer_trim has new option --keep-all-seqs which keeps sequences who's primers did not pass thresholds
 * bioy map_clusters takes a fasta file as an optional argument.  Optional arguments have been renamed and reorganized.
 * bugfix: formatting error in bioy map_clusters

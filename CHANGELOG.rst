==================
 Changes for bioy
==================

1.12-dev
========
 * new ``bioy ssearch_count`` output columns [tax_name, position, A, T, G, C, N, expected, naligns, nseqs, rank, id]
 * ``bioy classifier --copy-numbers`` uses column 16S instead of median

1.12
=======
 * ``bioy classifier --best-n-hits`` argument looks at each query sequence's Nth best blast hit, based on (GH 46)
   number of mismatches, and disregards all hits for that query with greater mismatches
 * ``bioy --split-condensed-assignments`` to retain assignment threshold assignments regardless of common 
   condense ids assignments (GH 53)
 * ``bioy classifier --include-ref-rank $RANK`` will pull in ${RANK}_id and ${RANK}_name columns from taxonomy into details.
 * Fixed ``bioy --split-condensed-assignments`` centroid bug (GH 58)

1.11
========
 * ``bioy --version`` outputs git tag release plus a count of additional commits above tag release (GH 44)
 * ``bioy classifier --hits-below-threshold`` switch will append hits below the threshold to ``--details-out`` file (GH 41)

1.10
=====
 * updated classifier to work with Python Pandas >=0.17.0
 * fixed classifier bug where assignment_id assignment was being unordered
 * classifier has new argument for emitting blast results normally discarded as below assignment threshold

1.9.4
=====
 * fixed a bug when duplicate values are in the specimen map cause duplicate entries in details (GH 38)

1.9.3
==========
 * added a workaround to get past a datatype inference bug in Pandas (GH 36)

1.9.2
==========
 * dynamically assigning rank thresholds and refreshed the available tax_ids in 
   data/rank_thresholds.csv (GH 32)
 * exchanged target_rank with condensed_rank (GH 31)
 * status messages now go to 100% (GH 26)
 * fixed bug when no blast hits (GH 34)

1.9.1
=====
 * bioy reverse_complement can take a raw dna string now. Fasta file input is specified using --is-file.
 * new csvjoin, csvconvert and csvdeduplicate that improve on some csvkit commands
 * new subcommand ncbi_fetch fetches sequences from NCBI, given sequence id

1.9
============

 * classify is now deprecated for the new faster classifier with test cases.  See bioy classifier --help for details.
 * classifier summary details now selects the largest clusters grouped by assignment_threshold to more closely align with the assignment pct_id range

 * Known bugs: Tax_ids of valid Blast hits (hits that meet their rank thresholds) may be assigned
              tax_ids of a higher threshold that *could* represent invalid tax_ids (tax_ids that may
              *not* have passed the rank threshold).
              Travis tests are failing even though they pass locally.  The Travis test environment is a work in progress.

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

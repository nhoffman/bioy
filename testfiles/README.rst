Input files for unit tests

2013-10-28: add testfiles/{F1_3,R1_36}
======================================

Take a subset of forward reads::

  bc=F1_3
  N=1000
  dest=~/src/bioy/testfiles/$bc
  cd /home/molmicro/analysis/043_ion_mcb/output
  mkdir -p $dest
  seqmagick convert --head $N $bc/primers/trimmed_rle.fasta $dest/trimmed_rle.fasta
  bzcat $bc/primers/trimmed_rle.csv.bz2 | head -$(expr $N + 1) | bzip2 > $dest/trimmed_rle.csv.bz2
  ~/src/bioy/bioy rldecode $dest/trimmed_rle.fasta $dest/trimmed_rle.csv.bz2 > $dest/trimmed.fasta
  usearch6 -cluster_fast $dest/trimmed.fasta -uc $dest/trimmed.uc -id 0.985 -quiet
  usearch6 -cluster_fast $dest/trimmed_rle.fasta -uc $dest/trimmed_rle.uc -id 0.985 -quiet

And some reverse reads - reverse complement these::

  bc=R1_36
  N=2000
  dest=~/src/bioy/testfiles/$bc
  mkdir -p $dest
  seqmagick convert --head $N $bc/primers/trimmed_rle.fasta $dest/temp_rle.fasta
  bzcat $bc/primers/trimmed_rle.csv.bz2 | head -$(expr $N + 1) | bzip2 > $dest/temp_rle.csv.bz2
  ~/src/bioy/bioy reverse_complement $dest/temp_rle.fasta $dest/temp_rle.csv.bz2 --out-fasta $dest/trimmed_rle.fasta --out-rle $dest/trimmed_rle.csv.bz2
  ~/src/bioy/bioy rldecode $dest/trimmed_rle.fasta $dest/trimmed_rle.csv.bz2 -o $dest/trimmed.fasta
  rm $dest/temp*
  usearch6 -cluster_fast $dest/trimmed.fasta -uc $dest/trimmed.uc -id 0.985 -quiet
  usearch6 -cluster_fast $dest/trimmed_rle.fasta -uc $dest/trimmed_rle.uc -id 0.985 -quiet

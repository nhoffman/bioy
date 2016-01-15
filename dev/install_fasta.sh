#!/bin/bash

# Install fasta tools - installs to the first element of $PATH by default.

set -e
set -v

if [[ -z $1 ]]; then
    echo "Usage $(basename $0) bindir [version] [srcdir]"
    exit 1
fi

# if [[ -z $1 ]]; then
#     bindir=$(echo $PATH | tr ':' '\n' | head -1)
# else
#     bindir=$1
# fi

_here=$(pwd)

bindir=${1-/usr/local/bin}
version=${2-36.3.6d}
srcdir=${3-src}

if $prefix/bin/fasta36 -help 2>&1 | grep -q ${version:0:-1}; then
    echo "fasta version $version is already installed in $prefix/bin"
    exit 0
fi

mkdir -p $srcdir
cd $srcdir
wget -N http://faculty.virginia.edu/wrpearson/fasta/fasta36/fasta-${version}.tar.gz
tar -xf fasta-${version}.tar.gz
cd fasta-${version}/src
# make -f ../make/Makefile.mpi_icc_sse2 all
make -f ../make/Makefile.linux64_sse2 all

# return starting dir in case $bindir is a relative path
cd $_here
cp $srcdir/fasta-${version}/bin/* $bindir

# confirm success
$prefix/bin/fasta36 -help 2>&1 | grep ${version:0:-1}

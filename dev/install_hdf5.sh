#!/bin/bash
# Usage: install_hdf5.sh [--prefix <prefix>]

# Install the hdf5 binaries and shared libraries

set -e

# options configurable from the command line
GREP_OPTIONS=--color=never
PREFIX=/usr/local
SRC=src
HDF5_VERSION=1.8.13

if [[ $1 == '-h' || $1 == '--help' ]]; then
    echo "Install the hdf5 binaries and shared libraries"
    echo "Options:"
    echo "--prefix    - path to directory containing bin/ lib/ etc [$PREFIX]"
    echo "--src       - path to source directory [$SRC]"
    echo "--version   - hdf5 library version [$HDF5_VERSION]"
    exit 0
fi

while true; do
    case "$1" in
	--prefix ) PREFIX="$(readlink -f $2)"; shift 2 ;;
	--src ) SRC="$2"; shift 2 ;;
	--version ) HDF5_VERSION="$2"; shift 2 ;;
	* ) break ;;
    esac
done

mkdir -p "$SRC"

hdf5=hdf5-${HDF5_VERSION}-linux-x86_64-shared

cd $SRC
test -f ${hdf5}.tar.gz || \
    wget http://www.hdfgroup.org/ftp/HDF5/current/bin/linux-x86_64/${hdf5}.tar.gz

echo "Installing to $PREFIX"

tar -xf ${hdf5}.tar.gz
for subdir in include lib share bin; do
    mkdir -p "$PREFIX/$subdir"
    cp -r $hdf5/$subdir/* "$PREFIX/$subdir"
done


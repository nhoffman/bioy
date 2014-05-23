#!/bin/bash
# Usage: bin/bootstrap.sh [Options]

# Create a virtualenv, and install requirements to it.

# override the default python interpreter using
# `PYTHON=/path/to/python bin/bootstrap.sh`

if [[ -n $VIRTUAL_ENV ]]; then
    echo "You can't run this script inside an active virtualenv"
    exit 1
fi

set -e

# options configurable from the command line
GREP_OPTIONS=--color=never
VENV=$(basename $(pwd))-env
PYTHON=$(which python)
PY_VERSION=$($PYTHON -c 'import sys; print "{}.{}.{}".format(*sys.version_info[:3])')
# WHEELSTREET=$(readlink -f ~/wheelstreet)
REQFILE=requirements.txt

if [[ $1 == '-h' || $1 == '--help' ]]; then
    echo "Create a virtualenv and install all pipeline dependencies"
    echo "Options:"
    echo "--venv            - path of virtualenv [$VENV]"
    echo "--python          - path to the python interpreter [$PYTHON]"
    echo "--wheelstreet     - path to directory containing python wheels; wheel files will be"
    echo "                    in a subdirectory named according to the python interpreter version;"
    echo "                    uses WHEELSTREET if defined."
    echo "                    (a suggested location is ~/wheelstreet) [$WHEELSTREET]"
    echo "--requirements    - a file listing python packages to install [$REQFILE]"
    exit 0
fi

while true; do
    case "$1" in
	--venv ) venv="$2"; shift 2 ;;
	--python ) PYTHON="$2"; shift 2 ;;
	--wheelstreet ) WHEELSTREET=$(readlink -f "$2"); shift 2 ;;
	--requirements ) REQFILE="$2"; shift 2 ;;
	* ) break ;;
    esac
done

VENV_VERSION=1.11.4
WHEELHOUSE=

if [[ -n $WHEELSTREET ]]; then
    WHEELHOUSE=$WHEELSTREET/$PY_VERSION
    test -d $WHEELHOUSE || (echo "cannot access $WHEELHOUSE"; exit 1)
fi

mkdir -p src

# Create the virtualenv using a specified version of the virtualenv
# source. This also provides setuptools and pip. Inspired by
# http://eli.thegreenplace.net/2013/04/20/bootstrapping-virtualenv/

# create virtualenv if necessary
if [ ! -f ${VENV:?}/bin/activate ]; then
    # download virtualenv source if necessary
    if [ ! -f src/virtualenv-${VENV_VERSION}/virtualenv.py ]; then
	VENV_URL='https://pypi.python.org/packages/source/v/virtualenv'
	(cd src && \
	    wget -N ${VENV_URL}/virtualenv-${VENV_VERSION}.tar.gz && \
	    tar -xf virtualenv-${VENV_VERSION}.tar.gz)
    fi
    $PYTHON src/virtualenv-${VENV_VERSION}/virtualenv.py $VENV
    # $PYTHON src/virtualenv-${VENV_VERSION}/virtualenv.py --relocatable $VENV
else
    echo "found existing virtualenv $VENV"
fi

# install required libraries
# hdf5: tables (PyTables)

# export HDF5_DIR=$(readlink -f $VENV)
# dev/install_hdf5.sh --prefix $HDF5_DIR --src $(readlink -f src)

# # make the above libraries available when the virtualenv is activated
# if grep -q LD_LIBRARY_PATH $VENV/bin/activate; then
#     true
# else
#     cat >> $VENV/bin/activate <<EOF

# # added by $0
# LD_LIBRARY_PATH="$VIRTUAL_ENV/lib"
# export LD_LIBRARY_PATH
# EOF
# fi

source $VENV/bin/activate

# install python packages from pipy or wheels
grep -v -E '^#|git+|^-e' $REQFILE | while read pkg; do
    if [[ -z $WHEELHOUSE ]]; then
	pip install --allow-external argparse $pkg
    else
	pip install --allow-external argparse --use-wheel --find-links=$WHEELHOUSE $pkg
    fi
done

# install packages from git repos if necessary
if [[ ! -z $(grep git+ $REQFILE | grep -v -E '^#') ]]; then
    pip install -r <(grep git+ $REQFILE | grep -v -E '^#')
else
    echo "no packages to install from git repositories"
fi

# correct any more shebang lines
# virtualenv --relocatable $VENV

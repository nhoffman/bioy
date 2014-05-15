#!/bin/bash

set -e

GREP_OPTIONS=--color=never

WHEELSTREET=$(readlink -f ~/wheelstreet)
REQFILE=requirements.txt
PYTHON=$(which python)

if [[ $1 == '-h' || $1 == '--help' ]]; then
    echo "Create a virtualenv and install all pipeline dependencies"
    echo "Options:"
    echo "--python          - path to an alternative python interpreter [$PYTHON]"
    echo "--wheelstreet     - path to directory containing python wheels; wheel files will be"
    echo "                    in a subdirectory named according to the python interpreter version"
    echo "                    (eg '2.7.5') [$WHEELSTREET]"
    echo "--requirements    - an alternative requiremets file [$REQFILE]"
    exit 0
fi

while true; do
    case "$1" in
	--python ) PYTHON="$2"; shift 2 ;;
	--wheelstreet ) WHEELSTREET=$(readlink -f "$2"); shift 2 ;;
	--requirements ) REQFILE="$2"; shift 2 ;;
	* ) break ;;
    esac
done

if [[ ! -f "$REQFILE" ]]; then
    echo "Cannot find requirements file named $REQFILE"
    exit 1
fi

# get the python version
TAG=$($PYTHON -c 'import sys; print "{}.{}.{}".format(*sys.version_info[:3])')
PY_VER="${TAG:0:3}"
PY_MAJOR="${PY_VER:0:1}"
VENV_VERSION=1.11.4

WHEELHOUSE="$WHEELSTREET/$TAG"  # wheels for this python version
CACHE="$WHEELSTREET/cache/$TAG"
VENV="$WHEELSTREET/venv"
SRC="$WHEELSTREET/src"

mkdir -p "$CACHE"
mkdir -p "$WHEELHOUSE"

# create virtualenv if necessary
if [ ! -f $VENV/bin/activate ]; then
    # download virtualenv source if necessary
    if [ ! -f $CACHE/virtualenv-${VENV_VERSION}/virtualenv.py ]; then
	VENV_URL='https://pypi.python.org/packages/source/v/virtualenv'
	(cd $CACHE && \
	    wget -N ${VENV_URL}/virtualenv-${VENV_VERSION}.tar.gz && \
	    tar -xf virtualenv-${VENV_VERSION}.tar.gz)
    fi
    $PYTHON $CACHE/virtualenv-${VENV_VERSION}/virtualenv.py $VENV
    # $PYTHON $CACHE/virtualenv-${VENV_VERSION}/virtualenv.py --relocatable $VENV
else
    echo "found existing virtualenv $VENV"
fi

source $VENV/bin/activate

# install wheel package to $VENV
pip install --download-cache $CACHE wheel

# put the requirements file in the wheelhouse
cp $REQFILE $WHEELHOUSE

# install and build the wheels; if a package is installed from pipy in
# the form of a wheel, copy it from the cache into the wheelhouse.

wheelname(){
    python -c "print \"$1\".split('%2F')[-1]"
}

get_wheels_from_cache(){
    echo "Copying wheels from $CACHE"
    if [[ -f $CACHE/*.whl ]]; then
	for whl in $CACHE/*.whl; do
	    mv $whl $WHEELHOUSE/$(wheelname $whl)
	done
    fi
}

grep -v -E '^#|git\+' $REQFILE | while read pkg; do
    # --find-links avoids rebuilding existing wheels
    pip wheel $pkg \
    	--download-cache $CACHE \
    	--use-wheel \
    	--find-links=$WHEELHOUSE \
    	--wheel-dir=$WHEELHOUSE

    get_wheels_from_cache

    # install to provide dependencies for subsequent packages
    pip install --use-wheel --find-links=$WHEELHOUSE --no-index $pkg
done

get_wheels_from_cache

echo 'Created wheels:'
tree -F -L 2 $WHEELSTREET

echo 'install wheels using'
echo "pip install --use-wheel --find-links=$WHEELHOUSE --no-index -r $REQFILE"

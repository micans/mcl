#!/bin/bash

set -euo pipefail

# This is a quite clumsy demo script with recipes for installing mcl.
# Requirement: some directory with two tarballs, one for cimfomfa, one for mcl.
# Release versions are hardcoded in this script.
# Usage:
# ./installmcl cm THETARBALLDIRECTORY     # install both cimfomfa and mcl
# Or
# ./installmcl c THETARBALLDIRECTORY      # first install cimfomfa
# ./installmcl m THETARBALLDIRECTORY      # then install mcl

# Change versions below if you have different versions.
cff=21-101
mcl=21-257

# Change if you want to install somewhere else
INSTALL=$HOME/phloobah


# Now the rest of this script should have enough to run.



modes=${1?Need mode; c for cimfomfa and/or m for mcl}
SOURCE=${2?Please supply directory to grab tarballs from}

mcltar=mcl-$mcl.tar.gz
cfftar=cimfomfa-$cff.tar.gz

if [[ $modes =~ c ]]; then
  thedir=./${cfftar%.tar.gz}
  rm -rf $thedir
  cp $SOURCE/$cfftar . && tar xzf $cfftar
  ( cd $thedir
    ./configure --prefix=$INSTALL
    make
    make install
  )
fi

if [[ $modes =~ m ]]; then
  cp $SOURCE/$mcltar . && tar xzf $mcltar
  thedir=./${mcltar%.tar.gz}
  ( cd $thedir
    ./configure CFLAGS=-I$INSTALL/include LDFLAGS=-L$INSTALL/lib --prefix=$INSTALL
    make
    make install
  )
fi



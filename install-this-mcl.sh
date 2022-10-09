#!/bin/bash

set -euo pipefail

# This script is a simple script for downloading + compiling mcl
cff=22-273
mcl=22-282

# Change if you want to install somewhere else
INSTALL=$HOME/local

# Now the rest of this script should have enough to run.
mcltar=mcl-$mcl.tar.gz
cfftar=cimfomfa-$cff.tar.gz

if command -v wget > /dev/null; then 
   webbit=wget
elif command -v curl > /dev/null; then 
   webbit="curl -O"
else
   echo "Explain to me how to download stuff please"
   false
fi

$webbit http://micans.org/mcl/src/$mcltar
$webbit http://micans.org/mcl/src/$cfftar

if true; then
  thedir=./${cfftar%.tar.gz}
  rm -rf $thedir
  tar xzf $cfftar
  ( cd $thedir
    ./configure --prefix=$INSTALL --disable-shared
    make
    make install
  )
fi

if true; then
  tar xzf $mcltar
  thedir=./${mcltar%.tar.gz}
  ( cd $thedir
    ./configure CFLAGS=-I$INSTALL/include LDFLAGS=-L$INSTALL/lib --prefix=$INSTALL --enable-rcl
    make
    make install
  )
fi


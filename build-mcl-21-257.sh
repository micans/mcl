#!/bin/bash

set -euo pipefail

# Change versions below if you have different versions.
cff=21-257
mcl=21-257

# Change if you want to install somewhere else
INSTALL=$HOME/local


# Now the rest of this script should have enough to run.

mcltar=mcl-$mcl.tar.gz
cfftar=cimfomfa-$cff.tar.gz

wget http://micans.org/mcl/dev/$mcltar
wget http://micans.org/mcl/dev/$cfftar

if true; then
  thedir=./${cfftar%.tar.gz}
  rm -rf $thedir
  tar xzf $cfftar
  ( cd $thedir
    ./configure --prefix=$INSTALL
    make
    make install
  )
fi

if true; then
  tar xzf $mcltar
  thedir=./${mcltar%.tar.gz}
  ( cd $thedir
    ./configure CFLAGS=-I$INSTALL/include LDFLAGS=-L$INSTALL/lib --prefix=$INSTALL
    make
    make install
  )
fi


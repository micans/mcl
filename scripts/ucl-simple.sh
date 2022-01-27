#!/bin/bash

set -euo pipefail

   # Given a set of input clusterings (everything but the first argument),
   # construct the matrix with a contribution of 1 for the edge (i,j) each time
   # i and j co-occur in a clustering.  This matrix is restricted to
   # the edges present in the network itself (the first argument).
   # The output is normalised to a scale of [0-1000].

   # mcxi understands a simple stack language with some primitive
   # facilities for matrix manipulations in the mcl framework.
   # A useful thing about it is it can handle very large matrices
   # with millions of nodes as well as mcl can.
   # This script can be used to compare consensus clustering methods (such as
   # rcl) with the very simplistic approach implemented here.
   #
matrix=${1?Need <mxfile> <clsfile>+}
shift 1

   # prepend a slash to each input cluster file,
   # this is how mcxi expects strings.
   #
clusters=("${@/#/\/}")
ncluster=${#clusters[*]}

   # Below, project onto matrix (=network) using
   # hadamard product in the loop.
   # lm  load matrix
   # tp  transpose
   # hdm hadamard product
   # 1444 is a sentinel to check end of loop with 'type'
   # .x dim pop just checks whether this is a network and not a clustering,
   # it will complain and exit if src/dst dimensions are not the same.
   # 
mcxi <<EOC
0 vb
/$matrix lm ch dup .x def .sum def
   .x dim pop
   1444 ${clusters[@]}
   { type /str eq }
   { lm tp mul .x hdm .sum add .sum def }
   while
   .sum 1000.0 $ncluster.0 1 add div mul
   /out.ucl wm
EOC

echo done


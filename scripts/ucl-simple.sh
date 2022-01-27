#!/bin/bash

   # ucl - Unrestricted Contingency Linkage
   # A lowest common denominator approach to consensus clustering;
   # this just creates a graph based on (i,j) cluster co-occurrences.
   # It can be used e.g. to compare, as a baseline, to
   # the rcl approach (restricted contigency linkage).

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

set -euo pipefail

matrix=${1?Need <mxfile> <clsfile>+}
shift 1

   # prepend a slash to each input cluster file,
   # this is how mcxi expects strings.
   #
clusters=("${@/#/\/}")
ncluster=${#clusters[*]}


   # Below, project onto matrix (=network) using
   # hadamard product in the loop.
   # lm  : load matrix
   # ch  : make characteristic (all entries set to 1)
   # tp  : transpose
   # hdm : hadamard product
   # .x dim pop just checks whether this is a network and not a clustering,
   # it will complain and exit if src/dst dimensions are not the same.
   # 'type' returns null if the stack is exhausted.
   # 
mcxi <<EOC
0 vb
/$matrix lm ch dup .x def .sum def
   .x dim pop
   ${clusters[@]}
   { type /str eq }
   { lm tp mul .x hdm .sum add .sum def }
   while
   .sum 1000.0 $ncluster.0 1 add div mul
   /out.ucl wm
EOC

echo $? done


#!/bin/bash

set -euo pipefail

matrixfile=
rclfile=
tabfile=
ucl=false
LEVELS=
RESOLUTION=

while getopts :m:n:t:l:r:Uh opt
do
    case "$opt" in
    m)
      matrixfile=$OPTARG
      ;;
    n)
      rclfile=$OPTARG
      ;;
    r)
      RESOLUTION=$OPTARG
      ;;
    l)
      LEVELS=$OPTARG
      ;;
    t)
      tabfile=$OPTARG
      ;;
    U)
      ucl=true
      ;;
    h)
      cat <<EOU
Options described below. Below 'networks', 'clusterings' and 'matrix' are files
in mcl matrix format. Use srt2cls.sh to convert seurat clusterings to mcl format.
Suggested usage:
1) compute the rcl object with

      rcl.sh -m <network> -t <tab> -n <rcl-file> [-U] <LIST-OF-CLUSTER-FILE-NAMES>

2) gain an understanding of the dynamic range of the clusters present in the tree
   encoded in the rcl object by computing cluster sizes of thresholded trees

      rcl.sh -n <rcl-file> -l <LOW/STEP/HIGH>

   e.g. -l 200/100/600
   Note that the edge weight range in the rcl objects is [0-1000].

3) derive resolution-based balanced clusterings from the rcl object.

      rcl.sh -n <rcl-file> -r "N1 N2 N3 .."

   e.g. -r "500 1000 1500 2000 2500"
   The largest clusters obtained will be roughly near the resolution
   in that they cannot be split into smaller clusters at least that size

Options:
-m  <file>   Input network/matrix file
-t  <file>   Tab file with index - label mapping, format INDEX<TAB>LABEL
-n  <file>   Name for rcl/ucl object
-U           Compute the Unrestricted Contingency Linkage object
-l  LOW/STEP/HIGH    e.g. 200/100/600 to show threshold cluster sizes
-r  "N1 N2 N3 .."    e.g. "500 1000 1500 2000 2500" to compute resolution clusterings
EOU
      exit
      ;;
    :) echo "Flag $OPTARG needs argument"
        exit 1;;
    ?) echo "Flag $OPTARG unknown"
        exit 1;;
   esac
done


if [[ -z $rclfile ]]; then
   rclfile=rcl
   if $ucl; then rclfile=ucl; fi
fi

if [[ -z $matrixfile && ! -f $rclfile ]]; then
   echo "Please supply network and clusterings to compute rcl object (see -h)"
   false
fi
if [[ -z $tabfile && ! -f $rclfile.tab ]]; then
   echo "Please supply the tab file mapping indexes to labels (-t)"
   false
fi

if [[ ! -f $rclfile.tab ]]; then
   cp $tabfile $rclfile.tab
fi


if [[ ! -f $rclfile ]]; then

   shift $(($OPTIND-1))
   if (( $# < 2 )); then
      echo "Please supply a few clusterings"
      false
   fi

   if $ucl; then

      # prepend a slash to each input cluster file,
      # this is how mcxi expects strings. Much love to bash syntax.
      # No spaces in file names please.
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

      echo "-- Computing UCL object"
      mcxi <<EOC
0 vb
/$matrixfile lm ch dup .x def .sum def
   .x dim pop
   ${clusters[@]}
   { type /str eq }
   { lm tp mul .x hdm .sum add .sum def }
   while
   .sum 1000.0 $ncluster.0 1 add div mul
   /$rclfile wm
EOC

   else
      echo "-- Computing RCL object"
      clm vol --progress -imx $matrixfile -write-rcl $rclfile -o $rclfile.vol "$@"
   fi

   echo "-- Computing single linkage join order for network $rclfile"
   clm close --sl -imx $rclfile -tab $rclfile.tab -o $rclfile.join-order -write-sl-list $rclfile.node-values

else
   echo "-- $rclfile exists, querying its dimensions"
   if ! mcx query -imx $rclfile --dim; then
      echo "Suggest removing $rclfile"
   fi
   if [[ -z $LEVELS && -z $RESOLUTION ]]; then
      echo "Suggest using -l LOW/STEP/HIGH e.g. -l 200/100/600 (scale 0-1000)"
      echo "or -r \"N1 N2 N3 ..\" to compute resolution clusters, e.g. -r \"500 1000 1500 2000 2500\""
   fi
fi



if [[ ! -z $LEVELS ]]; then
   echo "-- cluster sizes resulting from simple thresholding with levels $LEVELS"
   clm close -imx $rclfile -levels $LEVELS
fi

   # From the tree, compute a balanced clustering according to a
   # resolution parameter. Given resolution R, the script will
   # continue to descend the tree as long as there is a split
   # in the nodes below into two clusters each of size >= R.
   # The info files contain
   # 1) size of clusters
   # 2) a measure of cohesiveness of the cluster (higher is more cohesive)
   # 3) the node indexes for each cluster
   # Then transform the info files into clusterings that mcl
   # and siblings understand; these are then dumped with labels
   # in a line-based format.
   #

if [[ ! -z $RESOLUTION ]]; then
   echo "-- computing balanced clusterings with resolution parameters $RESOLUTION"
   export MCLXIOVERBOSITY=2
   for r in $RESOLUTION; do
      prefix="$rclfile.mix$r"
                            # this sort orders largest clusters first.
      rcl-mix.pl $r $rclfile.join-order | sort -nr > $prefix.info
      cut -f 3 $prefix.info | mcxload -235-ai - -o $prefix.cls
      mcxdump -icl $prefix.cls -tabr $rclfile.tab -o $prefix.labels
      mcxdump -imx $prefix.cls -tabr $rclfile.tab --no-values --transpose -o $prefix.txt
   done
   for r in $RESOLUTION; do
      file="$rclfile.mix$r.cls"
      printf "%-15s" $file; echo $(mcx query -imx $file | cut -f 2 | tail -n +2)
   done
fi



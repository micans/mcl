#!/bin/bash

# RCL - Restricted Contingency Linkage - consensus clustering.
# See github.com/micans/mcl and -h (or -H).

set -euo pipefail

matrixfile=
tabfile=
pfx=
LEVELS=
RESOLUTION=

do_ucl=false
do_gralog=false

while getopts :m:n:t:l:r:UShH opt
do
    case "$opt" in
    m)
      matrixfile=$OPTARG
      ;;
    n)
      pfx=$OPTARG
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
    S)
      do_gralog=true
      ;;
    U)
      do_ucl=true
      ;;
    h)
      cat <<EOU
1) compute the rcl object with
      rcl.sh -n NAME -m <network> -t <tabfile> <LIST-OF-CLUSTER-FILE-NAMES>
   NAME will be used as a prefix for various outputs; think of it as a project tag.

2) derive resolution-based clusterings from the rcl object.
      rcl.sh -n NAME [-S] -r "N1 N2 N3 .."
      e.g. -r "500 1000 1500 2000 2500"

Options:
-m  <file>   Input network/matrix file in mcl format (obtain e.g. with mcxload)
-t  <file>   Tab file with index - label mapping (obtain e.g. with mcxload)
-n  NAME     NAME will be used as prefix for various objects
-l  LOW/STEP/HIGH    e.g. 200/50/700 to show threshold cluster sizes
-r  "N1 N2 N3 .."    e.g. "500 1000 1500 2000 2500" to compute resolution clusterings
-S           Output cluster size distribution information (use in (2))
-U           Compute the Unrestricted Contingency Linkage object (use in (1))

Use -H to output a slightly longer description.
EOU
       exit
      ;;
    H)
      cat <<EOU | less
Options described below. Below 'networks', 'clusterings' and 'matrix' are files
in mcl matrix format. Networks can be loaded from label data with mcxload.  Use
srt2cls.sh to convert seurat clusterings to mcl format.  Both will need a label
index 'tab' file, this can be created with srt2tab.sh.  Note that to maintain
index correspondence with Seurat's 1-based indexing, srt2cls.sh and srt2tab.sh
introduce a dummy node in the mcl representations for index 0.

Suggested usage:
1) compute the rcl object with
      rcl.sh -n NAME -m <network> -t <tab> <LIST-OF-CLUSTER-FILE-NAMES>
   NAME will be used as a prefix for various outputs; think of it as a project tag.
   NAME is used in 2) to retrieve the right objects.

2) derive resolution-based clusterings from the rcl object.
      rcl.sh -n NAME [-S] -r "N1 N2 N3 .."
      e.g. -r "500 1000 1500 2000 2500"
   The largest clusters obtained will be above the resolution limit in that
   there is no sub-split into smaller clusters at least that size.  Cluster
   sizes can be below the resolution limit as such clusters may need to be
   split off in order to allow another allowable split to happen.

Optionally:
3) Investigate in more detail the dynamic range of the clusters present in the tree
   encoded in the rcl object by computing cluster sizes of thresholded trees
      rcl.sh -n NAME -l <LOW/STEP/HIGH>
      e.g. -l 200/50/700
   Note that the edge weight range in the rcl objects is [0-1000].
   From the output you may wish to zoom in
      e.g. -l 450/10/550
   if the cluster sizes in that range are what you are after.
   To save a bunch of such clusterings, use e.g.
      rcl.sh -n NAME -l 470/10/530/pfx

Options:
-m  <file>   Input network/matrix file in mcl format (obtain e.g. with mcxload)
-t  <file>   Tab file with index - label mapping (obtain e.g. with mcxload)
-n  NAME     NAME will be used as prefix for various objects
-l  LOW/STEP/HIGH    e.g. 200/50/700 to show threshold cluster sizes
-r  "N1 N2 N3 .."    e.g. "500 1000 1500 2000 2500" to compute resolution clusterings
-S           Output cluster size distribution information (use in (2))
-U           Compute the Unrestricted Contingency Linkage object (use in (1))
EOU
      exit
      ;;
    :) echo "Flag $OPTARG needs argument"
        exit 1;;
    ?) echo "Flag $OPTARG unknown"
        exit 1;;
   esac
done


if [[ -z $pfx ]]; then
   echo "Please specify -n NAME to tag this analysis (see -h)"
   false
fi

rclfile=$pfx.rcl

if [[ -z $matrixfile && ! -f $rclfile ]]; then
   echo "Please supply network and clusterings to compute rcl object (see -h)"
   false
fi
if [[ -z $tabfile && ! -f $pfx.tab ]]; then
   echo "Please supply the tab file mapping indexes to labels (-t)"
   false
fi

if [[ ! -f $pfx.tab ]]; then
   if ! ln $tabfile $pfx.tab 2> /dev/null; then
     cp $tabfile $pfx.tab
   fi
fi


if [[ ! -f $rclfile ]]; then

   export MCLXIOFORMAT=8

   shift $(($OPTIND-1))
   if (( $# < 2 )); then
      echo "Please supply a few clusterings"
      false
   fi

   if $do_ucl; then

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
      clm vol --progress -imx $matrixfile -write-rcl $rclfile -o $pfx.vol "$@"
   fi

   echo "-- Computing single linkage join order for network $rclfile"
   clm close --sl -imx $rclfile -tab $pfx.tab -o $pfx.join-order -write-sl-list $pfx.node-values

else
   echo "-- $rclfile exists, querying its dimensions"
   if ! mcx query -imx $rclfile --dim; then
      echo "Suggest removing $rclfile"
   fi
fi

if [[ -z $LEVELS && -z $RESOLUTION ]]; then
  echo "-- suggest rcl.sh -n $pfx -r \"N1 N2 N3 ..\" to compute resolution clusters"
  echo "-- e.g. rcl.sh -n $pfx -r \"500 1000 1500 2000 2500\""
  echo "-- vary N according to preference and data set size"
fi

export MCLXIOFORMAT=1

if [[ ! -z $LEVELS ]]; then
   echo "-- cluster sizes resulting from simple thresholding with levels $LEVELS"
   echo "-- (output may be truncated)"
   clm close -imx $rclfile -levels $LEVELS | cut -b 1-100
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
   echo "-- computing clusterings with resolution parameters $RESOLUTION"
   export MCLXIOVERBOSITY=2

   rcl-res.pl $pfx $RESOLUTION < $pfx.join-order
   echo "-- saving resolution cluster files"
   echo "-- displaying size of the 20 largest clusters"
   if $do_gralog; then echo "-- summarising cluster size distribution on a log scale"; fi

                                       # To help space the granularity output.
   export CLXDO_WIDTH=$((${#pfx}+14))  # .res .cls length 8, leave 6 for resolution

   for r in $RESOLUTION; do
      rfile=$pfx.res$r.info
      if [[ ! -f $rfile ]]; then
         echo "Expected file $rfile not found"
         false
      fi
      prefix="$pfx.res$r"
      cut -f 4 $rfile | mcxload -235-ai - -o $prefix.cls
      mcxdump -icl $prefix.cls -tabr $pfx.tab -o $prefix.labels
      mcxdump -imx $prefix.cls -tabr $pfx.tab --no-values --transpose -o $prefix.txt
      clxdo grabig 20 $prefix.cls
      if $do_gralog; then clxdo gralog $prefix.cls; echo "()"; fi
   done
   commalist=$(tr -s ' ' ',' <<< $RESOLUTION)

cat <<EOM

The following outputs were made.
One cluster-per line files with labels:
   $(eval echo $pfx.res{$commalist}.labels)
LABEL<TAB>CLUSID files:
   $(eval echo $pfx.res{$commalist}.txt)
mcl-edge matrix/cluster files (suitable input e.g. for 'clm dist'):
   $(eval echo $pfx.res{$commalist}.cls)
EOM

fi


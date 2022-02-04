#!/bin/bash

# RCL - Restricted Contingency Linkage - consensus clustering.
# See github.com/micans/mcl and -h (or -H).

set -euo pipefail

matrixfile=
tabfile=
pfx=
LEVELS=
RESOLUTION=
cpu=1

do_ucl=false
do_gralog=false
do_force=false


HELP_intro="
Options described below. Below 'networks', 'clusterings' and 'matrix' are files
in mcl matrix format. Networks can be loaded from label data with mcxload.  Use
srt2cls.sh to convert seurat clusterings to mcl format.  Both will need a label
index 'tab' file, this can be created with srt2tab.sh.  Note that to maintain
index correspondence with Seurat's 1-based indexing, srt2cls.sh and srt2tab.sh
introduce a dummy node in the mcl representations for index 0."

Help_1='
1) compute the rcl object with
      rcl.sh -n NAME -m <network> -t <tabfile> <LIST-OF-CLUSTER-FILE-NAMES>
   NAME will be used as a prefix for various outputs; think of it as a project tag.'
HELP_1b='
   NAME is used in 2) to retrieve the right objects.'

Help_2='
2) derive resolution-based clusterings from the rcl object.
      rcl.sh -n NAME [-S] -r "N1 N2 N3 .."
      e.g. -r "50 100 200 400 1000 2000"
   A logarithmic scale such as above is suggested.'

HELP_2b='
 In graphs/ in the mcl source distribution, you can run a small test case:
   rcl.sh -n TINY -m rcltiny.mci -t rcltiny.tab rcltiny.cls[123]
   RCL_RES_PLOT_LIMIT=2 rcl.sh -n TINY -r "1 2 3 4"'

HELP_3='
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
      rcl.sh -n NAME -l 470/10/530/pfx'

Help_options='Options:
-m  <file>   Input network/matrix file in mcl format (obtain e.g. with mcxload)
-t  <file>   Tab file with index - label mapping (obtain e.g. with mcxload)
-n  NAME     NAME will be used as prefix for various objects
-l  LOW/STEP/HIGH    e.g. 200/50/700 to show threshold cluster sizes
-r  "N1 N2 N3 .."    e.g. "50 100 200 400" to compute resolution clusterings
-p  <num>    Parallel/CPU, use this many CPUs for parallel RCL compute   
-S           Output cluster size distribution information (use in (2))
-U           Compute the Unrestricted Contingency Linkage object (use in (1))
-F           Force computation, ignore existing RCL object
-H           Expanded help, opened in less'


while getopts :m:n:t:l:r:p:UShFH opt
do
    case "$opt" in
    m)
      matrixfile=$OPTARG
      ;;
    p)
      cpu=$OPTARG
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
    F)
      do_force=true
      ;;
    S)
      do_gralog=true
      ;;
    U)
      do_ucl=true
      ;;
    h)
      cat <<EOU
$Help_1

$Help_2

$Help_options

Use -H to output a longer description
EOU
       exit
      ;;
    H)
   {  cat <<EOU
$HELP_intro

Suggested usage:$Help_1$HELP_1b
$Help_2
$HELP_2b
$HELP_3
$Help_options
EOU
   } | less
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


if $do_force || [[ ! -f $rclfile ]]; then

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
      if (( cpu == 1 )); then
        clm vol --progress -imx $matrixfile -write-rcl $rclfile -o $pfx.vol "$@"
      else
        maxp=$((cpu-1))
        list=$(eval "echo {0..$maxp}")
        echo "-- All $cpu processes are chasing .,=+<>()-"
        for id in $list; do
          clm vol --progress -imx $matrixfile -gi $id/$cpu -write-rcl $pfx.R$id -o pfx.V$id "$@" &
        done
        wait
        clxdo mxsum $(eval echo "$pfx.R{0..$maxp}") > $rclfile
        echo "-- Components summed in $rclfile"
      fi
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
  echo "-- suggest log scale, e.g. rcl.sh -n $pfx -r \"50 100 200 400 1000\""
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
   export RCL_RES_PLOT_LIMIT=${RCL_RES_PLOT_LIMIT-500}

   rcl-res.pl $pfx $RESOLUTION < $pfx.join-order

   echo "-- saving resolution cluster files and"
   echo "-- displaying size of the 20 largest clusters"
   if $do_gralog; then echo "-- summarising cluster size distribution on a log scale";
   else                echo "-- (use in addition to these parameters -S for size distribution summary)"; fi
   echo "-- in parentheses N identical clusters with previous level among 30 largest clusters"
                                       # To help space the granularity output.
   export CLXDO_WIDTH=$((${#pfx}+14))  # .res .cls length 8, leave 6 for resolution

   res_prev=""
   file_prev=""

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
      nshared="--"
      if [[ -n $res_prev ]]; then
         nshared=$(grep -Fcwf <(head -n 30 $rfile | cut -f 1) <(head -n 30 $file_prev | cut -f  1) || true)
      fi
      export CLXDO_GRABIG_TAG="($(printf "%2s" $nshared)) "
      clxdo grabig 20 $prefix.cls
      if $do_gralog; then clxdo gralog $prefix.cls; echo "()"; fi
      res_prev=$r
      file_prev=$rfile
   done
   commalist=$(tr -s ' ' ',' <<< $RESOLUTION)

   if [[ -f $pfx.digraph ]]; then
      if ! dot -Tpng <  $pfx.digraph  > $pfx.cls.png; then
        echo "-- Attempt to visualise hierachy with dot failed"
      else
        echo "-- Hierarchy visualised in $pfx.cls.png with display limit $RCL_RES_PLOT_LIMIT"
      fi
   fi

cat <<EOM

The following outputs were made.
One cluster-per line files with labels:
   $(eval echo $pfx.res{$commalist}.labels)
LABEL<TAB>CLUSID files:
   $(eval echo $pfx.res{$commalist}.txt)
mcl-edge matrix/cluster files (suitable input e.g. for 'clm dist' and others):
   $(eval echo $pfx.res{$commalist}.cls)
EOM

fi



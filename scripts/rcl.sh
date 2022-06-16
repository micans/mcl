#!/bin/bash

# Copyright 2022 Stijn van Dongen

# This program is free software; you can redistribute it and/or modify it
# under the terms of version 3 of the GNU General Public License as published
# by the Free Software Foundation. It should be shipped with MCL in the top
# level directory as the file COPYING.

#    _ _  __|_ _. __|_ _  _|   _ _  _ _|_. _  _  _  _    |. _ |  _  _  _
#   | (/__\ | | |(_ | (/_(_|  (_(_)| | | |(_|(/_| |(_\/  ||| ||<(_|(_|(/_
#                                           |        /               |
# RCL consensus clustering
# Author: Stijn van Dongen
#
# This RCL implementation uses programs/tools that are shipped with mcl.  It
# can be run on any set of clusterings from any method or program, but the
# network and clusterings have to be supplied in mcl matrix format.
#
# See github.com/micans/mcl#rcl and this script with no arguments.

set -euo pipefail

themode=              # first argument, mode 'setup' 'tree' 'select', 'mcl', or 'qc'
projectdir=           # second argument, for modes 'setup' 'tree' 'select'.
network=              # -n FNAME
tabfile=              # -t FNAME
cpu=1                 # -p NUM
RESOLUTION=           # -r, e.g. -r "100 200 400 800 1600 3200"
do_force=false        # -F (for mode 'tree')
mcldir=               # -m DIRNAME to run a bunch of mcl analyses for use as input later
INFLATION=            # -I, e.g. -I "1.3 1.35 1.4 1.45 1.5 1.55 1.6 1.65 1.7 1.8 1.9 2"
SELF=                 # -D
test_ucl=false        # -U

function usage() {
  local e=$1
  cat <<EOH
An rcl workflow requires three steps:

1) rcl.sh setup   TAG -n NETWORKFILENAME -t TABFILENAME  LIST-OF-CLUSTERING-FILES
2) rcl.sh tree    TAG [-F] [-p NCPU]
3) rcl.sh select  TAG -r "RESOLUTIONLIST"

TAG will be used as a project directory name in the current directory.
NETWORKFILENAME and TABFILENAME are usually created by mcxload.
Hard links to these files will be made in TAG, symbolic if this is not possible.
All rcl.sh commands are issued from outside and directly above directory TAG.
LIST-OF-CLUSTERING-FILES is stored in a file and retrieved when needed.
-F forces a run if a previous output exists.
For RESOLUTIONLIST a doubling is suggested, e.g. -r "50 100 200 400 800 1600 3200"
You may want to re-run with a modified list if the largest cluster size is either
too small or too large for your liking.
The history of your commands will be tracked in TAG/rcl.cline.
If 'dot' is available, a plot of results is left in TAG/rcl.hi.RESOLUTION.pdf
A table of all clusters of size above the smallest resolution is left in TAG/rcl.hi.RESOLUTION.txt

To make mcl clusterings to give to rcl:
rcl.sh mcl [-p NCPU] -n NETWORKFILENAME -m OUTPUTDIR -I "INFLATIONLIST"
This may take a while for large graphs.
In step 1) you can then use
   rcl.sh setup TAG -n NETWORKFILENAME -t TABFILENAME  OUTPUTDIR/out.*
INFLATIONLIST:
- for single cell use e.g. -I "1.3 1.35 1.4 1.45 1.5 1.55 1.6 1.65 1.7 1.8 1.9 2"
- for protein families you will probably want somewhat larger values

Additional modes:
rcl.sh qc  TAG      create (1) heat map of clustering discrepancies and (2) granularity plot.
rcl.sh qc2 TAG      create scatter plot of cluster sizes versus induced mean eccentricity of nodes.
EOH
  exit $e
}

MODES="setup tree select mcl qc qc2"

function require_mode() {
  local mode=$1
  [[ -z $mode ]] && echo "Please provide a mode" && usage 0
  if ! grep -qFw -- $mode <<< "$MODES"; then
    echo "Need a mode, one of { $MODES }"; usage 1
  fi
  themode=$mode
}
function require_tag() {
  local mode=$1 tag=$2
  if [[ -z $tag ]]; then
    echo "Mode $mode requires a project directory tag"
    false
  elif [[ $tag =~ ^- ]]; then
    echo "Project directory not allowed to start with a hyphen"
    false
  fi
  projectdir=$tag
}
function require_file() {
  local fname=$1 response=$2
  if [[ ! -f $fname ]]; then
    echo "Expected file $fname not found; $response"
    false
  fi
}
function require_imx() {
  local fname=$1 response=$2
  require_file "$fname" "$response"
  if ! mcx query -imx $fname --dim > /dev/null; then
    echo "Input is not in mcl network format; check file $fname"
    false
  fi
}
function test_exist() {
  local fname=$1
  if [[ -f $fname ]]; then
    if $do_force; then
      echo "Recreating existing file $fname"
    else
      echo "File $fname exists, use -F to force renewal"
      return 1
    fi
  fi
}

require_mode "${1-}"          # themode now set
shift 1
if grep -qFw $themode <<< "setup tree select qc qc2"; then
  require_tag $themode "${1-}"   # projectdir now set
  shift 1
fi

pfx=
if [[ -n $projectdir ]]; then
   mkdir -p $projectdir
   pfx=$projectdir/rcl
   echo -- "$themode $projectdir $@" >> $pfx.cline
fi


while getopts :n:m:p:r:t:I:FDU opt
do
    case "$opt" in
    n) network=$OPTARG ;;
    m) mcldir=$OPTARG ;;
    p) cpu=$OPTARG ;;
    r) RESOLUTION=$OPTARG ;;
    t) tabfile=$OPTARG ;;
    I) INFLATION=$OPTARG ;;
    F) do_force=true ;;
    D) SELF="--self" ;;
    U) test_ucl=true ;;
    :) echo "Flag $OPTARG needs argument" exit 1 ;;
    ?) echo "Flag $OPTARG unknown" exit 1 ;;
   esac
done

shift $((OPTIND-1))

function require_opt() {
  n_req=0
  local mode=$1 option=$2 value=$3 description=$4
  if [[ -z $value ]]; then
    echo "Mode $mode requires $description for option $option"
    (( ++n_req ))
  fi
}


  ##
  ##  M C L

if [[ $themode == 'mcl' ]]; then
  require_opt mcl -n "$network" "a network in mcl format"
  require_opt mcl -m "$mcldir" "a directory output name"
  require_opt mcl -I "$INFLATION" "a set of inflation values between quotes"
  if (( n_req )); then exit 1; fi
  mkdir -p $mcldir
  echo "-- Running mcl for inflation values in ($INFLATION)"
  for I in $(tr -s ' ' '\n' <<< "$INFLATION" | sort -rn); do
    echo -n "-- inflation $I start .."
    mcl $network -I $I -t $cpu --i3 -odir $mcldir
    echo " done"
  done 2> $mcldir/log.mcl
  exit 0


  ##
  ##  S E T U P

elif [[ $themode == 'setup' ]]; then
  require_opt setup -n "$network" "a network in mcl format"
  require_opt setup -t "$tabfile" "a tab file mapping indexes to labels"
  if (( n_req )); then exit 1; fi
  require_imx "$network" "Is $network an mcl matrix file?"

  ndim=$(grep -o "[0-9]\+$" <<< $(mcx query --dim -imx "$network"))
  ntab=$(wc -l < "$tabfile")

  if [[ $ntab != $ndim ]]; then
    echo "Dimension mismatch between network ($ndim) and tab file ($ntab)"
    false
  fi

  (( $# < 2 )) && echo "Please supply a few clusterings" && false
  ls "$@" > $pfx.lsocls
  for f in $(cat $pfx.lsocls); do
    require_imx "$f" "Is clustering $f an mcl matrix file?"
  done
  echo "-- Supplied clusterings are in mcl format"

  if ! ln -f $tabfile $pfx.tab 2> /dev/null; then
    cp $tabfile $pfx.tab
  fi
  if ! ln -f $network $pfx.input 2> /dev/null; then
    ln -sf $network $pfx.input
  fi
  wc -l < $pfx.tab > $pfx.nitems
  echo "Project directory $projectdir is ready ($pfx.input)"


  ##
  ##  T R E E

elif [[ $themode == 'tree' ]]; then
  require_imx "$pfx.input" "did you run rcl.sh setup $projectdir?"
  require_file "$pfx.lsocls" "cluster file $pfx.lsocls is missing, weirdly"

  rclfile=$pfx.rcl
  test_exist "$rclfile"

  mapfile -t cls < $pfx.lsocls
  echo "-- Computing RCL graph on ${cls[@]}"

  if $test_ucl; then
    ucl-simple.sh "$pfx.input" "${cls[@]}"
    mv -f out.ucl $rclfile
    echo "Ran UCL succesfully, output $rclfile was made accordingly"
  elif (( cpu == 1 )); then
    clm vol --progress $SELF -imx $pfx.input -write-rcl $rclfile -o $pfx.vol "${cls[@]}"
  else
    maxp=$((cpu-1))
    list=$(eval "echo {0..$maxp}")
    echo "-- All $cpu processes are chasing .,=+<>()-"
    for id in $list; do
      clm vol --progress $SELF -imx $pfx.input -gi $id/$cpu -write-rcl $pfx.R$id -o pfx.V$id "${cls[@]}" &
    done
    wait
    clxdo mxsum $(eval echo "$pfx.R{0..$maxp}") > $rclfile
    echo "-- Components summed in $rclfile"
  fi
  echo "-- Computing single linkage join order for network $rclfile"
  clm close --sl -sl-rcl-cutoff ${RCL_CUTOFF-0} -imx $rclfile -tab $pfx.tab -o $pfx.join-order -write-sl-list $pfx.node-values
  echo "RCL network and linkage both ready, you can run rcl.sh select $projectdir"


  ##
  ##  R E S

elif [[ $themode == 'select' ]]; then
  minres=-
  require_opt select -r "$RESOLUTION" "a list of resolution values between quotes"
  (( n_req )) && false
  require_imx "$pfx.rcl" "did you run rcl.sh tree $projectdir?"
  echo "-- computing clusterings with resolution parameters $RESOLUTION"
  export MCLXIOVERBOSITY=2
  # export RCL_RES_PLOT_LIMIT=${RCL_RES_PLOT_LIMIT-500}

  rcl-select.pl $pfx $RESOLUTION < $pfx.join-order

  echo "-- saving resolution cluster files and"
  echo "-- displaying size of the 20 largest clusters"
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
    res_prev=$r
    file_prev=$rfile
  done
  commalist=$(tr -s $'\t ' ',' <<< $RESOLUTION)
  hyphenlist=$(tr -s $'\t ' '-' <<< $RESOLUTION)
  resmapfile=$pfx.hi.$hyphenlist.resdot

  if [[ -f $resmapfile ]]; then
    rlist=($RESOLUTION)
    for minres in ${rlist[@]}; do
      resdotfile=$pfx.hi.$minres.dot
      respdffile=$pfx.hi.$minres.pdf
      restxtfile=$pfx.hi.$minres.txt
      rcl-dot-resmap.pl ${RCL_DOT_RESMAP_OPTIONS-} --minres=$minres --label=size < $resmapfile > $resdotfile
      if ! dot -Tpdf -Gsize=10,10\! < $resdotfile > $respdffile; then
        echo "-- dot did not run, pdf not produced"
      else
        echo "-- map of output produced in $respdffile"
      fi
    done
  else
    echo "-- Expected file $resmapfile not present"
  fi

cat <<EOM

The following outputs were made.
One cluster-per line files with labels:
   $(eval echo $pfx.res{$commalist}.labels)
LABEL<TAB>CLUSID files:
   $(eval echo $pfx.res{$commalist}.txt)
mcl-edge matrix/cluster files (suitable input e.g. for 'clm dist' and others):
   $(eval echo $pfx.res{$commalist}.cls)
A table with all clusters at different levels up to size $minres, including their nesting structure:
   $restxtfile
EOM


  ##
  ##  Q C

elif [[ $themode == 'qc' ]]; then

       #       ____  ___)              ______          ______                                 #
      #       (, /   /                (, /    )       (, /    )                              #
     #          /---/   _  __   _       /---(   _       /    / __  _   _   _____   _        #
    #        ) /   (___(/_/ (__(/_   ) / ____)_(/_    _/___ /_/ (_(_(_(_/_(_) / (_/_)_     #
   #        (_/                     (_/ (           (_/___ /         .-/                  #
  #                                                                 (_/                  #
   # Some things still hardcoded, e.g. cluster size range.
   # This looks absolutely awful, with perl madness and inline R scripts.
   # Remember that things like daisies and clover exist.

  require_file "$pfx.lsocls" "(created in step 'tree')"
  require_file "$pfx.nitems" "(created in step 'setup')"

  export RCLPLOT_PARAM_SCALE=${RCLPLOT_PARAM_SCALE:-2}
  mapfile -t cls < $pfx.lsocls

# Heatmap:
  out_heat_txt=$pfx.qc-heat.txt
  out_heat_pdf=${out_heat_txt%.txt}.pdf
  ( echo -e "d\tP1\tP2"; clm dist "${cls[@]}" | rcldo.pl distwrangle $(cat $pfx.nitems) ) > $out_heat_txt # \

  R --slave --quiet --silent --vanilla <<EOR
library(ggplot2, warn.conflicts=FALSE)
library(viridis, warn.conflicts=FALSE)

mytheme = theme(plot.title = element_text(hjust = 0.5), plot.margin=grid::unit(c(5,5,5,5), "mm"),
   axis.text.x=element_text(size=rel(1.0), angle=0), axis.text.y=element_text(size=rel(1.0), angle=0), 
   legend.position="bottom", legend.text=element_text(size=rel(0.6)), legend.title=element_text(size=rel(0.8)),
   text=element_text(family="serif"))
  
h <- read.table("$out_heat_txt",header=T, colClasses=c("double", "character", "character"))
h\$d <- 100 * h\$d

p1 <- ggplot(h, aes(x=P1, y=P2, fill=d)) +  
     geom_tile() + mytheme + ggtitle("Cluster discrepancies") +
     guides(shape = guide_legend(override.aes = list(size = 0.3))) +
     labs(x="Granularity parameter", y="Granularity parameter") +
     geom_text(aes(label=round(d)), color="white", size=4) + scale_fill_viridis(option = "A",
       name = "Distance to GCS as percentage of nodes", direction = -1, limits=c(0,100),
       guide = guide_colourbar(direction = "horizontal",
       draw.ulim = F,
       title.hjust = 0.5,
       barheight = unit(3, "mm"),
       label.hjust = 0.5, title.position = "top"))
ggsave("$out_heat_pdf", width=${RCLPLOT_X:-6}, height=${RCLPLOT_Y:-5}); 
EOR
  echo "-- $out_heat_pdf created"

# Granularity:
  out_gra_txt="$pfx.qc-gra.txt"
  out_gra_pdf="${out_gra_txt%.txt}.pdf"
  for fname in "${cls[@]}"; do
    clxdo gra $fname | tr -s ' ' '\n' | tail -n +2 | rcldo.pl cumgra $fname $(cat $pfx.nitems)
  done > "$out_gra_txt"
  R --slave --quiet --silent --vanilla <<EOR
library(ggplot2, warn.conflicts=FALSE)

mytheme = theme(plot.title = element_text(hjust = 0.5),
  plot.margin=grid::unit(c(4,4,4,4), "mm"), legend.spacing.y=unit(0, 'mm'), legend.key.size = unit(5, "mm"),
               text=element_text(family="serif"))

a <- read.table("$out_gra_txt", colClasses=c("factor", "double", "double"))
xLabels <- c(1,10,100,1000,10000)
  ggplot(a) + geom_line(aes(x = V2, y = V3, group=V1, colour=V1)) +
  mytheme + scale_color_viridis_d() +
  geom_point(aes(x = V2, y = V3, group=V1, colour=V1)) +
  labs(col="${RCLPLOT_GRA_COLOUR:-Granularity parameter}", x="Cluster size x", y=expression("Fraction of nodes in clusters of size "<="x")) +
  scale_x_continuous(breaks = c(0,10,20,30,40),labels= xLabels) +
  ggtitle("${RCLPLOT_GRA_TITLE:-Granularity signatures across inflation}")
ggsave("$out_gra_pdf", width=${RCLPLOT_X:-6}, height=${RCLPLOT_Y:-4})
EOR
  echo "-- $out_gra_pdf created"


elif [[ $themode == 'qc2' ]]; then

  require_file "$pfx.lsocls" "(created in step 'tree')"
  export MCLXIOVERBOSITY=2
  mapfile -t cls < $pfx.lsocls

  # if $do_force || [[ ! -f $pfx.qc2all.txt ]]; then
  if test_exist $pfx.qc2all.txt; then

    for cls in ${cls[@]}; do

      export tag=$(rcldo.pl clstag $cls)
      ecc_file=$pfx.ecc.$tag.txt
      cls_dump=$pfx.clsdump.$tag

      echo "-- computing cluster-wise eccentricity for $cls [$tag]"
      mcx alter -icl $cls -imx $pfx.input --block | mcx diameter -imx - -t $cpu --progress | tail -n +2 > $ecc_file

      mcxdump -imx $cls --transpose --no-values > $cls_dump

      if ! diff -q <(cut -f 1 $ecc_file) <(cut -f 1 $cls_dump); then
        echo "Difference in domains!"
        false
      fi

      paste $ecc_file  <(cut -f 2 $cls_dump) \
        | sort -nk 3 \
        | datamash --group 3 mean 2 count 3 \
        | perl -ne 'chomp; print "$_\t$ENV{tag}\n"' \
        > $pfx.qc2.$tag.txt

    done
    cat $pfx.qc2.*.txt | tac > $pfx.qc2all.txt
  else
    echo "-- reusing $pfx.qc2all.txt"
  fi

  out_qc2_pdf=$pfx.qc2all.pdf

  R --slave --quiet --silent --vanilla <<EOR
library(ggplot2, warn.conflicts=FALSE)
library(viridis, warn.conflicts=FALSE)
d <- read.table("$pfx.qc2all.txt")
d <- d[d\$V3 > 1,]          # remove singletons.
d\$V4 <- as.factor(d\$V4)
mytheme = theme(plot.title = element_text(hjust = 0.5),
  plot.margin=grid::unit(c(4,4,4,4), "mm"), legend.spacing.y=unit(0, 'mm'), legend.key.size = unit(5, "mm"),
               text=element_text(family="serif"))
ggplot() + mytheme +
geom_point(data=d, aes(x=V2, y=log10(V3), colour=V4)) +
expand_limits(x=c(0,${RCL_ECC_X:-10}), y=c(0,${RCL_ECC_Y:-5})) +
labs(colour = "Cls", x="Average eccentricity", y="Cluster size (log10)") +
scale_color_viridis(discrete=TRUE, labels=as.numeric(levels(d\$V4))/10**${RCLPLOT_PARAM_SCALE:-0}) +
ggtitle("${RCLPLOT_ECC_TITLE:-Cluster size / eccentricity}")
ggsave("$out_qc2_pdf", width=${RCLPLOT_X:-5}, height=${RCLPLOT_Y:-5})
EOR
  echo "-- file $out_qc2_pdf created"

fi



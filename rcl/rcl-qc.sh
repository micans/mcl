#!/bin/bash

# Copyright 2022 Stijn van Dongen

# This program is free software; you can redistribute it and/or modify it
# under the terms of version 3 of the GNU General Public License as published
# by the Free Software Foundation. It should be shipped with MCL in the top
# level directory as the file COPYING.

#  _ _  __|_ _. __|_ _  _|   _ _  _ _|_. _  _  _  _    |. _ |  _  _  _   _  _
# | (/__\ | | |(_ | (/_(_|  (_(_)| | | |(_|(/_| |(_\/  ||| ||<(_|(_|(/_ (_|(_
#                                         |        /               |      |/ 
# RCL consensus clustering QC
# Author: Stijn van Dongen
#
# See rcl.sh, github.com/micans/mcl#rcl and this script with no arguments.

set -euo pipefail

themode=              # first argument, mode 'qc' 'qc2' 'heatannot' 'quickmark'
projectdir=           # second argument, for all modes.
infix='infix'         # -x infix, a secondary tag
cpu=1                 # -p NUM
ANNOTATION=           # -a FNAME annotation file
CLUSTERING=           # -c FNAME clustering file (mcl format)
HIERARCHY=            # -h FNAME hierarchy file  (rcl.hm*)
COUNTMATRIX=          # -d FNAME data file
GENENAMES=            # -g FNAME row names for COUNTMATRIX
do_force=false        # -F (for mode 'tree')

function usage() {
  local e=$1
  cat <<EOH
This program provides qc/reporting modes additional to rcl.sh, where TAG is an rcl project
directory. It has some in-line R code to create plots (not pretty, but there are reasons, e.g. file
transformations are needed before the R code can be invoked).  Should you find it useful and wish to
change the R code, you can simply make those changes in (bespoke copies of) this program.  Several
options can be changed from the command line:

  RCLPLOT_YFTSIZE=6 RCLHM_XTITLE="My Title" rcl-qc.sh heatannot TAG -a file1 -c file2

Options that can not be set this way currently are colour schemes, some dimensions as well as many
other things.  To find out what options can be set I'm afraid the best solution is to look inside
the script. Suggestions/patches / wholesale better ideas are welcome.

rcl-qc.sh qc  TAG      create (1) heat map of clustering discrepancies and (2) granularity plot.

rcl-qc.sh qc2 TAG      create scatter plot of cluster sizes versus induced mean eccentricity of nodes.
                       This is compute intensive and benefits from multiple CPUs (-p option).

rcl-qc.sh heatannot TAG -a annotationfile -h rclhierarchyfile (output from rcl select, TAG/rcl.hm*)
rcl-qc.sh heatannot TAG -a annotationfile -c clustering (in mcl matrix format)

                       This mode creates a heatmap from an annotation file, where each
                       node is scored for the same list of traits. Scores are added for each trait and
                       each cluster by adding all scores for that trait for all nodes in the cluster.

                       Annotation file should be rectangular tab-separated, with a header line
                       for the trait/class names and with row names the same as the labels in TAG/rcl.tab.
                       Thus row names represent nodes/cells/barcodes, and A[i,j] is the score for
                       class j for node i.

                       Additional options:
                       -x infix     adds infix to the output files (use to distinguish versions/parameters)

rcl-qc.sh quickmark TAG -d countmatrixfile -g generownamesfile -h rclhierarchyfile
rcl-qc.sh quickmark TAG -d countmatrixfile -g generownamesfile -c clustering

                       Additional options:
                       -x infix     adds infix to the output files (use to distinguish versions/parameters)

R libraries required:
qc and qc2: ggplot2 viridis
heatannot: circlize ComplexHeatmap DECIPHER

MCL sibling programs required:
mcx mcxdump rcldo.pl clxdo
EOH
  exit $e
}

MODES="qc qc2 heatannot quickmark"

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
function test_absence () {
  local fname=$1
  local action=${2:?Programmer error, need second argument}  # return or exit
  if [[ -f $fname ]]; then
    if $do_force; then
      echo "Recreating existing file $fname"
    else
      echo "File $fname exists, use -F to force renewal"
      if [[ $action == 'return' ]]; then
        return 1
      elif [[ $action == 'exit' ]]; then
        exit 1
      fi
    fi
  fi
}

require_mode "${1-}"          # themode now set
shift 1
if grep -qFw $themode <<< "qc qc2 heatannot quickmark"; then
  require_tag $themode "${1-}"   # projectdir now set
  shift 1
fi

pfx=$projectdir/rcl
require_file "$pfx.lsocls" "cluster file $pfx.lsocls is missing"
require_file "$pfx.nitems" "count file $pfx.nitems is missing"
require_file "$pfx.tab" "tab file $pfx.nitems is missing"

echo -- "$themode $projectdir $@" >> $pfx.qcline


while getopts :a:c:d:g:h:p:x:F opt
do
    case "$opt" in
    a) ANNOTATION=$OPTARG ;;
    c) CLUSTERING=$OPTARG ;;
    d) COUNTMATRIX=$OPTARG ;;
    g) GENENAMES=$OPTARG ;;
    h) HIERARCHY=$OPTARG ;;
    p) cpu=$OPTARG ;;
    x) infix=$OPTARG ;;
    F) do_force=true ;;
    :) echo "Flag $OPTARG needs argument"; exit 1 ;;
    ?) echo "Flag $OPTARG unknown"; exit 1 ;;
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
  ##  Q C

if [[ $themode == 'qc' ]]; then

       #       ____  ___)              ______          ______                                 #
      #       (, /   /                (, /    )       (, /    )                              #
     #          /---/   _  __   _       /---(   _       /    / __  _   _   _____   _        #
    #        ) /   (___(/_/ (__(/_   ) / ____)_(/_    _/___ /_/ (_(_(_(_/_(_) / (_/_)_     #
   #        (_/                     (_/ (           (_/___ /         .-/                  #
  #                                                                 (_/                  #
   # Some things still hardcoded, e.g. cluster size range.
   # This looks absolutely awful, with perl madness and inline R scripts.
   # Remember that things like daisies and clover exist.
   # Additionally It Should Work

  require_file "$pfx.lsocls" "(created in step 'tree')"
  require_file "$pfx.nitems" "(created in step 'setup')"

  export RCLPLOT_PARAM_SCALE=${RCLPLOT_PARAM_SCALE:-2}
  mapfile -t cls < $pfx.lsocls

# Heatmap:
  out_heat_txt=$pfx.qc-heat.txt
  out_heat_pdf=${out_heat_txt%.txt}.pdf
  ( echo -e "d\tP1\tP2"; clm dist "${cls[@]}" | rcldo.pl distwrangle $(cat $pfx.nitems) ) > $out_heat_txt # \

  R --no-echo --quiet --silent --vanilla <<EOR
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
    clxdo gra $fname | tr -s ' ' '\n' | tail -n +2 | rcldo.pl granul $fname $(cat $pfx.nitems)
  done > "$out_gra_txt"
  R --no-echo --quiet --silent --vanilla <<EOR
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
  if ! test_absence  $pfx.qc2all.txt return; then
    echo "-- reusing $pfx.qc2all.txt"
  else

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
  fi

  out_qc2_pdf=$pfx.qc2all.pdf

  R --no-echo --quiet --silent --vanilla <<EOR
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


elif [[ $themode == 'heatannot' ]]; then

require_opt $themode -a "$ANNOTATION" "an annotation table with header line, row names as in tab file"
if (( n_req )); then exit 1; fi
require_file "$ANNOTATION" "an annotation table with header line, row names as in tab file"

mybase=$projectdir/hm${infix:+.$infix}

  # cannot use labelcls to unify outputs - rcldo does lots of stuff for hierarchy including Newick tree.
  # 
dispatch_type=none
dispatch_file=none

if [[ -n $CLUSTERING ]]; then     # the old cls mode.
  dispatch_type=heatannotcls
  dispatch_file=$mybase.thecls.txt
  mcxdump -imx $CLUSTERING --dump-rlines --no-values > $dispatch_file
elif [[ -n $HIERARCHY ]]; then
  dispatch_type=heatannot
  dispatch_file=$HIERARCHY
else
  echo "Need -c cls or -h hierarchy argument" && false
fi

how=created
if test_absence $mybase.sum.txt return; then
  # rcldo.pl $dispatch_type $ANNOTATION $dispatch_file $projectdir/rcl.tab > $mybase.txt
  rcldo.pl $dispatch_type $ANNOTATION $dispatch_file $projectdir/rcl.tab $mybase
else
  how=reused
fi
echo "-- file $mybase.sum.txt $how"

# exit 0

  R --no-echo --quiet --silent --vanilla <<EOR
suppressPackageStartupMessages(library(circlize, warn.conflicts=FALSE))
suppressPackageStartupMessages(library(ComplexHeatmap, warn.conflicts=FALSE))
suppressPackageStartupMessages(library(DECIPHER, warn.conflicts=FALSE))           # For newick read
col_mv    = colorRamp2(${RCLHM_MVSCALE:-1} * rev(c(0, 1, 2, 3, 4, 5)), c("darkred", "orange", "lightgoldenrod", "white", "lightblue", "darkblue"))
col_logit = colorRamp2(c(-2, -1, 0, 1, 2, 4, 6)-4, c("darkblue", "lightblue", "white", "lightgoldenrod", "orange", "red", "darkred"))
col_freq  = colorRamp2(${RCLHM_FRSCALE:-1} * c(-5,-3,0,1,2,3,5), c("darkblue", "lightblue", "white", "lightgoldenrod", "orange", "red", "darkred"))
col_zs    = colorRamp2(${RCLHM_ZSSCALE:-1} * c(-2, -1, 0, 1, 2), c("darkblue", "lightblue", "white", "orange", "darkred"))

g  <- read.table("$mybase.sum.txt", header=T, sep="\t", row.names=4)
each_cls_sz  <- g[3:nrow(g),"Size"]
universe_sz <- g\$Size[1]

g2 <- as.matrix(g[3:nrow(g),4:ncol(g)])
term_uv_sum  <- as.numeric(g[1,4:ncol(g)])    # vector over terms, sum score over universe.
term_uv_sdev <- as.numeric(g[2,4:ncol(g)])    # vector over terms, sdev for each term.
term_uv_mean <- term_uv_sum / universe_sz     # vector over terms, mean score over universe.

g3 <- apply(g2, 2, function(x) { x / (each_cls_sz/universe_sz)})
g4 <- apply(g3, 1, function(y) { log(y / term_uv_sum) })
h4 <- t(apply(g2, 2, function(x) { ((x / each_cls_sz)) })) # { -log10(x / (each_cls_sz)) }))

z3 <- apply(g2, 2, "/", each_cls_sz)          # average score for cluster.
z4 <- apply(z3, 1, function(y) { (y - term_uv_mean ) / term_uv_sdev })

# h4: mean value in cluster.
# g4: logarithm of frequency(cluster)/frequency(universe)
# z4: z-score transformed data.

r  <- FALSE             # Either FALSE or a dendrogram ...
if ("$dispatch_type" == 'heatannot') {
  r  <- ReadDendrogram("$mybase.nwk")
  if (nobs(r) != ncol(g4)) {
    stop(sprintf("Dendrogram has %d elements, table $mybase.sum.txt has %d elements", nobs(r), ncol(g4)))
  }
  # used to have/need this: complete dendrogram/matrix ordering understanding still pending.
  # g4 <- g4[,order(order.dendrogram(r))]
}

ls_transform = c("rf", "mv", "zs")
  # rf relative frequency
  # mv mean value
  # zs Z score
for (transform in ls_transform) {
pdf(sprintf("$mybase.%s.pdf", transform), width = ${RCLPLOT_X:-8}, height = ${RCLPLOT_Y:-8})

inf_expected <- ifelse(transform == "rf", "expected", "unexpected")
nna  <- sum(is.na(g4))
ninf <- sum(is.infinite(g4))
if (nna > 0)  { cat(sprintf("<> %d NA resulting for transform %s\n", nna, transform), file=stderr()) }
if (ninf > 0) { cat(sprintf("<> %d infinite instances for transform %s (%s)\n", ninf, transform, inf_expected), file=stderr()) }

g4[is.na(g4)] <- 0
h4[is.na(h4)] <- 0
g4[is.infinite(g4) & g4 > 0] <-  10
g4[is.infinite(g4) & g4 < 0] <- -10
h4[is.infinite(h4) & h4 > 0] <-  10
h4[is.infinite(h4) & h4 < 0] <- -10

myclr <- col_freq
obj <- g4
if (transform == "mv") { myclr <- col_mv; obj <- h4 }
if (transform == "zs") { myclr <- col_zs; obj <- z4 }

  ## the first value is the first residual cluster, usually much larger than the rest.
clr_size = colorRamp2(c(0, median(g\$Size[-1]), max(g\$Size[-c(1,2)])), c("white", "lightgreen", "darkgreen"))
clr_type = colorRamp2(c(0, median(term_uv_sum), max(term_uv_sum)), c("white", "plum1", "purple4"))
size_ha  = HeatmapAnnotation("Cluster size" = each_cls_sz, col=list("Cluster size"=clr_size))
type_ha  = HeatmapAnnotation("Annot" = term_uv_sum, col=list("Annot"=clr_type), which='row')
ht <- Heatmap(obj, name = "${RCLHM_NAME:-Heat}",
  column_title = "${RCLHM_XTITLE:-Clusters}", # row_title = "${RCLHM_YTITLE:-Annotation}",
  cluster_rows = ${RCLHM_ROWCLUSTER:-FALSE},
  cluster_columns= r,
  top_annotation = size_ha,
  right_annotation = type_ha,
  col=myclr,
  row_names_gp = gpar(fontsize = ${RCLPLOT_YFTSIZE:-8}),
  column_names_gp = gpar(fontsize = ${RCLPLOT_XFTSIZE:-8}, fontfamily="Courier"),
  show_column_names = ${RCLHM_CNAMES:-FALSE},
  show_row_names = ${RCLHM_RNAMES:-TRUE},
  row_labels=lapply(rownames(obj), function(x) { substr(x, 1, ${RCLPLOT_YLABELMAX:-20}) }))

options(repr.plot.width = ${RCLPLOT_HM_X:-20}, repr.plot.height = ${RCLPLOT_HM_Y:-16}, repr.plot.res = 100)
ht = draw(ht)
invisible(dev.off())

cat(sprintf("-- relative frequency file $mybase.%s.pdf created\n", transform), file=stderr())
}
EOR


elif [[ $themode == 'quickmark' ]]; then

  if [[ -z ${RCL_SCRIPT_HOME+x} ]]; then
    echo "-- Please set RCL_SCRIPT_HOME to the path where rcl-qm.R is found"
    false
  fi

  mybase=$projectdir/qm${infix:+.$infix}
  mycachebase=$projectdir/qmcache

  echo "-- Using $mybase as file prefix"

  require_opt $themode -d "$COUNTMATRIX" "Count matrix in matrix market format (saved by writeMM)"
  require_opt $themode -g "$GENENAMES" "Gene names - should correpond to row names of count matrix"
  if (( n_req )); then exit 1; fi

  # ? todo check file $mybase.cls
  if [[ -n $CLUSTERING ]]; then     # the old cls mode.
    mcxdump -imx $CLUSTERING -tabr $pfx.tab --dump-rlines --no-values > $mybase.cls
  elif [[ -n $HIERARCHY ]]; then
    rcldo.pl labelcls $pfx.tab < $HIERARCHY > $mybase.cls
  else
    echo "Need -c cls or -h hierarchy argument" && false
  fi

  echo "Starting quickMarkers() analysis [https://github.com/constantAmateur/SoupX/]"
  R --no-echo --silent --quiet --vanilla --args $COUNTMATRIX $GENENAMES $mybase.cls $mybase.annot.txt $mycachebase.qm.txt $mycachebase.data.txt ${RCL_QM_N:-2} ${RCL_QM_TFIDF:-2.0} < $RCL_SCRIPT_HOME/rcl-qm.R
  echo "Output in $mybase.annot.txt"

fi



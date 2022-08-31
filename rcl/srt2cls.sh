#!/bin/bash

# We keep the Seurat 1-based indexing and we want Seurat and mcl indexing to
# refer to the same nodes (so we can easily computed distances etc).  To this
# end we introduce a dummy node with index 0 for mcl; this should already be
# present in the tab file (index-label mapping) given to this script. The first
# line should be '0\tdummy' (the script checks this).
#
# Subsequently we need to add a cluster with this dummy node to the seurat
# clustering when we encode it in an mcl matrix file.

# We do not in all cases maintain cluster label correspondence; one of the
# seurat algorithms labels the clusters starting at 1, this will be shifted to 0.


set -euo pipefail

fn=${1:?Please supply seurat cluster text file}
nodetab=${2:?Please supply mcl tab file for the network nodes}

if (( $(wc -l < $fn) + 1 != $(wc -l < $nodetab) )); then
  echo "Not found: expected different of 1 between file line counts"
  false
fi
if [[ $(head -n 1 $nodetab) != 0$'\t'__dummy__ ]]; then
  echo "Dummy node not found in tab file"
  false
fi

dummyclsid="srtdummyclsid"

(cut -f 2 $fn | sort -un; echo $dummyclsid) | nl -v0 -nln -w1 > .tmp.$fn.clstab

mcxload --transpose --stream-split \
  -abc <(cat $fn; echo -e "dummy\t$dummyclsid") \
  -strict-tabc $nodetab \
  -strict-tabr .tmp.$fn.clstab \
  -o $fn.cls

rm .tmp.$fn.clstab



#!/bin/bash
  
set -euo pipefail

mx=${1?Need matrix file}
dom=${2?Need cluster file}
pfx=${3-diam}

clm close -imx $mx -dom $dom --write-block | mcx diameter --progress | tail -n +2 > $pfx.ecc
mcxdump -imx $dom --no-values --transpose | sort -nk 1 > $pfx.cls

if ! diff -q <(cut -f 1 $pfx.ecc) <(cut -f 1 $pfx.cls); then
   echo "index range not identical"
   false
fi

datamash -g 2 count 4 mean 4 < <(paste $pfx.cls $pfx.ecc | sort -nk 2) | sort -rnk 2


#!/bin/bash

set -euo pipefail

if (( $# == 0 )); then
  cat <<EOU
Usage: srt2tab.sh CELLNAMEFILE > TABFILE
where CELLNAMEFILE was saved from the slot cell.names from the Seurat neighbour object
EOU
  exit 0
fi

( echo __dummy__; cat $1 | tr -d '"') | nl -v0 -nln -w1


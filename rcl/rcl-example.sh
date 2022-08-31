#!/bin/bash
  
set -euo pipefail

TAG=fox        # everything will be made in a directory fox.

if [[ ! -f data/falkner.mci || ! -f data/falkner.tab ]]; then
   echo "Please make a directory 'data' and copy graphs/falkner.mci and graphs/falkner.tab from the mcl source to data"
   exit 0
fi

   # Make mcl matrix IO silent
   #
export MCLXIOVERBOSITY=2

   # This generates many clusterings at different levels of granularity.
   #
echo "1) Generating clusterings in $TAG"; sleep 1
rcl.sh mcl $TAG -n data/falkner.mci -I "1.3 1.4 1.5 1.6 1.8 2.0 2.2 2.5 3.0 4.0"

echo "2) Setting up directory $TAG with network and tab files"; sleep 1
rcl.sh setup $TAG -n data/falkner.mci -t data/falkner.tab

echo "3) Creating RCL matrix and single linkage tree"; sleep 1
   # -F forces overwrite; useful for small example, otherwise not necessary.
rcl.sh tree $TAG -F

echo "4) Creating resolution-based hierarchy"; sleep 1
rcl.sh select $TAG -r "5 10 20 40"



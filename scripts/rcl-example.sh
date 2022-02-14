#!/bin/bash
  
set -euo pipefail

UCL=""
TAG=rcl
if [[ ! -z ${MCLX_USEUCL+x} ]]; then UCL=-U; TAG=ucl; fi

if [[ ! -f falkner.mci || ! -f falkner.tab ]]; then
   echo "Please copy graphs/falkner.mci and graphs/falkner.tab from the mcl source to here"
   exit 0
fi

   # Make mcl matrix IO silent
   #
export MCLXIOVERBOSITY=2

   # This generates many clusterings at different levels of granularity.
   # -dae triggers degree-adjustment with the specified exponent.
   #
echo "1) Generating clusterings"
cached=0
for i in 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0 2.2 2.5 3.0 4.0 5.0 6.0; do
   cline1="mcl falkner.mci -q x -V all -I $i"
   cline2="mcl falkner.mci -q x -V all -I $i -dae  1 -aa .dad"
   ( [[ ! -f $($cline1 -az) ]] && $cline1 ) || (( ++cached ))
   ( [[ ! -f $($cline2 -az) ]] && $cline2 ) || (( ++cached ))
done
echo "Clustering done ($cached cached)"

rcl.sh $TAG $UCL -n falkner.mci -t falkner.tab out.falkner.mci*

rcl.sh $TAG -l 200/100/600

rcl.sh $TAG -r "5 10 20 40"



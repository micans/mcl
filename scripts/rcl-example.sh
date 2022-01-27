#!/bin/bash
  
set -euo pipefail

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
if (( $(ls | grep -c out.falkner.mci) < 3 * 14 )); then
   echo "1) Generating clusterings"
   for i in 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0 2.2 2.5 3.0 4.0 5.0 6.0; do
      mcl falkner.mci -q x -V all -I $i
      mcl falkner.mci -q x -V all -I $i -dae  1 -aa .dad
      mcl falkner.mci -q x -V all -I $i -dae  2 -aa .dae
   done
else
   echo "1) Clusterings already present"
fi


   # Output f.rcl contains the consensus information for nodes co-clustering,
   # weighted with a consistency criterion (restricted contingency).
   # Output f.vol contains the volatility of each node, which equals
   # the maximum scores among its edges in f.rcl (disregarding self-value)
   #
echo "2) computing rcl graph (restricted contingency linkage)"
clm vol -imx falkner.mci -write-rcl f.rcl -o f.vol out.falkner.mci*

   # This computes the minimum spanning tree for f.rcl,
   # equivalent with single link clustering.
   # The join information plus some stats are stored in f.join-order.
   #
echo "3) computing single linkage join order"
clm close --sl -imx f.rcl -tab falkner.tab -o f.join-order -write-sl-list f.node-vals

   # This step can be skipped; it shows some cluster granularities when
   # thresholding the consensus graph at different levels.
   # The graph edge weights fall within [0,1000], where 1000 indicates
   # nodes always cluster together.
   # To find the right dynamic range for the levels it can help to
   # look at a histogram of the edge weights in f.vol, e.g. 
   # the values obtained by
   #                        mcxdump -imx f.vol | cut -f 3 
   # For the resolution parameter in the next stage
   # such inspection is not necessary. The resolution parameter
   # simply indicates the rough range/size of clusters one is
   # interested in.
   #
echo "4) cluster sizes resulting from simple thresholding"
clm close -imx f.rcl -levels 50/100/450/f

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
echo "5) computing balanced clusterings with resolution parameter"
for r in 5 10 20 30; do
                         # this sort orders largest clusters first.
   rcl-mix.pl $r f.join-order | sort -nr > f.mix$r.info
   cut -f 3 f.mix$r.info | mcxload -235-ai - -o f.mix$r.cls
   mcxdump -icl f.mix$r.cls -tabr falkner.tab -o f.mix$r.labels
done


   # Show the cluster sizes of the resolution clusters
   #
echo "7) cluster sizes resulting from resolution-balanced clusterings"
for f in f.mix{5,10,20,30}.cls; do
   printf "%-15s" $f; echo $(mcx query -imx $f | cut -f 2 | tail -n +2)
done



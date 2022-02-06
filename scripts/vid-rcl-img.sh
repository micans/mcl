#!/bin/bash
  
dir=../link

i=${1?Need number}
i3=$(printf "%03d" $i)
prevfile=$(printf "h.%03d.ps" $((2*i-1)))
out1=$(printf "h.%03d.ps" $((2*i   )) )
out2=$(printf "h.%03d.ps" $((2*i+1 )) )

read La Lb Lc Ld Le < <(head -n 1 $dir/out.$i3)

# if node already linked, do not draw the vertex.
# draw newer edges first.

(
  cat $prevfile | grep -v showpage
  echo v${La#xyz} v${Lb#xyz} 0.0 e4
  echo "showpage"
) > $out1

(
   # print previous until and not including %mq
   # perl -ne "print if /^%/ .. /%mq/" $prevfile
   perl -pe "last if /%mq/" $prevfile

   # These are the induced edges (all between linked components)
   tail -n +2 $dir/out.$i3 | while read a b c; do
      echo v${a#xyz} v${b#xyz} $c e2
   done

   echo "%mq"

   # Print vertex if first-time linked.
   if (( Ld == 1 )); then
     echo v${La#xyz} $Lc n4
   fi
   if (( Le == 1 )); then
     echo v${Lb#xyz} $Lc n4
   fi

   # Print MST trunk edges
   echo v${La#xyz} v${Lb#xyz} $Lc e3

   # Print rest of previous entry
   perl -ne "print unless 1 .. /%mq/" $prevfile

) > $out2



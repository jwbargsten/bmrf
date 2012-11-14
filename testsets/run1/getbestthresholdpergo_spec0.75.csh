#!/bin/csh

foreach i (`awk '{print $1}' $1 |sort -u`)
 echo $i `awk '$1==a' a=$i $1 | awk '$7>=0.75' | sort -n -k 8 | tail -1`
end


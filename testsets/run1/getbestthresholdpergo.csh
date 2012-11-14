#!/bin/csh

foreach i (`awk '{print $1}' $1 |sort -u`)
 echo $i `awk '$1==a' a=$i $1 | sort -n -k 9 | tail -1`
end


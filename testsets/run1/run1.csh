#!/bin/csh

set R=/home/adj/software/R-2.10.1/bin/R

cd /home/adj/BMRF/plants/run1

$R --no-save < bmrf_batch.r >& bmrf_batch.out
gzip ytest.out

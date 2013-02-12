#!/usr/bin/env Rscript

setwd("reference_implementation_yiannis_kourmpetis")
source("bmrf.R"); 

net="../test_data/topless_2lvl_int.tsv"
ann="../test_data/go_labels.up.tsv"
ipr="../test_data/ipr_labels.tsv"
out="/tmp/ytest1.out"

set.seed(10)
bmrf(FILEnet=net, FILEann=ann, FILEclust=ipr, outprobsfile=out, no.size.restrictions=TRUE);

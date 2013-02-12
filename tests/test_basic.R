library(bmrf)

FILEnet="test_data/topless_2lvl_int.tsv"
FILEann="test_data/go_labels.up.tsv"
FILEclust="test_data/ipr_labels.tsv"
outprobsfile="/tmp/yntest.out"

bmrf <- read.bmrf(FILEnet, FILEann, FILEclust, minGOsize=6,minFDsize=6, maxGOsize=Inf, maxFDsize=Inf )
set.seed(10); p <- predict.bmrf(bmrf, file=outprobsfile, format="3col", verbose=TRUE)

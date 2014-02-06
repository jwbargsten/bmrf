library(bmrf)

FILEnet="test_data/topless_2lvl_int.tsv"
FILEann="test_data/go_labels.up.tsv"
FILEclust="test_data/ipr_labels.tsv"
outprobsfile="/tmp/yntest.out"

bmrf <- read.bmrf(net.file=FILEnet, go.file=FILEann, fd.file=FILEclust, minGOsize=3,minFDsize=3, maxGOsize=Inf, maxFDsize=Inf )
set.seed(10); p <- predict.bmrf(bmrf, file=outprobsfile, format="3col", verbose=TRUE)

outprobsfile="/tmp/yntest_nofd.out"
bmrf <- read.bmrf(net.file=FILEnet, go.file=FILEann, fd.file=NA, minGOsize=3,minFDsize=3, maxGOsize=Inf, maxFDsize=Inf )
set.seed(10); p <- predict.bmrf(bmrf, file=outprobsfile, format="3col", verbose=TRUE)

outprobsfile="/tmp/yntest.h5"
bmrf <- read.bmrf(net.file=FILEnet, go.file=FILEann, fd.file=FILEclust, minGOsize=3,minFDsize=3, maxGOsize=Inf, maxFDsize=Inf )
set.seed(10); p <- predict.bmrf(bmrf, file=outprobsfile, format="h5", verbose=TRUE)

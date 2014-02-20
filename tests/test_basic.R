library(bmrf)

FILEnet="test_data/topless_2lvl_int.tsv"
FILEann="test_data/go_labels.up.tsv"
FILEclust="test_data/ipr_labels.tsv"
outprobsfile="/tmp/yntest.out"

bmrf <- read.bmrf(net.file=FILEnet, go.file=FILEann, fd.file=FILEclust, minGOsize=3,minFDsize=3, maxGOsize=Inf, maxFDsize=Inf )
set.seed(10); p <- predict.bmrf(bmrf, file=outprobsfile, format="3col", verbose=TRUE)

bmrf <- read.bmrf(net.file=FILEnet, go.file=FILEann, fd.file=FILEclust, minGOsize=3,minFDsize=3, maxGOsize=Inf, maxFDsize=Inf )
bmrf <- transformNet(bmrf, method="degree.meandev")
set.seed(10); p <- predict.bmrf(bmrf, file=outprobsfile, format="3col", verbose=TRUE)

outprobsfile="/tmp/yntest_nofd.out"
bmrf <- read.bmrf(net.file=FILEnet, go.file=FILEann, fd.file=NA, minGOsize=3,minFDsize=3, maxGOsize=Inf, maxFDsize=Inf )
set.seed(10); p <- predict.bmrf(bmrf, file=outprobsfile, format="3col", verbose=TRUE)

outprobsfile="/tmp/yntest.h5"
bmrf <- read.bmrf(net.file=FILEnet, go.file=FILEann, fd.file=FILEclust, minGOsize=3,minFDsize=3, maxGOsize=Inf, maxFDsize=Inf )
set.seed(10); p <- predict.bmrf(bmrf, file=outprobsfile, format="h5", verbose=TRUE)

FILEnet="test_data/yeast_biogrid-string.net.gz"
FILEann="test_data/yeast_biogrid-string.exp.gz"
FILEclust="test_data/yeast_biogrid-string.ipr_dom.gz"
outprobsfile="/tmp/yntest.out"
bmrf <- read.bmrf(net.file=FILEnet, net.header=TRUE, go.file=FILEann, go.header=TRUE, fd.file=FILEclust, fd.header=TRUE, minGOsize=3,minFDsize=3, maxGOsize=Inf, maxFDsize=Inf )

bmrf.t <- clean(transformNet(bmrf, method="degree.component"))

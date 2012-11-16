library(devtools)
load_all("./")
FILEnet="testsets/run1/ExpMatAra.conv.filt2.gz"
FILEann="testsets/run1/extract_exp_tair_BP.lis.uppropagated.filt2.gz"
FILEclust="testsets/run1/TAIR10.ipr.filt2.gz"
outprobsfile="/tmp/ytest.out"

bmrf <- read.bmrf(FILEnet, FILEann, FILEclust)
#p <- predict.bmrf(bmrf)

FILEnet="testsets/run1/ExpMatAra.conv.filt2"
FILEann="testsets/run1/extract_exp_tair_BP.lis.uppropagated.filt2"
FILEclust="testsets/run1/TAIR10.ipr.filt2"
outprobsfile="/tmp/ytest.out"

bmrf <- read.bmrf(FILEnet, FILEann, FILEclust)

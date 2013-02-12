#!/usr/bin/env Rscript


source("bmrf.R"); 

args <- commandArgs(T)

bmrf(FILEnet=args[1], FILEann=args[2], FILEclust=args[3], outprobsfile=args[4]);

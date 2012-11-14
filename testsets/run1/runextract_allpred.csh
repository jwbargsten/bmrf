#!/bin/csh

cd /home/adj/BMRF/plants/run1

#python extract_allpred.py extract_exp_tair_BP.lis.uppropagated.filt2  ytest.out.gz >& extract_allpred.out 
python getthreshold.py extract_allpred.out >& getthreshold.out

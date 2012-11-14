import string
import sys

gofile=sys.argv[1]

def read(file):
 f=open(file,'r')
 data=f.readlines()
 f.close()
 dict={}
 for i in data:
  s=string.split(i)
  prot,go,exp,pred=s
  if dict.has_key(go):
   tmp=dict[go]
  else:
   tmp=[]
  tmp.append([exp,float(pred),prot])
  dict[go]=tmp
 return dict

gdict=read(gofile)
for g in gdict.keys():
 predlist=gdict[g]
 for cun in range(1,100):
  cu=cun/100.0
  tp,fp,tn,fn=0,0,0,0
  for ijk in predlist:
   exp,pred,prot=ijk
   if exp=="0":
    if pred>cu:
     fp=fp+1
   else:
    if pred<cu:
     fn=fn+1
    else:
     tp=tp+1
  small=0.000000000000000000001
  spec=float(tp)/(tp+fp+small)
  cov=float(tp)/(tp+fn+small)
  f=2*spec*cov/(spec+cov+small)
  print "%s %.3f %i %i %i %i %.3f %.3f %.3f"%(g,cu,tp,fp,tn,fn,spec,cov,f)

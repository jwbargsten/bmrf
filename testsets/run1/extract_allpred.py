import string
import sys
import gzip

annotallfile=sys.argv[1]
predfile=sys.argv[2]

def read_plist(file):
 f=open(file,'r')
 data=f.readlines()
 f.close()
 plist=[]
 dict={}
 for i in data:
  s=string.split(i)
  prot,go=s[0],s[1]
  if plist.count(prot)==0:
   plist.append(prot)
   dict[prot]=[go]
  else:
   tmp=dict[prot]
   tmp.append(go)
   dict[prot]=tmp
 return plist,dict

def read_pred(file):
 f=gzip.open(file,'r')
 count=0
 while 1:
  count=count+1
  t=f.readline()
  if not t:
   break
  s=string.split(t)
  prot=s[0]
  go=s[1]
  pred=s[2]
  if allannot.has_key(prot):
   val=allannot[prot].count(go)
   print prot,go,val,pred
 f.close()

allplist,allannot=read_plist(annotallfile)
read_pred(predfile)
